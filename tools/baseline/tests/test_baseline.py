#!/usr/bin/env python3
"""Offline tests; do not access CERN/EOS or submit jobs."""
import array
import json
import os
import shutil
import subprocess
import sys
import tarfile
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

HELPER_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HELPER_DIR))
import baseline_common as common
import baseline_worker as worker
import prepare_baseline as prepare

try:
    import ROOT
    ROOT.gROOT.SetBatch(True)
except ImportError:
    ROOT = None


class WorkflowTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.directory = Path(self.temporary.name)
        self.config = common.load_json(HELPER_DIR / 'baseline_config.json')

    def tearDown(self):
        self.temporary.cleanup()

    def test_bitmap_all_4096_combinations(self):
        accepted = []
        for bits in range(4096):
            passed = ((bits & common.BASELINE_COMMON_MASK) == common.BASELINE_COMMON_MASK
                      and (bits & common.CHANNEL_MASK) != 0)
            if passed:
                accepted.append(bits)
        self.assertEqual(accepted, [3070, 3071, 3582, 3583, 4094, 4095])
        # Signal region is retained whether blinding bit is zero or one.
        self.assertEqual(common.BASELINE_COMMON_MASK & 1, 0)

    def test_discovery_exact_depth_no_cluster_filter(self):
        for relative in ('cluster1.root', 'cluster2.root', '0000/uuid.root',
                         '0000/deeper/excluded.root', 'logs/log.txt'):
            path = self.directory / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.touch()
        settings = dict(input_dir=str(self.directory), input_host='root://unused', maxdepth=1)
        self.assertEqual(len(prepare.discover_files(settings)), 2)
        settings['maxdepth'] = 2
        self.assertEqual(len(prepare.discover_files(settings)), 3)

    def test_failed_remote_listing_not_silently_empty(self):
        with patch.object(prepare.os, 'scandir', side_effect=OSError), \
                patch.object(prepare, 'run', side_effect=RuntimeError('network')):
            with self.assertRaisesRegex(RuntimeError, 'network'):
                prepare.list_entries('root://host', '/missing')

    def test_xrdfs_listing_fallback(self):
        with patch.object(prepare.os, 'scandir', side_effect=OSError), \
                patch.object(prepare, 'run', return_value=(
                    'dr-x 2026-10-07 12:00:00 0 /sample/0000\n'
                    '-r-- 2026-10-07 12:00:00 10 /sample/file.root')):
            self.assertEqual(prepare.list_entries('root://host', '/sample'),
                             [('/sample/0000', True), ('/sample/file.root', False)])

    def test_aliases_and_outputs_are_unique_and_recognizable(self):
        settings = self.config['years']['2022']['ggH']
        path = settings['input_dir'] + '/0000/uuid.root'
        record = common.make_record('ggH', '2022', settings, path)
        self.assertIn('GluGluHtoZG', record['input_name'])
        self.assertIn('Run3Summer22NanoAODv12', record['input_name'])
        self.assertTrue(record['source'].startswith('root://eoscms.cern.ch//eos/cms/'))
        self.assertEqual(record, common.make_record('ggH', '2022', settings, path))
        other = common.make_record('ggH', '2022', settings,
                                   settings['input_dir'] + '/0001/uuid.root')
        self.assertNotEqual(record['output_name'], other['output_name'])
        with self.assertRaises(ValueError):
            common.make_record('ggH', '2022', settings, '/outside/uuid.root')

    def make_manifest(self):
        manifest = dict(year='2022', samples={})
        for sample in ('HJ', 'ggH'):
            settings = self.config['years']['2022'][sample]
            files = [common.make_record(sample, '2022', settings,
                                        settings['input_dir'] + '/{}.root'.format(i))
                     for i in range(3)]
            manifest['samples'][sample] = dict(settings=settings, files=files, process_indices=[1])
        return manifest

    def test_dag_has_normalization_parents_and_one_job_pilot(self):
        manifest = self.make_manifest()
        prepare.write_workflow(self.directory, self.config, manifest, '/afs/user/tmp/proxy')
        dag = (self.directory / 'workflow.dag').read_text()
        self.assertIn('PARENT NORMALIZE_HJ CHILD PROCESS_HJ_1', dag)
        self.assertIn('PARENT NORMALIZE_ggH CHILD PROCESS_ggH_1', dag)
        self.assertNotIn('PROCESS_HJ_0', dag)
        self.assertNotIn('PROCESS_HJ_2', dag)
        for sample in ('HJ', 'ggH'):
            submission = (self.directory / ('process_' + sample + '.sub')).read_text()
            self.assertIn('when_to_transfer_output = ON_SUCCESS', submission)
            self.assertIn('on_exit_hold = ', submission)
            self.assertIn('transfer_output_files = $(pico_name)', submission)
            self.assertIn('/pico/' + sample + '2022pico/', submission)
            self.assertIn('normalization_' + sample + '.json', submission)
            self.assertNotIn('transfer_output_files = raw', submission)
        # Limit concerns conversion only, not the normalization's manifest.
        self.assertEqual(len(manifest['samples']['HJ']['files']), 3)

    def test_package_and_safe_extraction(self):
        checkout = self.directory / 'checkout'
        (checkout / 'data').mkdir(parents=True)
        (checkout / 'data/correction.json').write_text('{}')
        binary = checkout / 'converter'
        binary.write_bytes(b'\x7fELFfixture')
        binary.chmod(0o755)
        bundle = self.directory / 'nano2pico_payload.tar.gz'
        prepare.package_converter(checkout, binary, bundle)
        previous = Path.cwd()
        try:
            os.chdir(self.directory)
            workspace = self.directory / 'work'
            workspace.mkdir()
            worker.extract_payload(workspace)
            self.assertEqual((workspace / 'nano2pico/run/process_nano.exe').read_bytes(), binary.read_bytes())
            with tarfile.open(bundle, 'w:gz') as archive:
                member = tarfile.TarInfo('../escaped')
                member.size = 0
                archive.addfile(member)
            with self.assertRaisesRegex(RuntimeError, 'Unsafe payload'):
                worker.extract_payload(workspace)
            self.assertFalse((self.directory / 'escaped').exists())
        finally:
            os.chdir(previous)

    def test_shell_launcher_rejected(self):
        launcher = self.directory / 'process_nano.exe'
        launcher.write_text('#!/bin/sh\necho launcher\n')
        launcher.chmod(0o755)
        with self.assertRaisesRegex(RuntimeError, 'ELF executable'):
            prepare.resolve_binary(self.directory, str(launcher))

    def test_normalization_keys_match_converter(self):
        sample = self.make_manifest()['samples']['HJ']
        payload = common.normalization_payload(sample, -13.0, [2.0] * 9)
        values = payload[common.NANO_INPUT_DIR][sample['settings']['dataset']]
        self.assertEqual(values['genEventSumw'], -13.0)
        self.assertEqual(len(values), 10)

    def test_prepare_pilot_skips_existing_but_normalizes_entire_set(self):
        helper = self.directory / 'baseline'
        helper.mkdir()
        for name in ('baseline_config.json', 'run_baseline.sh', 'setup_runtime.sh',
                     'baseline_worker.py', 'baseline_common.py'):
            (helper / name).write_bytes((HELPER_DIR / name).read_bytes())
        settings = self.config['years']['2022']['HJ']
        paths = [settings['input_dir'] + '/{}.root'.format(i) for i in range(3)]
        first = common.make_record('HJ', '2022', settings, paths[0])
        with patch.object(prepare, 'HELPER_DIR', helper), \
                patch.object(prepare, 'check_proxy', return_value='/afs/user/tmp/proxy'), \
                patch.object(prepare, 'check_checkout', return_value={}), \
                patch.object(prepare, 'resolve_binary', return_value=HELPER_DIR / 'run_baseline.sh'), \
                patch.object(prepare, 'discover_files', return_value=paths), \
                patch.object(prepare, 'list_entries', return_value=[
                    (settings['output_dir'] + '/' + first['output_name'], False)]), \
                patch.object(prepare, 'run') as commands, \
                patch.object(prepare, 'package_converter'):
            prepare.main(['2022', '--sample', 'HJ', '--campaign', 'pilot',
                          '--nano2pico-dir', str(self.directory),
                          '--limit', '1', '--skip-existing'])
            manifest = common.load_json(helper / 'runs/pilot/manifest.json')
            self.assertEqual(len(manifest['samples']['HJ']['files']), 3)
            self.assertEqual(manifest['samples']['HJ']['process_indices'], [1])
            self.assertFalse((helper / 'runs/pilot/normalization_HJ.json').exists())
            commands.assert_called_once_with([
                'xrdfs', 'root://eosuser.cern.ch', 'mkdir', '-p', settings['output_dir'] + '/logs'])
            # No condor command was executed by preparation.
            self.assertNotIn('condor', commands.call_args[0][0][0])

    def test_source_hash_guard(self):
        source = self.directory / 'src/zgamma_producer.cpp'
        source.parent.mkdir()
        source.write_text('Unreviewed baseline')
        with self.assertRaisesRegex(RuntimeError, 'differs from the reviewed'):
            prepare.check_checkout(self.directory, self.config)

    def test_failure_wrapper_uploads_log_and_preserves_exit_status(self):
        # Only this test runtime is mocked; the real wrapper/pipeline is run.
        shutil.copy2(HELPER_DIR / 'run_baseline.sh', self.directory / 'run_baseline.sh')
        (self.directory / 'setup_runtime.sh').write_text('return 0\n')
        (self.directory / 'baseline_worker.py').write_text(
            'import sys\nprint("Synthetic conversion error")\nsys.exit(7)\n')
        manifest = self.make_manifest()
        (self.directory / 'manifest.json').write_text(json.dumps(manifest))
        fake_proxy = self.directory / 'fake_proxy'
        fake_proxy.write_text('Test placeholder; not a credential')
        bin_dir = self.directory / 'bin'
        bin_dir.mkdir()
        transfer = bin_dir / 'xrdcp'
        transfer.write_text('#!/bin/bash\nprintf "%s\\n" "$@" > "$XRD_CAPTURE"\nexit "${XRD_STATUS:-0}"\n')
        transfer.chmod(0o755)
        environment = dict(os.environ, X509_USER_PROXY=str(fake_proxy),
                           BASELINE_JOB_ID='123.0', XRD_CAPTURE=str(self.directory / 'transfer_args'),
                           PATH=str(bin_dir) + ':' + os.environ['PATH'])
        for transfer_status in ('0', '44'):
            environment['XRD_STATUS'] = transfer_status
            result = subprocess.run(['bash', 'run_baseline.sh', 'process', 'manifest.json',
                                     'HJ', '1', 'normalization_HJ.json'],
                                    cwd=str(self.directory), env=environment,
                                    text=True, capture_output=True)
            self.assertEqual(result.returncode, 7)
            self.assertIn('Synthetic conversion error', result.stdout)
            destination = (self.directory / 'transfer_args').read_text()
            self.assertIn('/pico/HJ2022pico/logs/baseline_failure_HJ_process_1_123.0.log', destination)
            if transfer_status == '44':
                self.assertIn('original exit=7', result.stderr)

    def test_tmp_proxy_rejected_before_any_remote_command(self):
        proxy = self.directory / 'fake_proxy'
        proxy.write_text('Test placeholder; not a credential')
        with patch.dict(os.environ, {'X509_USER_PROXY': str(proxy)}), \
                patch.object(prepare, 'run') as command:
            with self.assertRaisesRegex(RuntimeError, 'Keep the proxy on AFS'):
                prepare.check_proxy()
            command.assert_not_called()


@unittest.skipUnless(ROOT, 'PyROOT unavailable')
class RootTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.directory = Path(self.temporary.name)

    def tearDown(self):
        self.temporary.cleanup()

    def raw_fixture(self, path, bitmaps, nll=1, nphoton=1):
        output = ROOT.TFile(str(path), 'RECREATE')
        tree = ROOT.TTree('tree', 'fixture')
        buffers = {name: array.array('i', [0]) for name in ('nll', 'nphoton', 'zg_cutBitMap')}
        for name, buffer in buffers.items():
            tree.Branch(name, buffer, name + '/I')
        weight, lumi, other = (array.array('d', [value]) for value in (-2.5, -1.25, 42.0))
        for name, buffer in [('weight', weight), ('w_lumi', lumi), ('other_branch', other)]:
            tree.Branch(name, buffer, name + '/D')
        for bits in bitmaps:
            buffers['nll'][0] = nll
            buffers['nphoton'][0] = nphoton
            buffers['zg_cutBitMap'][0] = bits
            tree.Fill()
        tree.Write()
        output.Close()

    def test_skim_every_bitmap_and_preserve_weight_and_other_branches(self):
        raw, selected = self.directory / 'raw.root', self.directory / 'pico.root'
        self.raw_fixture(raw, range(4096))
        self.assertEqual(worker.skim_baseline(ROOT, raw, selected, {'sample': 'test'}), (4096, 6))
        output = ROOT.TFile.Open(str(selected))
        tree = output.Get('tree')
        self.assertEqual([int(event.zg_cutBitMap) for event in tree], [3070, 3071, 3582, 3583, 4094, 4095])
        for event in tree:
            self.assertEqual(event.weight, -2.5)
            self.assertEqual(event.other_branch, 42.0)
        self.assertEqual(json.loads(output.Get('baseline_provenance').GetTitle())['baseline_events'], 6)
        output.Close()

    def test_zero_passing_events_is_valid(self):
        raw, selected = self.directory / 'raw.root', self.directory / 'pico.root'
        self.raw_fixture(raw, [0, 2558, 3582], nphoton=0)
        self.assertEqual(worker.skim_baseline(ROOT, raw, selected, {}), (3, 0))
        output = ROOT.TFile.Open(str(selected))
        self.assertTrue(output.Get('tree').GetBranch('other_branch'))
        output.Close()

    def test_missing_baseline_branch_fails(self):
        raw = self.directory / 'raw.root'
        output = ROOT.TFile(str(raw), 'RECREATE')
        ROOT.TTree('tree', 'missing branches').Write()
        output.Close()
        with self.assertRaisesRegex(RuntimeError, 'Missing pico branch'):
            worker.skim_baseline(ROOT, raw, self.directory / 'selected.root', {})

    def test_full_sample_multi_run_signed_normalization(self):
        sample = dict(settings={'dataset': 'GluGluHtoZG'}, files=[])
        for file_index, runs in enumerate(([10.0, -2.0], [5.0])):
            filename = self.directory / ('nano{}.root'.format(file_index))
            output = ROOT.TFile(str(filename), 'RECREATE')
            events = ROOT.TTree('Events', 'events')
            value = array.array('i', [1])
            events.Branch('event', value, 'event/I')
            for _ in range(3):
                events.Fill()
            events.Write()
            metadata = ROOT.TTree('Runs', 'runs')
            sumw = array.array('d', [0])
            scales = array.array('d', [1.0] * 9)
            metadata.Branch('genEventSumw', sumw, 'genEventSumw/D')
            metadata.Branch('LHEScaleSumw', scales, 'LHEScaleSumw[9]/D')
            for weight in runs:
                sumw[0] = weight
                metadata.Fill()
            metadata.Write()
            output.Close()
            sample['files'].append(dict(source=str(filename)))
        norm_path = self.directory / 'norm.json'
        worker.normalize(ROOT, sample, norm_path)
        norm = common.load_json(norm_path)[common.NANO_INPUT_DIR]['GluGluHtoZG']
        self.assertEqual(norm['genEventSumw'], 13.0)
        self.assertEqual(norm['LHEScaleSumw0'], 3.0)

    def test_process_sequence_and_only_final_output(self):
        config = common.load_json(HELPER_DIR / 'baseline_config.json')
        settings = config['years']['2022']['HJ']
        record = common.make_record('HJ', '2022', settings, settings['input_dir'] + '/input.root')
        checkout = self.directory / 'checkout'
        (checkout / 'data').mkdir(parents=True)
        binary = checkout / 'converter'
        binary.write_bytes(b'\x7fELFfixture')
        binary.chmod(0o755)
        prepare.package_converter(checkout, binary, self.directory / 'nano2pico_payload.tar.gz')
        manifest = dict(year='2022', nano2pico={'binary_sha256': common.sha256_file(binary)},
                        samples={'HJ': dict(settings=settings, files=[record], process_indices=[0])})
        norm_path = self.directory / 'normalization_HJ.json'
        norm_path.write_text(json.dumps(common.normalization_payload(manifest['samples']['HJ'], 10, [1]*9)))
        calls = []

        def fake_command(command, check, cwd=None):
            calls.append((command, cwd))
            if command[0] == 'xrdcp':
                output = ROOT.TFile(command[-1], 'RECREATE')
                events = ROOT.TTree('Events', 'events')
                event = array.array('i', [0])
                events.Branch('event', event, 'event/I')
                for _ in range(4):
                    events.Fill()
                events.Write()
                output.Close()
            else:
                expected_name = command[command.index('-f') + 1]
                self.raw_fixture(Path(cwd) / 'out/zgamma/raw_pico' / ('raw_pico_' + expected_name),
                                 [0, 3070, 3582, 3583])

        previous = Path.cwd()
        try:
            os.chdir(self.directory)
            with patch.object(worker.subprocess, 'run', side_effect=fake_command):
                worker.process(ROOT, manifest, 'HJ', 0, norm_path)
            self.assertTrue((self.directory / record['output_name']).is_file())
            converter_command = calls[1][0]
            self.assertEqual(converter_command[converter_command.index('-i') + 1], common.NANO_INPUT_DIR)
            self.assertNotIn('--skim', converter_command)
            self.assertNotIn('--nent', converter_command)
            self.assertIn('out/zgamma', converter_command)
            self.assertEqual(len(list(self.directory.glob('*.root'))), 1)
            # No direct EOS write: HTCondor uploads the explicitly listed pico.
            self.assertEqual(len(calls), 2)
        finally:
            os.chdir(previous)


if __name__ == '__main__':
    unittest.main(verbosity=2)
