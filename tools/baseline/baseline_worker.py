#!/usr/bin/env python3
"""Condor normalization and per-file conversion/baseline tasks."""
import argparse
import json
import math
import subprocess
import tarfile
from pathlib import Path

from baseline_common import (BASELINE_CUT, NANO_INPUT_DIR, load_json,
                             normalization_payload, sha256_file)


def open_root(ROOT, filename):
    root_file = ROOT.TFile.Open(str(filename), 'READ')
    if not root_file or root_file.IsZombie() or root_file.TestBit(ROOT.TFile.kRecovered):
        raise RuntimeError('Invalid/incomplete ROOT file: ' + str(filename))
    return root_file


def normalize(ROOT, sample, output):
    """Same sums/schema as upstream find_normalization.py, across our manifest.

    Unlike upstream, missing/invalid metadata is fatal instead of producing
    zero normalizations. Never normalize only the events passing the baseline.
    """
    sumw = 0.0
    scales = [0.0] * 9
    entries = 0
    for index, record in enumerate(sample['files']):
        root_file = open_root(ROOT, record['source'])
        try:
            runs = root_file.Get('Runs')
            events = root_file.Get('Events')
            if not runs or runs.GetEntries() < 1 or not events:
                raise RuntimeError('Missing Runs/Events: ' + record['source'])
            for branch in ('genEventSumw', 'LHEScaleSumw'):
                if not runs.GetBranch(branch):
                    raise RuntimeError('Missing ' + branch + ': ' + record['source'])
            entries += int(events.GetEntries())
            for run in runs:
                weight = float(run.genEventSumw)
                variations = list(run.LHEScaleSumw)
                if not math.isfinite(weight) or len(variations) != 9:
                    raise RuntimeError('Invalid sumw / non-nine LHE scale metadata: ' + record['source'])
                sumw += weight
                for i, value in enumerate(variations):
                    value = float(value)
                    if not math.isfinite(value):
                        raise RuntimeError('Nonfinite LHE scale sum: ' + record['source'])
                    scales[i] += value
        finally:
            root_file.Close()
        if (index + 1) % 100 == 0 or index + 1 == len(sample['files']):
            print('Normalization: {}/{} files'.format(index + 1, len(sample['files'])), flush=True)
    if not math.isfinite(sumw) or sumw == 0 or any(not math.isfinite(x) or x == 0 for x in scales):
        raise RuntimeError('Zero or nonfinite global normalization')
    Path(output).write_text(json.dumps(normalization_payload(sample, sumw, scales), indent=2) + '\n')
    print('Full sample: {} events, sum(genWeight)={:.12g}'.format(entries, sumw))


def extract_payload(destination):
    with tarfile.open('nano2pico_payload.tar.gz', 'r:gz') as bundle:
        for member in bundle.getmembers():
            target = (destination / member.name).resolve()
            if destination.resolve() not in target.parents or not (member.isfile() or member.isdir()):
                raise RuntimeError('Unsafe payload member: ' + member.name)
        if hasattr(tarfile, 'data_filter'):
            bundle.extractall(destination, filter='data')
        else:
            # CMSSW's Python 3.9 predates extraction filters. The checks above
            # still reject traversal, links, special devices, and FIFOs.
            bundle.extractall(destination)


def skim_baseline(ROOT, raw_path, output_path, provenance):
    """Copy all pico branches, selecting the upstream MC baseline bitmap."""
    if Path(output_path).exists():
        raise RuntimeError('Refusing to overwrite ' + str(output_path))
    source = open_root(ROOT, raw_path)
    output = None
    try:
        tree = source.Get('tree')
        if not tree:
            raise RuntimeError('Missing pico tree')
        for branch in ('nll', 'nphoton', 'zg_cutBitMap', 'weight', 'w_lumi'):
            if not tree.GetBranch(branch):
                raise RuntimeError('Missing pico branch: ' + branch)
        raw_entries = int(tree.GetEntries())
        expected = int(tree.GetEntries(BASELINE_CUT))
        output = ROOT.TFile.Open(str(output_path), 'RECREATE')
        if not output or output.IsZombie():
            raise RuntimeError('Cannot create output pico')
        output.cd()
        selected = tree.CopyTree(BASELINE_CUT)
        if not selected or int(selected.GetEntries()) != expected:
            raise RuntimeError('Baseline tree copy failed')
        if selected.Write('tree', ROOT.TObject.kOverwrite) <= 0:
            raise RuntimeError('Pico write failed')
        provenance = dict(provenance, baseline_cut=BASELINE_CUT,
                          raw_events=raw_entries, baseline_events=expected)
        ROOT.TNamed('baseline_provenance', json.dumps(provenance, sort_keys=True)).Write()
        output.Close()
        output = None
    finally:
        if output:
            output.Close()
        source.Close()
    check = open_root(ROOT, output_path)
    try:
        saved = check.Get('tree')
        if not saved or int(saved.GetEntries()) != expected:
            raise RuntimeError('Saved pico validation failed')
        if int(saved.GetEntries(BASELINE_CUT)) != expected:
            raise RuntimeError('Saved pico contains events failing the baseline')
    finally:
        check.Close()
    print('Baseline: {} / {} raw pico events -> {}'.format(expected, raw_entries, output_path))
    return raw_entries, expected


def process(ROOT, manifest, sample_name, index, norm_file):
    sample = manifest['samples'][sample_name]
    record = sample['files'][index]
    if index not in sample['process_indices']:
        raise RuntimeError('File is not scheduled in this campaign')
    norm = load_json(norm_file)
    if sample['settings']['dataset'] not in norm.get(NANO_INPUT_DIR, {}):
        raise RuntimeError('Normalization directory/tag does not match process_nano')
    if Path(record['output_name']).exists():
        raise RuntimeError('Output already exists in scratch')
    workspace = Path('work').resolve()
    workspace.mkdir(exist_ok=False)
    extract_payload(workspace)
    repo = workspace / 'nano2pico'
    executable = repo / 'run/process_nano.exe'
    if sha256_file(executable) != manifest['nano2pico']['binary_sha256']:
        raise RuntimeError('Converter checksum does not match the manifest')
    input_dir = repo / NANO_INPUT_DIR
    input_dir.mkdir(parents=True)
    local_input = input_dir / record['input_name']
    subprocess.run(['xrdcp', '--nopbar', record['source'], str(local_input)], check=True)
    nano = open_root(ROOT, local_input)
    try:
        events = nano.Get('Events')
        if not events:
            raise RuntimeError('Missing NanoAOD Events tree')
        nano_entries = int(events.GetEntries())
    finally:
        nano.Close()
    output_dir = repo / 'out/zgamma'
    (output_dir / 'raw_pico').mkdir(parents=True)
    # CWD supplies relative correction-data paths. -i and JSON keys are identical.
    # No --nent or --skim: process every input event before baseline filtering.
    subprocess.run([str(executable), '-i', NANO_INPUT_DIR, '-f', record['input_name'],
                    '-o', 'out/zgamma', '--norm', str(Path(norm_file).resolve())],
                   cwd=str(repo), check=True)
    raw_path = output_dir / 'raw_pico' / ('raw_pico_' + record['input_name'])
    provenance = dict(year=manifest['year'], sample=sample_name, source=record['source'],
                      nano2pico=manifest['nano2pico'],
                      normalization_sha256=sha256_file(norm_file))
    raw_count, selected = skim_baseline(ROOT, raw_path, record['output_name'], provenance)
    if raw_count != nano_entries:
        raise RuntimeError('Unexpected pre-baseline event loss: Nano={} raw={}'.format(nano_entries, raw_count))
    print('Ready for EOS transfer: {} ({} baseline events)'.format(record['output_name'], selected))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=['normalize', 'process'])
    parser.add_argument('manifest')
    parser.add_argument('sample', choices=['HJ', 'ggH'])
    parser.add_argument('task', help='Normalization output filename or input-file index')
    parser.add_argument('norm_file', nargs='?')
    args = parser.parse_args(argv)
    manifest = load_json(args.manifest)
    if args.sample not in manifest['samples']:
        parser.error('Sample is not in this manifest')
    import ROOT
    ROOT.gROOT.SetBatch(True)
    if args.mode == 'normalize':
        if args.norm_file is not None:
            parser.error('normalize accepts one output filename')
        normalize(ROOT, manifest['samples'][args.sample], args.task)
    else:
        if args.norm_file is None:
            parser.error('process needs a normalization JSON filename')
        process(ROOT, manifest, args.sample, int(args.task), args.norm_file)


if __name__ == '__main__':
    try:
        main()
    except (OSError, ValueError, RuntimeError, subprocess.CalledProcessError) as error:
        raise SystemExit('ERROR: ' + str(error))
