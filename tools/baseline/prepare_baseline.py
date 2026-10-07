#!/usr/bin/env python3
"""Prepare, but never submit, a CERN HTCondor NanoAOD -> baseline pico DAG."""
import argparse
import json
import os
import shutil
import subprocess
import tarfile
from pathlib import Path

from baseline_common import SAFE_NAME, load_json, make_record, sha256_file

HELPER_DIR = Path(__file__).resolve().parent


def run(command):
    return subprocess.check_output(command, text=True).strip()


def list_entries(host, directory):
    """Return (absolute path, is_dir), using a mount or checked XRootD listing."""
    try:
        with os.scandir(directory) as entries:
            return sorted((entry.path, entry.is_dir(follow_symlinks=False))
                          for entry in entries)
    except OSError:
        output = run(['xrdfs', host, 'ls', '-l', directory])
        entries = []
        for line in output.splitlines():
            columns = line.split()
            if len(columns) < 2 or columns[0][0] not in '-d':
                raise RuntimeError('Unrecognized xrdfs listing: ' + line)
            path = columns[-1]
            if not path.startswith(directory.rstrip('/') + '/'):
                raise RuntimeError('Unexpected path in xrdfs listing: ' + path)
            entries.append((path, columns[0].startswith('d')))
        return sorted(entries)


def discover_files(settings):
    pending = [(settings['input_dir'], 1)]
    found = []
    while pending:
        directory, depth = pending.pop()
        for path, is_dir in list_entries(settings['input_host'], directory):
            if is_dir and depth < settings['maxdepth']:
                pending.append((path, depth + 1))
            elif not is_dir and path.endswith('.root'):
                found.append(path)
    return sorted(set(found))


def check_proxy():
    value = os.environ.get('X509_USER_PROXY', '')
    if not value:
        raise RuntimeError('Export X509_USER_PROXY to a valid proxy on AFS first.')
    proxy = Path(value).resolve(strict=True)
    # /tmp belongs to one login host and cannot be read by the remote schedd.
    if not str(proxy).startswith('/afs/'):
        raise RuntimeError('Keep the proxy on AFS, not /tmp: ' + str(proxy))
    if not os.access(str(proxy), os.R_OK):
        raise RuntimeError('Proxy is not readable: ' + str(proxy))
    remaining = int(run(['voms-proxy-info', '--file', str(proxy), '--timeleft']))
    if remaining <= 0:
        raise RuntimeError('Proxy has expired; renew it before preparing/submitting.')
    if remaining < 86400:
        print('WARNING: proxy expires in less than 24 hours; renew before submission.')
    return str(proxy)


def resolve_binary(checkout, explicit=None):
    if explicit:
        binary = Path(explicit).expanduser().resolve(strict=True)
    else:
        candidates = list(checkout.glob('kernel/*/run/process_nano.exe'))
        direct = checkout / 'run/process_nano.exe'
        if direct.is_file() and direct.read_bytes()[:4] == b'\x7fELF':
            candidates.append(direct)
        if len(candidates) != 1:
            raise RuntimeError('Build nano2pico first, or choose the real ELF binary '
                               'with --binary. Found: ' + str(candidates))
        binary = candidates[0].resolve()
    with binary.open('rb') as stream:
        if stream.read(4) != b'\x7fELF':
            raise RuntimeError('Expected an ELF executable, not run/ shell launcher: '
                               + str(binary))
    if not os.access(str(binary), os.X_OK):
        raise RuntimeError('Binary is not executable: ' + str(binary))
    newer = [path for folder in ('src', 'inc')
             for path in (checkout / folder).rglob('*')
             if path.is_file() and path.suffix in ('.cpp', '.cxx', '.hpp')
             and path.stat().st_mtime > binary.stat().st_mtime]
    if newer:
        raise RuntimeError('Source is newer than the binary; rebuild nano2pico. '
                           'Example: ' + str(newer[0]))
    return binary


def check_checkout(checkout, config):
    source = checkout / 'src/zgamma_producer.cpp'
    if sha256_file(source) != config['zgamma_source_sha256']:
        raise RuntimeError('zgamma_producer.cpp differs from the reviewed '
                           + config['upstream_ref'] + ' implementation. '
                           'Review baseline bit meanings before updating the source hash.')
    if not (checkout / 'data/zgamma/2022').is_dir():
        raise RuntimeError('Missing nano2pico 2022 correction data')
    return dict(commit=run(['git', '-C', str(checkout), 'rev-parse', 'HEAD']),
                dirty=bool(run(['git', '-C', str(checkout), 'status', '--porcelain'])),
                zgamma_source_sha256=sha256_file(source),
                process_source_sha256=sha256_file(checkout / 'src/process_nano.cxx'))


def package_converter(checkout, binary, destination):
    with tarfile.open(destination, 'w:gz', dereference=True) as bundle:
        bundle.add(binary, arcname='nano2pico/run/process_nano.exe')
        bundle.add(checkout / 'data', arcname='nano2pico/data')


def quote_submit(value):
    value = str(value)
    if any(character in value for character in '\n\r"'):
        raise ValueError('Unsupported character in submit value')
    return '"' + value + '"'


def base_submit(config, campaign, proxy):
    return '\n'.join([
        'universe = vanilla',
        'notification = Never',
        'initialdir = ' + quote_submit(campaign),
        'executable = ' + quote_submit(campaign / 'run_baseline.sh'),
        'transfer_executable = True',
        'should_transfer_files = YES',
        'when_to_transfer_output = ON_SUCCESS',
        'success_exit_code = 0',
        'on_exit_hold = (ExitBySignal == True) || (ExitCode != 0)',
        'request_cpus = ' + str(config['request_cpus']),
        'request_memory = ' + str(config['request_memory_mb']) + 'MB',
        'request_disk = ' + str(config['request_disk_mb']) + 'MB',
        '+JobFlavour = ' + quote_submit(config['job_flavour']),
        'requirements = HasSingularity',
        '+SingularityImage = ' + quote_submit(config['container_image']),
        'use_x509userproxy = True',
        'x509userproxy = ' + quote_submit(proxy),
        'environment = "BASELINE_JOB_ID=$(Cluster).$(Process)"',
    ]) + '\n'


def write_workflow(campaign, config, manifest, proxy):
    common = base_submit(config, campaign, proxy)
    helper_files = ['manifest.json', 'baseline_config.json', 'setup_runtime.sh',
                    'baseline_worker.py', 'baseline_common.py']
    inputs = ', '.join(str(campaign / name) for name in helper_files)
    dag = ['# No automatic retries: inspect failures before resubmitting.', '']
    for sample, info in manifest['samples'].items():
        norm_name = 'normalization_{}.json'.format(sample)
        norm_sub = common + '\n'.join([
            'arguments = normalize manifest.json ' + sample + ' ' + norm_name,
            'transfer_input_files = ' + inputs,
            'transfer_output_files = ' + norm_name,
            'log = logs/normalize_{}.$(Cluster).$(Process).log'.format(sample),
            'output = logs/normalize_{}.$(Cluster).$(Process).out'.format(sample),
            'error = logs/normalize_{}.$(Cluster).$(Process).err'.format(sample),
            'queue 1', '',
        ])
        (campaign / ('normalize_' + sample + '.sub')).write_text(norm_sub)
        dag.append('JOB NORMALIZE_{0} normalize_{0}.sub'.format(sample))
        process_sub = common + '\n'.join([
            'arguments = process manifest.json ' + sample + ' $(file_index) ' + norm_name,
            'transfer_input_files = ' + inputs + ', ' + str(campaign / 'nano2pico_payload.tar.gz')
            + ', ' + str(campaign / norm_name),
            'transfer_output_files = $(pico_name)',
            'output_destination = root://eosuser.cern.ch/' + info['settings']['output_dir'] + '/',
            'MY.XRDCP_CREATE_DIR = True',
            # Scheduler event logs stay on AFS. stdout/stderr are transferred to
            # the EOS output destination by the CERN XRootD output plugin.
            'log = logs/process_{}.$(file_index).$(Cluster).$(Process).log'.format(sample),
            'output = logs/process_{}.$(file_index).$(Cluster).$(Process).out'.format(sample),
            'error = logs/process_{}.$(file_index).$(Cluster).$(Process).err'.format(sample),
            'queue 1', '',
        ])
        (campaign / ('process_' + sample + '.sub')).write_text(process_sub)
        for index in info['process_indices']:
            node = 'PROCESS_{}_{}'.format(sample, index)
            dag += [
                'JOB {} process_{}.sub'.format(node, sample),
                'VARS {} file_index="{}" pico_name="{}"'.format(
                    node, index, info['files'][index]['output_name']),
                'PARENT NORMALIZE_{} CHILD {}'.format(sample, node),
            ]
        dag.append('')
    (campaign / 'workflow.dag').write_text('\n'.join(dag) + '\n')


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('year', choices=['2022'])
    parser.add_argument('--sample', choices=['HJ', 'ggH', 'both'], default='both')
    parser.add_argument('--campaign', required=True, help='Unique label, e.g. pilot_2022')
    parser.add_argument('--limit', type=int, help='Process only this many files per sample; '
                        'normalization still uses ALL input files.')
    parser.add_argument('--skip-existing', action='store_true',
                        help='Skip existing EOS picos; still normalize the full input set.')
    parser.add_argument('--nano2pico-dir', help='Override the configured existing checkout')
    parser.add_argument('--binary', help='Explicit path to compiled ELF converter')
    args = parser.parse_args(argv)
    if not SAFE_NAME.fullmatch(args.campaign) or args.campaign in ('.', '..'):
        parser.error('Campaign must contain only letters, digits, _, -, or .')
    if args.limit is not None and args.limit < 1:
        parser.error('--limit must be positive')
    config = load_json(HELPER_DIR / 'baseline_config.json')
    checkout = Path(args.nano2pico_dir or config['nano2pico_dir']).expanduser().resolve(strict=True)
    campaign = HELPER_DIR / 'runs' / args.campaign
    if campaign.exists():
        parser.error('Campaign already exists; use a new label or resume its DAG: ' + str(campaign))
    proxy = check_proxy()
    provenance = check_checkout(checkout, config)
    binary = resolve_binary(checkout, args.binary)
    provenance['binary_sha256'] = sha256_file(binary)
    samples = ['HJ', 'ggH'] if args.sample == 'both' else [args.sample]
    manifest = dict(year=args.year, reviewed_commit=config['reviewed_commit'],
                    nano2pico=provenance, samples={})
    for sample in samples:
        settings = config['years'][args.year][sample]
        paths = discover_files(settings)
        if not paths:
            raise RuntimeError('No ROOT inputs found for ' + sample)
        files = [make_record(sample, args.year, settings, path) for path in paths]
        names = [record['output_name'] for record in files]
        if len(set(names)) != len(names):
            raise RuntimeError('Output name collision')
        # These are the user-requested output directories only, never the inputs.
        run(['xrdfs', 'root://eosuser.cern.ch', 'mkdir', '-p', settings['output_dir'] + '/logs'])
        existing = {Path(path).name for path, is_dir in
                    list_entries('root://eosuser.cern.ch', settings['output_dir']) if not is_dir}
        indices = [i for i, record in enumerate(files) if record['output_name'] not in existing]
        conflicts = len(files) - len(indices)
        if conflicts and not args.skip_existing:
            raise RuntimeError('{} existing {} picos: use --skip-existing; no overwrite is allowed.'
                               .format(conflicts, sample))
        if args.limit is not None:
            indices = indices[:args.limit]
        if not indices:
            print(sample + ': all files already have outputs; no jobs prepared.')
            continue
        manifest['samples'][sample] = dict(settings=settings, files=files, process_indices=indices)
        print('{}: {} total inputs, {} conversion/baseline jobs'.format(sample, len(files), len(indices)))
    if not manifest['samples']:
        return 0
    campaign.mkdir(parents=True)
    (campaign / 'logs').mkdir()
    for name in ('run_baseline.sh', 'setup_runtime.sh', 'baseline_config.json',
                 'baseline_worker.py', 'baseline_common.py'):
        shutil.copy2(HELPER_DIR / name, campaign / name)
    (campaign / 'run_baseline.sh').chmod(0o755)
    (campaign / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    package_converter(checkout, binary, campaign / 'nano2pico_payload.tar.gz')
    write_workflow(campaign, config, manifest, proxy)
    print('\nPrepared, NOT submitted: ' + str(campaign / 'workflow.dag'))
    print('cd ' + str(campaign))
    print('condor_submit_dag -no_submit workflow.dag  # DAG generation check only')
    print('condor_submit_dag -maxjobs 100 workflow.dag')
    return 0


if __name__ == '__main__':
    try:
        raise SystemExit(main())
    except (OSError, ValueError, RuntimeError, subprocess.CalledProcessError) as error:
        raise SystemExit('ERROR: ' + str(error))
