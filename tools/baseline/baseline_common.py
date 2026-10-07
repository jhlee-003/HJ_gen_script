"""Shared, ROOT-independent configuration and baseline definition."""
import hashlib
import json
import re
from pathlib import Path

# src/zgamma_producer.cpp: bits 11, 8, 7, 6, 5, 4, 3, 2, 1.
# Bits 10/9 are the alternative ee/mumu channels. Bit 0 blinds data; IGNORE it.
BASELINE_COMMON_MASK = 0b100111111110
CHANNEL_MASK = 0b011000000000
BASELINE_CUT = (
    'nll >= 1 && nphoton >= 1 && '
    '(zg_cutBitMap & {0}) == {0} && '
    '(zg_cutBitMap & {1}) != 0'
).format(BASELINE_COMMON_MASK, CHANNEL_MASK)
NANO_INPUT_DIR = 'NanoAODv12/nano/2022/mc'
SAFE_NAME = re.compile(r'^[A-Za-z0-9_.-]+$')


def load_json(path):
    with open(path) as stream:
        return json.load(stream)


def sha256_file(path):
    digest = hashlib.sha256()
    with open(path, 'rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def make_record(sample, year, settings, path):
    path = str(Path(path))
    if not path.startswith(settings['input_dir'].rstrip('/') + '/'):
        raise ValueError('Input is outside the configured sample directory: ' + path)
    url = settings['input_host'].rstrip('/') + '/' + path
    identity = hashlib.sha256(url.encode()).hexdigest()[:20]
    # Central EOS files are UUID.root. These aliases are REQUIRED for nano2pico
    # to recognize the MC process, year, and NanoAOD version. No source renames.
    local_name = (settings['dataset'] + '__Run3Summer22NanoAODv12__'
                  '130X_mcRun3_2022_realistic_v5__' + identity + '.root')
    output_name = 'pico_baseline_{}{}_{}.root'.format(sample, year, identity)
    if not SAFE_NAME.fullmatch(local_name):
        raise ValueError('Unsafe input alias')
    return dict(source=url, source_path=path, input_name=local_name,
                output_name=output_name)


def normalization_payload(sample, sumw, scale_sums):
    """The JSON directory/tag keys must exactly match process_nano -i/-f."""
    weights = dict(genEventSumw=sumw)
    weights.update({'LHEScaleSumw{}'.format(i): value
                    for i, value in enumerate(scale_sums)})
    return {NANO_INPUT_DIR: {sample['settings']['dataset']: weights}}
