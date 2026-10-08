#!/usr/bin/env bash
set -euo pipefail
trap 'echo "ERROR: baseline preparation/submission failed (exit=$?, line=$LINENO)" >&2' ERR

if [[ "$#" -lt 3 || "$#" -gt 5 ]]; then
    echo "Usage: $0 YEAR HJ|ggH|both CAMPAIGN [FILE_LIMIT] [--dry-run|--submit]" >&2
    exit 2
fi
YEAR="$1"
SAMPLE="$2"
CAMPAIGN="$3"
shift 3
LIMIT=0
MODE=prepare
for ARG in "$@"; do
    case "$ARG" in
        --dry-run) [[ "$MODE" == prepare ]] || exit 2; MODE=dry-run ;;
        --submit) [[ "$MODE" == prepare ]] || exit 2; MODE=submit ;;
        *) [[ "$LIMIT" == 0 && "$ARG" =~ ^[1-9][0-9]*$ ]] || exit 2; LIMIT="$ARG" ;;
    esac
done
[[ "$YEAR" == 2022 ]] || { echo "ERROR: only 2022 is configured" >&2; exit 2; }
[[ "$SAMPLE" == HJ || "$SAMPLE" == ggH || "$SAMPLE" == both ]] || exit 2
[[ "$CAMPAIGN" =~ ^[A-Za-z0-9][A-Za-z0-9_-]*$ ]] || exit 2

HELPER_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HELPER_DIR/setup_runtime.sh"
test -x "$NANO2PICO_DIR/run/process_nano.exe"
# The remote schedd cannot read a login node's /tmp proxy.
[[ "${X509_USER_PROXY:-}" == /afs/* && -r "$X509_USER_PROXY" ]] ||
    { echo "ERROR: export X509_USER_PROXY to a readable AFS proxy" >&2; exit 2; }
TIMELEFT="$(voms-proxy-info --file "$X509_USER_PROXY" --timeleft)"
[[ "$TIMELEFT" =~ ^[0-9]+$ && "$TIMELEFT" -ge 1800 ]] ||
    { echo "ERROR: proxy needs at least 30 minutes remaining" >&2; exit 2; }
export X509_USER_PROXY
[[ ! -e "$HELPER_DIR/runs/$CAMPAIGN" ]] ||
    { echo "ERROR: choose a new campaign name; existing campaign is not overwritten" >&2; exit 2; }
mkdir -p "$HELPER_DIR/runs/$CAMPAIGN"

SAMPLES=("$SAMPLE")
[[ "$SAMPLE" != both ]] || SAMPLES=(HJ ggH)
for SAMPLE in "${SAMPLES[@]}"; do
    case "$SAMPLE" in
        HJ)
            INPUT_DIR=/eos/user/j/junhyuk/ggF_MiNNLO_NanoAOD/2022
            INPUT_HOST=root://eosuser.cern.ch
            MAXDEPTH=1
            DATASET=GluGluHtoZG_Zto2L_M-125_TuneCP5_13p6TeV_powhegMiNNLO-pythia8
            OUTPUT_DIR=/eos/user/j/junhyuk/ggF_MiNNLO_NanoAOD/pico/HJ2022pico
            ;;
        ggH)
            INPUT_DIR=/eos/cms/store/mc/Run3Summer22NanoAODv12/GluGluHtoZG_Zto2L_M-125_TuneCP5_13p6TeV_powheg-pythia8/NANOAODSIM/130X_mcRun3_2022_realistic_v5-v2
            INPUT_HOST=root://eoscms.cern.ch
            MAXDEPTH=2
            DATASET=GluGluHtoZG_Zto2L_M-125_TuneCP5_13p6TeV_powheg-pythia8
            OUTPUT_DIR=/eos/user/j/junhyuk/ggF_MiNNLO_NanoAOD/pico/ggH2022pico
            ;;
    esac
    RUN_DIR="$HELPER_DIR/runs/$CAMPAIGN/$SAMPLE"
    mkdir -p "$RUN_DIR/logs"
    # One preparation pass per sample. Normalize exactly the files selected
    # for conversion, using their pre-baseline Runs metadata. No NanoAODs
    # are copied to AFS.
    python3 - "$YEAR" "$SAMPLE" "$INPUT_DIR" "$INPUT_HOST" "$MAXDEPTH" \
        "$DATASET" "$OUTPUT_DIR" "$RUN_DIR" "$HELPER_DIR" "$LIMIT" <<'PREPARE'
import hashlib
import json
import math
import os
from pathlib import Path
import sys
import ROOT

year, sample, input_dir, host, depth, dataset, output_dir, run_dir, helper, limit = sys.argv[1:]
ROOT.gROOT.SetBatch(True)
inputs, output, run = Path(input_dir), Path(output_dir), Path(run_dir)
if not inputs.is_dir():
    raise SystemExit("Input directory is unavailable: " + input_dir)
# This walk is bounded to the requested maxdepth; never descend into deeper folders.
files = []
for directory, children, names in os.walk(inputs, onerror=lambda error: (_ for _ in ()).throw(error)):
    level = len(Path(directory).relative_to(inputs).parts)
    if level + 1 >= int(depth):
        children[:] = []
    files.extend(Path(directory) / name for name in names if name.endswith(".root"))
files.sort()
if not files:
    raise SystemExit("No NanoAOD ROOT files found")
output.mkdir(parents=True, exist_ok=True)
if any(output.glob("*.root")):
    raise SystemExit("Output directory already contains ROOT files; archive/remove them or use a fresh output path: " + output_dir)
sources = [host + "/" + str(path) for path in files]
if any(any(char.isspace() for char in url) for url in sources):
    raise SystemExit("Input paths with whitespace are not supported")
print(f"{sample}: found {len(files)} files (maxdepth {depth})", flush=True)
sources = sources[:int(limit)] if int(limit) else sources
print(f"{sample}: selected {len(sources)} files for conversion and normalization", flush=True)

sumw, scales, events = 0.0, [0.0] * 9, 0
for index, url in enumerate(sources, 1):
    source = ROOT.TFile.Open(url, "READ")
    if not source or source.IsZombie():
        raise SystemExit("Cannot open NanoAOD: " + url)
    try:
        runs, tree = source.Get("Runs"), source.Get("Events")
        if not runs or runs.GetEntries() < 1 or not tree:
            raise SystemExit("Missing Runs/Events: " + url)
        for name in ("genEventSumw", "LHEScaleSumw"):
            if not runs.GetBranch(name):
                raise SystemExit("Missing " + name + ": " + url)
        events += int(tree.GetEntries())
        for entry in runs:
            value, variations = float(entry.genEventSumw), list(entry.LHEScaleSumw)
            if not math.isfinite(value) or len(variations) != 9:
                raise SystemExit("Invalid Runs normalization metadata: " + url)
            sumw += value
            for i, variation in enumerate(variations):
                variation = float(variation)
                if not math.isfinite(variation):
                    raise SystemExit("Nonfinite LHE scale sum: " + url)
                scales[i] += variation
    finally:
        source.Close()
    if index % 100 == 0 or index == len(sources):
        print(f"Normalization: {index}/{len(sources)} files", flush=True)
if any(not math.isfinite(value) or value == 0 for value in [sumw] + scales):
    raise SystemExit("Zero/nonfinite sample normalization")
weights = {"genEventSumw": sumw}
weights.update({f"LHEScaleSumw{i}": value for i, value in enumerate(scales)})
(run / "normalization.json").write_text(json.dumps({"dataset": dataset, "weights": weights}, indent=2) + "\n")
(run / "inputs.txt").write_text("\n".join(sources) + "\n")

rows = []
for url in sources:
    # A short ID is only for unique filenames, not a provenance/version check.
    file_id = hashlib.sha256(url.encode()).hexdigest()[:20]
    alias = f"{dataset}__Run3Summer22NanoAODv12__130X_mcRun3_2022_realistic_v5__{file_id}.root"
    pico = f"pico_baseline_{sample}{year}_{file_id}.root"
    rows.append(f"{url}\t{url.rsplit('/', 1)[1]}\t{alias}\t{pico}")
(run / "jobs.tsv").write_text("\n".join(rows) + "\n")
header = f"YEAR = {year}\nSAMPLE = {sample}\nHELPER_DIR = {helper}\nRUN_DIR = {run}\nOUTPUT_DESTINATION = root://eosuser.cern.ch/{output_dir}/\n"
template = (Path(helper) / "baseline.sub").read_text()
queue = f"\nqueue source_url, source_name, nano_name, pico_name from {run}/jobs.tsv\n"
(run / "baseline.sub").write_text(header + template + queue)
print(f"{sample}: {events:,} events, signed sumw={sumw:.12g}; prepared {len(rows)} jobs", flush=True)
PREPARE
done

# Prepare both samples successfully before submitting either.
for SAMPLE in "${SAMPLES[@]}"; do
    RUN_DIR="$HELPER_DIR/runs/$CAMPAIGN/$SAMPLE"
    case "$MODE" in
        submit) condor_submit "$RUN_DIR/baseline.sub" ;;
        dry-run) condor_submit -dry-run "$RUN_DIR/baseline.classad" "$RUN_DIR/baseline.sub" ;;
        prepare) echo "Prepared: $RUN_DIR/baseline.sub (not submitted)" ;;
    esac
done
