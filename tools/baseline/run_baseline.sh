#!/usr/bin/env bash
set -euo pipefail
trap 'echo "ERROR: baseline worker failed (exit=$?, line=$LINENO)" >&2' ERR

if [[ "$#" != 5 ]]; then
    echo "Usage: $0 YEAR LOCAL_NANO NANO_NAME LOCAL_PICO NORMALIZATION_JSON" >&2
    exit 2
fi
YEAR="$1"
LOCAL_NANO="$2"
NANO_NAME="$3"
LOCAL_PICO="$4"
NORMALIZATION="$5"
[[ "$YEAR" == 2022 ]] || { echo "ERROR: unsupported year: $YEAR" >&2; exit 2; }
[[ "$NANO_NAME" != */* && "$NANO_NAME" == *.root &&
   "$LOCAL_PICO" != */* && "$LOCAL_PICO" == *.root ]] ||
    { echo "ERROR: Nano alias/output must be ROOT basenames" >&2; exit 2; }
test -s "$LOCAL_NANO"
test -s "$NORMALIZATION"
[[ ! -e "$LOCAL_PICO" ]] || { echo "ERROR: output already exists: $LOCAL_PICO" >&2; exit 2; }

HELPER_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HELPER_DIR/setup_runtime.sh"
test -x "$NANO2PICO_DIR/run/process_nano.exe"
LOCAL_NANO="$(realpath "$LOCAL_NANO")"
NORMALIZATION="$(realpath "$NORMALIZATION")"
LOCAL_PICO="$PWD/$LOCAL_PICO"

#-----------------------------Prepare job--------------------------------
WORK_DIR="$(mktemp -d "$PWD/baseline_work.XXXXXX")"
# Only this worker-created scratch directory is removed.
trap 'rm -rf -- "$WORK_DIR"' EXIT
WORK_NANO_DIR="$WORK_DIR/NanoAODv12/nano/$YEAR/mc"
N2P_OUT_DIR="$WORK_DIR/zgamma"
RAW_PICO="$N2P_OUT_DIR/raw_pico/raw_pico_$NANO_NAME"
BASELINE_PICO="$WORK_DIR/baseline.root"
mkdir -p "$WORK_NANO_DIR" "$N2P_OUT_DIR/raw_pico"
ln -- "$LOCAL_NANO" "$WORK_NANO_DIR/$NANO_NAME" 2>/dev/null ||
    cp -- "$LOCAL_NANO" "$WORK_NANO_DIR/$NANO_NAME"

# The selected-input sums do not change: only adapt the JSON directory key to
# this job's --in_dir, as required by process_nano.
python3 - "$NORMALIZATION" "$WORK_NANO_DIR" "$NANO_NAME" "$WORK_DIR/norm.json" <<'NORM'
import json
import sys
source, input_dir, nano_name, output = sys.argv[1:]
with open(source) as handle:
    norm = json.load(handle)
dataset = nano_name.split("__", 1)[0]
if dataset != norm["dataset"]:
    raise SystemExit("Normalization dataset does not match NanoAOD alias")
with open(output, "w") as handle:
    json.dump({input_dir: {dataset: norm["weights"]}}, handle, allow_nan=False)
NORM

#-----------------------------Convert to pico----------------------------
# Run from the existing checkout so relative correction-data paths work.
cd "$NANO2PICO_DIR"
./run/process_nano.exe \
    --in_dir "$WORK_NANO_DIR" \
    --in_file "$NANO_NAME" \
    --out_dir "$N2P_OUT_DIR" \
    --nent -1 \
    --norm "$WORK_DIR/norm.json"
test -s "$RAW_PICO"

#-----------------------------Baseline selection-------------------------
python3 - "$RAW_PICO" "$BASELINE_PICO" <<'SKIM'
import sys
import ROOT
ROOT.gROOT.SetBatch(True)
source = ROOT.TFile.Open(sys.argv[1], "READ")
if not source or source.IsZombie():
    raise SystemExit("Cannot open raw pico")
tree = source.Get("tree")
if not tree:
    raise SystemExit("Raw pico is missing tree 'tree'")
for name in ("use_event", "zg_cutBitMap", "weight", "w_lumi"):
    if not tree.GetBranch(name):
        raise SystemExit("Missing pico branch: " + name)
cut = "use_event && (zg_cutBitMap==3582 || zg_cutBitMap==3583 || zg_cutBitMap==3070 || zg_cutBitMap==3071)"
output = ROOT.TFile.Open(sys.argv[2], "RECREATE")
if not output or output.IsZombie():
    raise SystemExit("Cannot create baseline pico")
output.cd()
selected = tree.CopyTree(cut)
if not selected or selected.Write("tree", ROOT.TObject.kOverwrite) <= 0:
    raise SystemExit("Baseline tree write failed")
print(f"Raw pico events: {tree.GetEntries()}")
print(f"Baseline events: {selected.GetEntries()}")
output.Close()
source.Close()
SKIM

#-----------------------------Save file----------------------------------
# Condor transfers this single file to the sample's EOS output directory.
test -s "$BASELINE_PICO"
mv -- "$BASELINE_PICO" "$LOCAL_PICO"
echo "Saved $LOCAL_PICO"
