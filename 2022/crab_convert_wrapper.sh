#!/bin/bash
set -euo pipefail

if [[ $# -lt 1 || ! "$1" =~ ^[0-9]+$ ]]; then
  echo "ERROR: CRAB must supply a non-negative integer job index" >&2
  exit 2
fi
JOB_INDEX=$1
shift

SCRIPT=""
EVENTS=""
NAMES=""
THREADS=""

for arg in "$@"; do
  case "$arg" in
    script=*) SCRIPT="${arg#*=}" ;;
    events=*) EVENTS="${arg#*=}" ;;
    names=*) NAMES="${arg#*=}" ;;
    threads=*) THREADS="${arg#*=}" ;;
    *) echo "ERROR: unknown argument: $arg" >&2; exit 2 ;;
  esac
done

[[ -n "$SCRIPT" && -f "$SCRIPT" ]] || { echo "ERROR: missing production script" >&2; exit 2; }
[[ -n "$NAMES" && -f "$NAMES" ]] || { echo "ERROR: missing naming file" >&2; exit 2; }
[[ "$EVENTS" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: invalid events=$EVENTS" >&2; exit 2; }
[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: invalid threads=$THREADS" >&2; exit 2; }

source "./$NAMES"
: "${NANOAOD_NAME:?NANOAOD_NAME is missing from the naming file}"

# Keep a valid report if setup fails before the first production cmsRun.
cmsRun -j FrameworkJobReport.xml PSet.py
[[ -s FrameworkJobReport.xml ]] || { echo "ERROR: missing bootstrap job report" >&2; exit 90; }

JOBNUM=$((JOB_INDEX + 1))
export NTHREADS="$THREADS"
echo "Running HJ production for CRAB index $JOB_INDEX as job $JOBNUM ($EVENTS events, $NTHREADS threads)"
bash "./$SCRIPT" "$JOB_INDEX" "$EVENTS" "$NAMES"

PRODUCED_ROOT="${NANOAOD_NAME}__job-${JOBNUM}.root"
OUTROOT="${NANOAOD_NAME}__job.root"
[[ -s "$PRODUCED_ROOT" ]] || { echo "ERROR: missing final output: $PRODUCED_ROOT" >&2; exit 20; }
[[ ! -e "$OUTROOT" ]] || { echo "ERROR: output already exists: $OUTROOT" >&2; exit 20; }
mv -- "$PRODUCED_ROOT" "$OUTROOT"

# The last cmsRun report contains the per-job filename. Keep it consistent with
# the fixed filename declared in the CRAB configuration.
sed -i "s|${PRODUCED_ROOT}|${OUTROOT}|g" FrameworkJobReport.xml
[[ -s FrameworkJobReport.xml ]] || { echo "ERROR: missing production job report" >&2; exit 90; }

NEV="$(root -l -b -q -e "TFile f(\"$OUTROOT\"); auto t=(TTree*)f.Get(\"Events\"); if(!t){std::cout<<0; gSystem->Exit(0);} std::cout<<t->GetEntries(); gSystem->Exit(0);" 2>/dev/null | tail -n 1 | tr -d '[:space:]')"
[[ "$NEV" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: Events tree is empty in $OUTROOT" >&2; exit 21; }

echo "Events = $NEV"
echo "Output ready for CRAB transfer to T3_KR_KNU: $OUTROOT"
