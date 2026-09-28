#!/bin/bash
set -euo pipefail

if [[ $# -lt 1 || ! "$1" =~ ^[1-9][0-9]*$ ]]; then
  echo "ERROR: CRAB must supply a positive integer job ID" >&2
  exit 2
fi
JOBID=$1
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

# CRAB flattens inputFiles into the job directory and provides its proxy.
source "./$NAMES"
: "${NANOAOD_NAME:?NANOAOD_NAME is missing from the naming file}"

# Keep a valid report if setup fails before the first production cmsRun.
cmsRun -j FrameworkJobReport.xml PSet.py
[[ -s FrameworkJobReport.xml ]] || { echo "ERROR: missing bootstrap job report" >&2; exit 90; }

echo "Running HJ production for CRAB job $JOBID ($EVENTS events, $THREADS threads)"
bash "./$SCRIPT" "$JOBID" "$EVENTS" "$NAMES" "$THREADS"

OUTROOT="${NANOAOD_NAME}__job.root"
[[ -s "$OUTROOT" ]] || { echo "ERROR: missing final output: $OUTROOT" >&2; exit 20; }
[[ -s FrameworkJobReport.xml ]] || { echo "ERROR: missing production job report" >&2; exit 90; }
echo "Output ready for CRAB transfer to T3_KR_KNU: $OUTROOT"
