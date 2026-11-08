#!/bin/bash
set -euo pipefail

if [[ $# -ne 6 ]]; then
  echo "Usage: $0 CONDOR_PROCESS_ID NUMBER_OF_EVENTS ENV_FILE NUMBER_OF_THREADS JOB_OFFSET STAGEOUT_FILE" >&2
  exit 2
fi

PROCESS_ID=$1
NEVENTS=$2
NAMES=$3
NTHREADS=$4
JOB_OFFSET=$5
STAGEOUT_FILE=$6
SCRIPT="ProduceHJ_M-125_TuneCP5_2022.sh"

[[ "$PROCESS_ID" =~ ^[0-9]+$ ]] || { echo "ERROR: invalid Condor process ID: $PROCESS_ID" >&2; exit 2; }
[[ "$NEVENTS" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: invalid event count: $NEVENTS" >&2; exit 2; }
[[ "$NTHREADS" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: invalid thread count: $NTHREADS" >&2; exit 2; }
[[ "$JOB_OFFSET" =~ ^[0-9]+$ ]] || { echo "ERROR: invalid job offset: $JOB_OFFSET" >&2; exit 2; }
[[ "$STAGEOUT_FILE" =~ ^[A-Za-z0-9._-]+[.]root$ ]] || { echo "ERROR: invalid stageout filename: $STAGEOUT_FILE" >&2; exit 2; }
[[ -f "$SCRIPT" ]] || { echo "ERROR: missing production script: $SCRIPT" >&2; exit 2; }
[[ -f "$NAMES" ]] || { echo "ERROR: missing naming file: $NAMES" >&2; exit 2; }

JOB_INDEX=$((PROCESS_ID + JOB_OFFSET))
JOBNUM=$((JOB_INDEX + 1))

source "./$NAMES"
: "${NANOAOD_NAME:?NANOAOD_NAME is missing from the naming file}"

if [[ -z "${X509_USER_PROXY:-}" || ! -r "$X509_USER_PROXY" ]]; then
  echo "ERROR: Condor did not provide a readable X.509 proxy" >&2
  exit 3
fi

export NTHREADS
echo "Running HJ production for Condor process $PROCESS_ID as index $JOB_INDEX / job $JOBNUM"
bash "./$SCRIPT" "$JOB_INDEX" "$NEVENTS" "$NAMES"

PRODUCED_ROOT="${NANOAOD_NAME}__job-${JOBNUM}.root"
[[ -s "$PRODUCED_ROOT" ]] || { echo "ERROR: missing final output: $PRODUCED_ROOT" >&2; exit 20; }
[[ ! -e "$STAGEOUT_FILE" ]] || { echo "ERROR: output already exists: $STAGEOUT_FILE" >&2; exit 20; }
mv -- "$PRODUCED_ROOT" "$STAGEOUT_FILE"
[[ -s "$STAGEOUT_FILE" ]] || { echo "ERROR: failed to prepare Condor output" >&2; exit 20; }

echo "Output ready for HTCondor stageout to EOS: $STAGEOUT_FILE"
