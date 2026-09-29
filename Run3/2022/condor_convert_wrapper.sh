#!/bin/bash
set -euo pipefail

if [[ $# -ne 5 ]]; then
  echo "Usage: $0 CONDOR_PROCESS_ID NUMBER_OF_EVENTS ENV_FILE NUMBER_OF_THREADS JOB_OFFSET" >&2
  exit 2
fi

PROCESS_ID=$1
NEVENTS=$2
NAMES=$3
NTHREADS=$4
JOB_OFFSET=$5
SCRIPT="ProduceHJ_M-125_TuneCP5_2022.sh"

[[ "$PROCESS_ID" =~ ^[0-9]+$ ]] || { echo "ERROR: invalid Condor process ID: $PROCESS_ID" >&2; exit 2; }
[[ "$NEVENTS" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: invalid event count: $NEVENTS" >&2; exit 2; }
[[ "$NTHREADS" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: invalid thread count: $NTHREADS" >&2; exit 2; }
[[ "$JOB_OFFSET" =~ ^[0-9]+$ ]] || { echo "ERROR: invalid job offset: $JOB_OFFSET" >&2; exit 2; }
[[ -f "$SCRIPT" ]] || { echo "ERROR: missing production script: $SCRIPT" >&2; exit 2; }
[[ -f "$NAMES" ]] || { echo "ERROR: missing naming file: $NAMES" >&2; exit 2; }

JOBID=$((PROCESS_ID + JOB_OFFSET + 1))

source "./$NAMES"
: "${NANOAOD_NAME:?NANOAOD_NAME is missing from the naming file}"

if [[ -z "${X509_USER_PROXY:-}" || ! -r "$X509_USER_PROXY" ]]; then
  echo "ERROR: Condor did not provide a readable X.509 proxy" >&2
  exit 3
fi

echo "Running HJ production for Condor process $PROCESS_ID as job $JOBID"
bash "./$SCRIPT" "$JOBID" "$NEVENTS" "$NAMES" "$NTHREADS"

SOURCE_OUTPUT="${NANOAOD_NAME}__job.root"
CONDOR_OUTPUT="condor_output.root"
[[ -s "$SOURCE_OUTPUT" ]] || { echo "ERROR: missing final output: $SOURCE_OUTPUT" >&2; exit 20; }
[[ ! -e "$CONDOR_OUTPUT" ]] || { echo "ERROR: output already exists: $CONDOR_OUTPUT" >&2; exit 20; }
mv -- "$SOURCE_OUTPUT" "$CONDOR_OUTPUT"
[[ -s "$CONDOR_OUTPUT" ]] || { echo "ERROR: failed to prepare Condor output" >&2; exit 20; }

echo "Output ready for Condor transfer: $CONDOR_OUTPUT"
