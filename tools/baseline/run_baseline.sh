#!/bin/bash
set -euo pipefail
HELPER_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
source "$HELPER_DIR/setup_runtime.sh"
if [[ -z "${X509_USER_PROXY:-}" || ! -r "$X509_USER_PROXY" ]]; then
  echo "ERROR: Condor did not delegate a readable X.509 proxy" >&2
  exit 3
fi
# The executable and data are transferred, not accessed via an AFS checkout.
# Capture a complete combined log. ON_SUCCESS deliberately skips output
# transfers for failures, so upload diagnostics explicitly before returning the
# original failure status; otherwise a held job could have no returned stderr.
set +e
python3 "$HELPER_DIR/baseline_worker.py" "$@" 2>&1 | tee baseline_worker.log
BASELINE_STATUS=${PIPESTATUS[0]}
set -e
if [[ "$BASELINE_STATUS" -ne 0 && $# -ge 4 ]]; then
  if FAILURE_URL=$(python3 - "$2" "$3" "$1" "$4" "${BASELINE_JOB_ID:-unknown}" <<'PY'
import json, re, sys
with open(sys.argv[1]) as stream:
    sample = json.load(stream)['samples'][sys.argv[2]]
suffix = '_'.join(sys.argv[2:])
if not re.fullmatch(r'[A-Za-z0-9_.-]+', suffix):
    raise SystemExit('Unsafe diagnostic filename')
print('root://eosuser.cern.ch/' + sample['settings']['output_dir']
      + '/logs/baseline_failure_' + suffix + '.log')
PY
  ); then
    xrdcp --nopbar baseline_worker.log "$FAILURE_URL" || {
      echo "WARNING: failed to upload diagnostic log; original exit=$BASELINE_STATUS" >&2
    }
  fi
fi
exit "$BASELINE_STATUS"
