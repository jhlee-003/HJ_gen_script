#!/usr/bin/env bash
# Source from lxplus or a CMS EL9 worker; also used by SCons via SET_ENV_PATH.
BASELINE_START_DIR="$PWD"
BASELINE_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASELINE_NOUNSET=0
[[ $- == *u* ]] && BASELINE_NOUNSET=1
set +u
if source /cvmfs/cms.cern.ch/cmsset_default.sh &&
   export SCRAM_ARCH=el9_amd64_gcc12 &&
   cd /afs/cern.ch/user/j/junhyuk/CMSSW_15_0_17/src &&
   BASELINE_RUNTIME="$(scramv1 runtime -sh)" &&
   eval "$BASELINE_RUNTIME"; then
    BASELINE_STATUS=0
else
    BASELINE_STATUS=$?
fi
cd "$BASELINE_START_DIR" || return 1
[[ "$BASELINE_NOUNSET" == 1 ]] && set -u
if [[ "$BASELINE_STATUS" != 0 ]]; then
    echo "ERROR: CMSSW setup failed (exit=$BASELINE_STATUS)" >&2
    return "$BASELINE_STATUS"
fi

export NANO2PICO_DIR=/afs/cern.ch/user/j/junhyuk/nano2pico_sequoia_v1
export PATH="/afs/cern.ch/user/j/junhyuk/.local/bin:$PATH"
export PYTHONUNBUFFERED=1
export SET_ENV_PATH="$BASELINE_SCRIPT_DIR/setup_runtime.sh"
export XrdSecPROTOCOL=gsi
unset BASELINE_START_DIR BASELINE_SCRIPT_DIR BASELINE_RUNTIME BASELINE_NOUNSET BASELINE_STATUS
return 0
