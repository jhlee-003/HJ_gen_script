#!/bin/bash
# Source this instead of nano2pico's UCSB-specific set_env.sh.
BASELINE_HELPER_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
# This helper lives in HJ_gen_script/tools/baseline (and campaign sandboxes),
# not nano2pico/scripts. Use the actual CERN home, not relative parent paths.
# SConstruct also sources it in an env -i shell where HOME may be unset.
BASELINE_CERN_HOME="/afs/cern.ch/user/j/junhyuk"
BASELINE_CMSSW_SRC="$BASELINE_CERN_HOME/CMSSW_15_0_17/src"
export SCRAM_ARCH="el9_amd64_gcc12"
if [[ ! -d "$BASELINE_CMSSW_SRC" ]]; then
  echo "ERROR: CMSSW project not found: $BASELINE_CMSSW_SRC" >&2
  return 1
fi
# CMS setup scripts reference variables which may be unset.
BASELINE_HAD_NOUNSET=0
[[ $- == *u* ]] && BASELINE_HAD_NOUNSET=1
set +u
source /cvmfs/cms.cern.ch/cmsset_default.sh || return 1
BASELINE_PREVIOUS_DIR=$PWD
cd "$BASELINE_CMSSW_SRC" || return 1
BASELINE_SCRAM_RUNTIME=$(scramv1 runtime -sh) || return 1
eval "$BASELINE_SCRAM_RUNTIME"
cd "$BASELINE_PREVIOUS_DIR" || return 1
# Use the libraries from the selected CMSSW runtime, not SL7 overrides.
export PATH="$BASELINE_CERN_HOME/.local/bin:$PATH"
export PYTHONUNBUFFERED=1
# SConstruct starts a clean shell and sources this file to get its environment.
export SET_ENV_PATH="$BASELINE_HELPER_DIR/setup_runtime.sh"
[[ "$BASELINE_HAD_NOUNSET" == 1 ]] && set -u
python3 -c 'import ROOT' || return 1
unset BASELINE_HELPER_DIR BASELINE_CERN_HOME BASELINE_CMSSW_SRC
unset BASELINE_PREVIOUS_DIR BASELINE_SCRAM_RUNTIME BASELINE_HAD_NOUNSET
return 0
