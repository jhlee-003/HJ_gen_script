#!/bin/bash
# Source this instead of nano2pico's UCSB-specific set_env.sh.
BASELINE_HELPER_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
# A clean CC7 build shell may have no Python 3 until scram runtime is loaded.
# Read these three simple string fields with Bash, without eval or Python.
BASELINE_CONFIG_TEXT=$(< "$BASELINE_HELPER_DIR/baseline_config.json")
for BASELINE_KEY in cmssw_version scram_arch correctionlib_dir; do
  BASELINE_FIELD_PATTERN='"'"$BASELINE_KEY"'"[[:space:]]*:[[:space:]]*"([^"\\]+)"'
  if [[ ! "$BASELINE_CONFIG_TEXT" =~ $BASELINE_FIELD_PATTERN ]]; then
    echo "ERROR: missing or invalid runtime field: $BASELINE_KEY" >&2
    return 1
  fi
  case "$BASELINE_KEY" in
    cmssw_version) BASELINE_CMSSW_VERSION=${BASH_REMATCH[1]} ;;
    scram_arch) BASELINE_SCRAM_ARCH=${BASH_REMATCH[1]} ;;
    correctionlib_dir) BASELINE_CORRECTIONLIB_DIR=${BASH_REMATCH[1]} ;;
  esac
done
if [[ ! "$BASELINE_CMSSW_VERSION" =~ ^CMSSW_[A-Za-z0-9_]+$ ||
      ! "$BASELINE_SCRAM_ARCH" =~ ^[A-Za-z0-9_]+$ ||
      ! "$BASELINE_CORRECTIONLIB_DIR" =~ ^/cvmfs/cms[.]cern[.]ch/[A-Za-z0-9_./-]+$ ]]; then
  echo "ERROR: invalid CMSSW runtime settings" >&2
  return 1
fi
export SCRAM_ARCH="$BASELINE_SCRAM_ARCH"
BASELINE_CMSSW_SRC="/cvmfs/cms.cern.ch/$SCRAM_ARCH/cms/cmssw/$BASELINE_CMSSW_VERSION/src"
if [[ ! -d "$BASELINE_CMSSW_SRC" || ! -d "$BASELINE_CORRECTIONLIB_DIR" ]]; then
  echo "ERROR: configured CMSSW/correctionlib runtime is missing from CVMFS" >&2
  return 1
fi
# CMS setup scripts reference variables which may be unset.
BASELINE_HAD_NOUNSET=0
[[ $- == *u* ]] && BASELINE_HAD_NOUNSET=1
set +u
source /cvmfs/cms.cern.ch/cmsset_default.sh || return 1
BASELINE_PREVIOUS_DIR=$PWD
cd "$BASELINE_CMSSW_SRC" || return 1
BASELINE_SCRAM_RUNTIME=$(scram runtime -sh) || return 1
eval "$BASELINE_SCRAM_RUNTIME"
cd "$BASELINE_PREVIOUS_DIR" || return 1
export LD_LIBRARY_PATH="$BASELINE_CORRECTIONLIB_DIR:${LD_LIBRARY_PATH:-}"
export PATH="$HOME/.local/bin:$PATH"
export PYTHONUNBUFFERED=1
# SConstruct starts a clean shell and sources this file to get its environment.
export SET_ENV_PATH="$BASELINE_HELPER_DIR/setup_runtime.sh"
[[ "$BASELINE_HAD_NOUNSET" == 1 ]] && set -u
python3 -c 'import ROOT' || return 1
return 0
