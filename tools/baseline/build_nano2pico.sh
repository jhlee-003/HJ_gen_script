#!/bin/bash
set -euo pipefail
HELPER_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
source "$HELPER_DIR/setup_runtime.sh"
NANO2PICO_DIR=$(python3 - "$HELPER_DIR/baseline_config.json" <<'PY'
import json, sys
with open(sys.argv[1]) as stream:
    print(json.load(stream)['nano2pico_dir'])
PY
)
if [[ ! -f "$NANO2PICO_DIR/SConstruct" ]]; then
  echo "ERROR: nano2pico checkout not found: $NANO2PICO_DIR" >&2
  exit 2
fi
command -v scons >/dev/null || {
  echo "ERROR: install SCons in this runtime: python3 -m pip install --user scons" >&2
  exit 2
}
cd "$NANO2PICO_DIR"
echo "Building nano2pico with $CMSSW_VERSION / $SCRAM_ARCH"
# No checkout, pull, source changes, or UCSB queue setup are performed here.
scons -j "${BASELINE_BUILD_CPUS:-2}"
