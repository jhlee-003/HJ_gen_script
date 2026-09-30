#!/bin/bash
set -euo pipefail

if [[ $# -ne 3 ]]; then
  echo "[Usage] $0 JOB_INDEX NUMBER_OF_EVENTS ENV_FILE" >&2
  echo "  JOB_INDEX starts from 0; output file numbering starts from 1." >&2
  echo "  Set NTHREADS in the environment to override the default of 2." >&2
  exit 2
fi

JOB_INDEX=$1
NEVENTS=$2
ENV_INPUT=$3
NTHREADS=${NTHREADS:-2}

[[ "$JOB_INDEX" =~ ^[0-9]+$ ]] || { echo "ERROR: JOB_INDEX must be a non-negative integer" >&2; exit 2; }
[[ "$NEVENTS" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: NUMBER_OF_EVENTS must be a positive integer" >&2; exit 2; }
[[ "$NTHREADS" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: NTHREADS must be a positive integer" >&2; exit 2; }
[[ -f "$ENV_INPUT" ]] || { echo "ERROR: ENV_FILE does not exist: $ENV_INPUT" >&2; exit 2; }

WORKDIR=$(pwd -P)
ENV_DIR=$(cd -- "$(dirname -- "$ENV_INPUT")" && pwd -P)
ENV_FILE="$ENV_DIR/$(basename -- "$ENV_INPUT")"
source "$ENV_FILE"
: "${BASE_TAG:?}" "${Fragment_filename:?}" "${AOD_NAME:?}" "${MINIAOD_NAME:?}" "${NANOAOD_NAME:?}"

JOBNUM=$((JOB_INDEX + 1))
TAG="${BASE_TAG}__job-${JOBNUM}"

# CRAB flattens inputFiles into the worker directory. A local checkout keeps the
# fragment under config/, so support both layouts without moving the source file.
if [[ -f "$WORKDIR/$Fragment_filename" ]]; then
  FRAGMENT_PATH="$WORKDIR/$Fragment_filename"
elif [[ -f "$WORKDIR/config/$Fragment_filename" ]]; then
  FRAGMENT_PATH="$WORKDIR/config/$Fragment_filename"
else
  echo "ERROR: fragment does not exist: $Fragment_filename" >&2
  exit 3
fi

mkdir -p job_scripts
JOB_SCRIPT="$WORKDIR/job_scripts/${TAG}_cmd.sh"

cat > "$JOB_SCRIPT" <<'EndOfJobScript'
#!/bin/bash
set -euo pipefail

if [[ $# -ne 6 ]]; then
  echo "ERROR: internal job script received the wrong number of arguments" >&2
  exit 2
fi

JOBNUM=$1
NEVENTS=$2
NTHREADS=$3
WORKDIR=$4
ENV_FILE=$5
FRAGMENT_PATH=$6

cd "$WORKDIR"
source "$ENV_FILE"
: "${BASE_TAG:?}" "${Fragment_filename:?}" "${AOD_NAME:?}" "${MINIAOD_NAME:?}" "${NANOAOD_NAME:?}"

TAG="${BASE_TAG}__job-${JOBNUM}"
SIM_FILE="${AOD_NAME}__job-${JOBNUM}__SIM.root"
HLT_FILE="${AOD_NAME}__job-${JOBNUM}__HLT.root"
AOD_FILE="${AOD_NAME}__job-${JOBNUM}.root"
MINIAOD_FILE="${MINIAOD_NAME}__job-${JOBNUM}.root"
NANOAOD_FILE="${NANOAOD_NAME}__job-${JOBNUM}.root"

# Use a fresh runtime subshell for each CMSSW release while keeping all outputs
# in the scheduler's working directory.
run_in_cmssw() (
  local release=$1
  local arch=$2
  shift 2

  set +u
  if [[ -n "${CMSSW_BASE:-}" ]]; then
    eval "$(scram unsetenv -sh)"
  fi
  export SCRAM_ARCH="$arch"
  source /cvmfs/cms.cern.ch/cmsset_default.sh
  cd "$WORKDIR"
  if [[ ! -d "$release/src" ]]; then
    scram project CMSSW "$release"
  fi
  cd "$release/src"
  eval "$(scram runtime -sh)"
  set -u
  cd "$WORKDIR"
  "$@"
)

run_lhe_gen_sim() {
  echo "--------------------------------LHE,GEN,SIM-------------------------------------"
  echo "CMSSW_12_4_11_patch3 with el8_amd64_gcc10"

  mkdir -p "$CMSSW_BASE/src/Configuration/GenProduction/python"
  cp -- "$FRAGMENT_PATH" "$CMSSW_BASE/src/Configuration/GenProduction/python/$Fragment_filename"
  (cd "$CMSSW_BASE/src" && scram b -j "$NTHREADS")

  local customise="process.source.numberEventsInLuminosityBlock=cms.untracked.uint32(200); process.source.firstRun=cms.untracked.uint32(${JOBNUM}); from IOMC.RandomEngine.RandomServiceHelper import RandomNumberServiceHelper; randSvc=RandomNumberServiceHelper(process.RandomNumberGeneratorService); randSvc.populate()"
  cmsDriver.py "Configuration/GenProduction/python/$Fragment_filename" \
    --era Run3 \
    --customise Configuration/DataProcessing/Utils.addMonitoring \
    --beamspot Realistic25ns13p6TeVEarly2022Collision \
    --step LHE,GEN,SIM \
    --geometry DB:Extended \
    --conditions 124X_mcRun3_2022_realistic_v12 \
    --customise_commands "$customise" \
    --datatier GEN-SIM,LHE \
    --eventcontent RAWSIM,LHE \
    --python_filename "${TAG}__LHE__cfg.py" \
    --fileout "file:$SIM_FILE" \
    --number "$NEVENTS" \
    --nThreads "$NTHREADS" \
    --no_exec --mc
  cmsRun -j FrameworkJobReport.xml "${TAG}__LHE__cfg.py"
}

run_digi_reco() {
  echo "---------------------------------DIGIPREMIX------------------------------------"
  echo "CMSSW_12_4_11_patch3 with el8_amd64_gcc10"

  cmsDriver.py \
    --era Run3 \
    --customise Configuration/DataProcessing/Utils.addMonitoring \
    --procModifiers premix_stage2,siPixelQualityRawToDigi \
    --datamix PreMix \
    --step DIGI,DATAMIX,L1,DIGI2RAW,HLT:2022v12 \
    --geometry DB:Extended \
    --conditions 124X_mcRun3_2022_realistic_v12 \
    --datatier GEN-SIM-RAW \
    --eventcontent PREMIXRAW \
    --customise_commands "from IOMC.RandomEngine.RandomServiceHelper import RandomNumberServiceHelper; RandomNumberServiceHelper(process.RandomNumberGeneratorService).populate()" \
    --python_filename "${TAG}__DIGIPREMIX__cfg.py" \
    --fileout "file:$HLT_FILE" \
    --filein "file:$SIM_FILE" \
    --number -1 \
    --nThreads "$NTHREADS" \
    --pileup_input "dbs:/Neutrino_E-10_gun/Run3Summer21PrePremix-Summer22_124X_mcRun3_2022_realistic_v11-v2/PREMIX" \
    --no_exec --mc
  cmsRun -j FrameworkJobReport.xml "${TAG}__DIGIPREMIX__cfg.py"
  rm -f -- "$SIM_FILE" "${SIM_FILE%.root}_inLHE.root"

  echo "------------------------------------RECO----------------------------------------"
  cmsDriver.py \
    --era Run3 \
    --customise Configuration/DataProcessing/Utils.addMonitoring \
    --procModifiers siPixelQualityRawToDigi \
    --step RAW2DIGI,L1Reco,RECO,RECOSIM \
    --geometry DB:Extended \
    --conditions 124X_mcRun3_2022_realistic_v12 \
    --datatier AODSIM \
    --eventcontent AODSIM \
    --python_filename "${TAG}__AOD__cfg.py" \
    --fileout "file:$AOD_FILE" \
    --filein "file:$HLT_FILE" \
    --number -1 \
    --nThreads "$NTHREADS" \
    --no_exec --mc
  cmsRun -j FrameworkJobReport.xml "${TAG}__AOD__cfg.py"
  rm -f -- "$HLT_FILE"
}

run_mini_nano() {
  echo "----------------------------------MiniAODv4-------------------------------------"
  echo "CMSSW_13_0_13 with el8_amd64_gcc11"

  cmsDriver.py \
    --era Run3,run3_miniAOD_12X \
    --customise Configuration/DataProcessing/Utils.addMonitoring \
    --step PAT \
    --geometry DB:Extended \
    --conditions 130X_mcRun3_2022_realistic_v5 \
    --datatier MINIAODSIM \
    --eventcontent MINIAODSIM \
    --python_filename "${TAG}__MiniAODv4__cfg.py" \
    --fileout "file:$MINIAOD_FILE" \
    --filein "file:$AOD_FILE" \
    --number -1 \
    --nThreads "$NTHREADS" \
    --no_exec --mc
  cmsRun -j FrameworkJobReport.xml "${TAG}__MiniAODv4__cfg.py"
  rm -f -- "$AOD_FILE"

  echo "-----------------------------------NanoAODv12-----------------------------------"
  cmsDriver.py \
    --scenario pp \
    --era Run3 \
    --customise Configuration/DataProcessing/Utils.addMonitoring \
    --step NANO \
    --conditions 130X_mcRun3_2022_realistic_v5 \
    --datatier NANOAODSIM \
    --eventcontent NANOAODSIM \
    --python_filename "${TAG}__NanoAODv12__cfg.py" \
    --fileout "file:$NANOAOD_FILE" \
    --filein "file:$MINIAOD_FILE" \
    --number -1 \
    --nThreads "$NTHREADS" \
    --no_exec --mc
  cmsRun -j FrameworkJobReport.xml "${TAG}__NanoAODv12__cfg.py"

  [[ -s "$NANOAOD_FILE" ]] || { echo "ERROR: missing NanoAOD output: $NANOAOD_FILE" >&2; exit 20; }
  root -l -b -q -e "TFile f(\"$NANOAOD_FILE\"); if(f.IsZombie()){gSystem->Exit(21);} auto t=(TTree*)f.Get(\"Events\"); if(!t || t->GetEntries()==0){gSystem->Exit(21);} std::cout << \"Events = \" << t->GetEntries() << std::endl; gSystem->Exit(0);"
  rm -f -- "$MINIAOD_FILE"
}

run_in_cmssw CMSSW_12_4_11_patch3 el8_amd64_gcc10 run_lhe_gen_sim
run_in_cmssw CMSSW_12_4_11_patch3 el8_amd64_gcc10 run_digi_reco
run_in_cmssw CMSSW_13_0_13 el8_amd64_gcc11 run_mini_nano

echo "Production completed: $NANOAOD_FILE"
EndOfJobScript

chmod +x "$JOB_SCRIPT"
echo "Made job_scripts/${TAG}_cmd.sh"
"$JOB_SCRIPT" "$JOBNUM" "$NEVENTS" "$NTHREADS" "$WORKDIR" "$ENV_FILE" "$FRAGMENT_PATH"
