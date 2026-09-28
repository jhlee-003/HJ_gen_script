#!/bin/bash
set -euo pipefail

if [[ $# -ne 4 ]]; then
  echo "Usage: $0 CRAB_JOB_ID NUMBER_OF_EVENTS ENV_FILE NUMBER_OF_THREADS" >&2
  exit 2
fi
[[ "$1" =~ ^[1-9][0-9]*$ && "$2" =~ ^[1-9][0-9]*$ && "$4" =~ ^[1-9][0-9]*$ ]] || {
  echo "ERROR: job ID, event count, and thread count must be positive integers" >&2
  exit 2
}
[[ -f "$3" ]] || { echo "ERROR: naming file not found: $3" >&2; exit 2; }

JOBNUM=$1
NEVENTS=$2
NTHREADS=$4
source "$3"
: "${BASE_TAG:?}" "${Fragment_filename:?}" "${AOD_NAME:?}" "${MINIAOD_NAME:?}" "${NANOAOD_NAME:?}"

WORKDIR=$(pwd -P)
TAG="${BASE_TAG}__job-${JOBNUM}"
# CRAB inputFiles arrive as basenames; config/ also supports a local pilot.
if [[ -f "$Fragment_filename" ]]; then
  FRAGMENT_PATH="$WORKDIR/$Fragment_filename"
elif [[ -f "config/$Fragment_filename" ]]; then
  FRAGMENT_PATH="$WORKDIR/config/$Fragment_filename"
else
  echo "ERROR: fragment not found: $Fragment_filename" >&2
  exit 3
fi

SIM_FILE="${TAG}__SIM.root"
HLT_FILE="${TAG}__HLT.root"
AOD_FILE="${AOD_NAME}__job-${JOBNUM}.root"
MINIAOD_FILE="${MINIAOD_NAME}__job-${JOBNUM}.root"
# CRAB adds its job suffix during stageout. Keep the declared output name fixed.
NANOAOD_FILE="${NANOAOD_NAME}__job.root"

# Each stage has a clean CMSSW environment; no copied proxy is needed.
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
  echo "LHE,GEN,SIM: CMSSW_12_4_11_patch3"
  # Central recipe: HIG-Run3Summer22wmLHEGS-00133; use our local HJ fragment.
  mkdir -p "$CMSSW_BASE/src/Configuration/GenProduction/python"
  cp -- "$FRAGMENT_PATH" "$CMSSW_BASE/src/Configuration/GenProduction/python/$Fragment_filename"
  (cd "$CMSSW_BASE/src" && scram b -j "$NTHREADS")

  # Fresh seeds for all modules, including the HJ external LHE producer.
  # Give each CRAB job a distinct run number; filenames contain no timestamps.
  local customise="from IOMC.RandomEngine.RandomServiceHelper import RandomNumberServiceHelper; RandomNumberServiceHelper(process.RandomNumberGeneratorService).populate(); process.source.firstRun=cms.untracked.uint32(${JOBNUM}); process.source.numberEventsInLuminosityBlock=cms.untracked.uint32(200)"
  cmsDriver.py "Configuration/GenProduction/python/$Fragment_filename" \
    --eventcontent RAWSIM,LHE \
    --customise Configuration/DataProcessing/Utils.addMonitoring \
    --datatier GEN-SIM,LHE \
    --conditions 124X_mcRun3_2022_realistic_v12 \
    --beamspot Realistic25ns13p6TeVEarly2022Collision \
    --customise_commands "$customise" \
    --step LHE,GEN,SIM \
    --geometry DB:Extended \
    --era Run3 \
    --python_filename "${TAG}__LHE__cfg.py" \
    --fileout "file:$SIM_FILE" \
    --number "$NEVENTS" \
    --nThreads "$NTHREADS" \
    --no_exec --mc
  cmsRun -j FrameworkJobReport.xml "${TAG}__LHE__cfg.py"
}

run_digi_reco() {
  echo "DIGI/HLT and RECO: CMSSW_12_4_11_patch3"
  # Central recipe: HIG-Run3Summer22DRPremix-00103.
  # Keep the central premix INPUT dataset; output transfer is handled by CRAB.
  cmsDriver.py \
    --eventcontent PREMIXRAW \
    --customise Configuration/DataProcessing/Utils.addMonitoring \
    --datatier GEN-SIM-RAW \
    --conditions 124X_mcRun3_2022_realistic_v12 \
    --step DIGI,DATAMIX,L1,DIGI2RAW,HLT:2022v12 \
    --procModifiers premix_stage2,siPixelQualityRawToDigi \
    --geometry DB:Extended \
    --datamix PreMix \
    --era Run3 \
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

  cmsDriver.py \
    --eventcontent AODSIM \
    --customise Configuration/DataProcessing/Utils.addMonitoring \
    --datatier AODSIM \
    --conditions 124X_mcRun3_2022_realistic_v12 \
    --step RAW2DIGI,L1Reco,RECO,RECOSIM \
    --procModifiers siPixelQualityRawToDigi \
    --geometry DB:Extended \
    --era Run3 \
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
  echo "MiniAODv4 and NanoAODv12: CMSSW_13_0_13"
  # Central recipe: HIG-Run3Summer22MiniAODv4-00116.
  cmsDriver.py \
    --eventcontent MINIAODSIM \
    --customise Configuration/DataProcessing/Utils.addMonitoring \
    --datatier MINIAODSIM \
    --conditions 130X_mcRun3_2022_realistic_v5 \
    --step PAT \
    --geometry DB:Extended \
    --era Run3,run3_miniAOD_12X \
    --python_filename "${TAG}__MiniAODv4__cfg.py" \
    --fileout "file:$MINIAOD_FILE" \
    --filein "file:$AOD_FILE" \
    --number -1 \
    --nThreads "$NTHREADS" \
    --no_exec --mc
  cmsRun -j FrameworkJobReport.xml "${TAG}__MiniAODv4__cfg.py"
  rm -f -- "$AOD_FILE"

  # Central recipe: HIG-Run3Summer22NanoAODv12-00116.
  # Use the requested flat NanoAOD output format.
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
  # The report describes the actual final file, with no subsequent renaming.
  cmsRun -j FrameworkJobReport.xml "${TAG}__NanoAODv12__cfg.py"
  [[ -s "$NANOAOD_FILE" ]] || { echo "ERROR: missing NanoAOD output" >&2; exit 20; }
  root -l -b -q -e "TFile f(\"$NANOAOD_FILE\"); if(f.IsZombie()){gSystem->Exit(21);} auto t=(TTree*)f.Get(\"Events\"); if(!t || t->GetEntries()==0){gSystem->Exit(21);} std::cout << \"Events = \" << t->GetEntries() << std::endl; gSystem->Exit(0);"
  rm -f -- "$MINIAOD_FILE"
}

run_in_cmssw CMSSW_12_4_11_patch3 el8_amd64_gcc10 run_lhe_gen_sim
run_in_cmssw CMSSW_12_4_11_patch3 el8_amd64_gcc10 run_digi_reco
run_in_cmssw CMSSW_13_0_13 el8_amd64_gcc11 run_mini_nano
echo "Production completed: $NANOAOD_FILE"
