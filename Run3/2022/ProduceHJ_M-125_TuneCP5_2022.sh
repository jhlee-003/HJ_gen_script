#!/bin/bash
set -euo pipefail

if [ $# -ne 3 ]; then
  echo "[Usage] $0 JOB_NUMBER NUMBER_OF_EVENTS ENV_FILE" >&2
  echo "  JOB_NUMBER starts from 0; output file numbering starts from 1." >&2
  echo "  ENV_FILE sets Fragment_filename, AOD_NAME, MINIAOD_NAME, NANOAOD_NAME, and BASE_TAG." >&2
  echo "  NTHREADS can be set in the environment and defaults to 2." >&2
  exit 2
fi

if [[ ! "$1" =~ ^[0-9]+$ ]]; then
  echo "ERROR: JOB_NUMBER must be a non-negative integer" >&2
  exit 2
fi
if [[ ! "$2" =~ ^[1-9][0-9]*$ ]]; then
  echo "ERROR: NUMBER_OF_EVENTS must be a positive integer" >&2
  exit 2
fi
if [ ! -f "$3" ]; then
  echo "ERROR: ENV_FILE does not exist: $3" >&2
  exit 2
fi

source "$3"
: "${Fragment_filename:?}" "${AOD_NAME:?}" "${MINIAOD_NAME:?}" "${NANOAOD_NAME:?}" "${BASE_TAG:?}"
NTHREADS=${NTHREADS:-2}
[[ "$NTHREADS" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: NTHREADS must be a positive integer" >&2; exit 2; }

echo "Set below variables with $3"
echo "Fragment_filename = $Fragment_filename"
echo "AOD_NAME = $AOD_NAME"
echo "MINIAOD_NAME = $MINIAOD_NAME"
echo "NANOAOD_NAME = $NANOAOD_NAME"
echo "BASE_TAG = $BASE_TAG"
echo "NTHREADS = $NTHREADS"

# Set variables. The scheduler index starts from zero.
JOBNUM=$(($1 + 1))
NEVENTS=$2
TAG="$BASE_TAG""__job-"${JOBNUM}

# CRAB transfers the fragment into the worker directory. A local checkout
# already keeps it under config/. Normalize both layouts to the DY structure.
mkdir -p config
[[ -e "$Fragment_filename" ]] && mv -- "$Fragment_filename" config/

if [ ! -f "config/${Fragment_filename}" ]; then
  echo "ERROR: config/${Fragment_filename} does not exist" >&2
  exit 3
fi

# The proxy is inherited from CRAB or HTCondor; no repository-local proxy file
# is created or required here.
mkdir -p job_scripts

cat <<EndOfTestFile > job_scripts/"$TAG"_cmd.sh
#!/bin/bash
set -euo pipefail


echo "--------------------------------LHE,GEN,SIM-------------------------------------"
echo "Setting up CMSSW"
export SCRAM_ARCH=el8_amd64_gcc10
[[ \$- == *u* ]] && U_WAS_ON=1 || U_WAS_ON=0; set +u; source /cvmfs/cms.cern.ch/cmsset_default.sh; ((U_WAS_ON)) && set -u
if [ -r CMSSW_12_4_11_patch3/src ] ; then
  echo release CMSSW_12_4_11_patch3 already exists
else
  scram p CMSSW CMSSW_12_4_11_patch3
fi
cd CMSSW_12_4_11_patch3/src
eval \`scram runtime -sh\`

# Setup custom fragment for CMSSW
mkdir -p Configuration/GenProduction/python
cp ../../config/${Fragment_filename} Configuration/GenProduction/python
scram b -j ${NTHREADS}
cd ../..


echo "Make cmssw configuration file"
Output_filename=$AOD_NAME"__job-"${JOBNUM}"__SIM".root
cmsDriver.py Configuration/GenProduction/python/$Fragment_filename \
  --era Run3 \
  --customise Configuration/DataProcessing/Utils.addMonitoring \
  --beamspot Realistic25ns13p6TeVEarly2022Collision \
  --step LHE,GEN,SIM \
  --geometry DB:Extended \
  --conditions 124X_mcRun3_2022_realistic_v12 \
  --customise_commands \
    "process.source.numberEventsInLuminosityBlock=cms.untracked.uint32(200); process.source.firstRun=cms.untracked.uint32(${JOBNUM}); from IOMC.RandomEngine.RandomServiceHelper import RandomNumberServiceHelper; randSvc = RandomNumberServiceHelper(process.RandomNumberGeneratorService); randSvc.populate()" \
  --datatier GEN-SIM,LHE \
  --eventcontent RAWSIM,LHE \
  --python_filename "$TAG"__LHE__cfg.py \
  --fileout file:\$Output_filename \
  -n $NEVENTS \
  --nThreads ${NTHREADS} \
  --no_exec \
  --mc

echo "Run cmssw with configuration file"
cmsRun -j FrameworkJobReport.xml "$TAG"__LHE__cfg.py


echo "---------------------------------DIGIPREMIX------------------------------------"
echo "Setting up CMSSW"
export SCRAM_ARCH=el8_amd64_gcc10
[[ \$- == *u* ]] && U_WAS_ON=1 || U_WAS_ON=0; set +u; source /cvmfs/cms.cern.ch/cmsset_default.sh; ((U_WAS_ON)) && set -u
if [ -r CMSSW_12_4_11_patch3/src ] ; then
  echo release CMSSW_12_4_11_patch3 already exists
else
  scram p CMSSW CMSSW_12_4_11_patch3
fi
cd CMSSW_12_4_11_patch3/src
eval \`scram runtime -sh\`
scram b -j ${NTHREADS}
cd ../..


echo "Make cmssw configuration file"
Input_filename=$AOD_NAME"__job-"${JOBNUM}"__SIM".root
Output_filename=$AOD_NAME"__job-"${JOBNUM}"__HLT".root
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
  --python_filename "$TAG"__DIGIPREMIX__cfg.py \
  --fileout file:\$Output_filename \
  --filein file:\$Input_filename \
  -n -1 \
  --nThreads ${NTHREADS} \
  --pileup_input "dbs:/Neutrino_E-10_gun/Run3Summer21PrePremix-Summer22_124X_mcRun3_2022_realistic_v11-v2/PREMIX" \
  --no_exec \
  --mc

echo "Run cmssw with configuration file"
cmsRun -j FrameworkJobReport.xml "$TAG"__DIGIPREMIX__cfg.py

echo "Clean up files"
rm -f \$Input_filename ${AOD_NAME}__job-${JOBNUM}__SIM_inLHE.root


echo "------------------------------------RECO----------------------------------------"
echo "Make cmssw configuration file"
Input_filename=$AOD_NAME"__job-"${JOBNUM}"__HLT".root
Output_filename=$AOD_NAME"__job-"${JOBNUM}.root
cmsDriver.py \
  --era Run3 \
  --customise Configuration/DataProcessing/Utils.addMonitoring \
  --procModifiers siPixelQualityRawToDigi \
  --step RAW2DIGI,L1Reco,RECO,RECOSIM \
  --geometry DB:Extended \
  --conditions 124X_mcRun3_2022_realistic_v12 \
  --datatier AODSIM \
  --eventcontent AODSIM \
  --python_filename "$TAG"__AOD__cfg.py \
  --fileout file:\$Output_filename \
  --filein file:\$Input_filename \
  -n -1 \
  --nThreads ${NTHREADS} \
  --no_exec \
  --mc

echo "Run cmssw with configuration file"
cmsRun -j FrameworkJobReport.xml "$TAG"__AOD__cfg.py

echo "Clean up files"
rm -f \$Input_filename


echo "----------------------------------MiniAODv4-------------------------------------"
echo "Setting up CMSSW"
export SCRAM_ARCH=el8_amd64_gcc11
[[ \$- == *u* ]] && U_WAS_ON=1 || U_WAS_ON=0; set +u; source /cvmfs/cms.cern.ch/cmsset_default.sh; ((U_WAS_ON)) && set -u
if [ -r CMSSW_13_0_13/src ] ; then
  echo release CMSSW_13_0_13 already exists
else
  scram p CMSSW CMSSW_13_0_13
fi
cd CMSSW_13_0_13/src
eval \`scram runtime -sh\`
scram b -j ${NTHREADS}
cd ../..


echo "Make cmssw configuration file"
Input_filename=$AOD_NAME"__job-"${JOBNUM}.root
Output_filename=$MINIAOD_NAME"__job-"${JOBNUM}.root
cmsDriver.py \
  --era Run3,run3_miniAOD_12X \
  --customise Configuration/DataProcessing/Utils.addMonitoring \
  --step PAT \
  --geometry DB:Extended \
  --conditions 130X_mcRun3_2022_realistic_v5 \
  --datatier MINIAODSIM \
  --eventcontent MINIAODSIM \
  --python_filename "$TAG"__MiniAODv4__cfg.py \
  --fileout file:\$Output_filename \
  --filein file:\$Input_filename \
  -n -1 \
  --nThreads ${NTHREADS} \
  --no_exec \
  --mc

echo "Run cmssw with configuration file"
cmsRun -j FrameworkJobReport.xml "$TAG"__MiniAODv4__cfg.py

echo "Clean up files"
rm -f \$Input_filename


echo "-----------------------------------NanoAODv12------------------------------------"
echo "Setting up CMSSW"
export SCRAM_ARCH=el8_amd64_gcc11
[[ \$- == *u* ]] && U_WAS_ON=1 || U_WAS_ON=0; set +u; source /cvmfs/cms.cern.ch/cmsset_default.sh; ((U_WAS_ON)) && set -u
if [ -r CMSSW_13_0_13/src ] ; then
  echo release CMSSW_13_0_13 already exists
else
  scram p CMSSW CMSSW_13_0_13
fi
cd CMSSW_13_0_13/src
eval \`scram runtime -sh\`
scram b -j ${NTHREADS}
cd ../..


echo "Make cmssw configuration file"
Input_filename=$MINIAOD_NAME"__job-"${JOBNUM}.root
Output_filename=$NANOAOD_NAME"__job-"${JOBNUM}.root
cmsDriver.py \
  --scenario pp \
  --era Run3 \
  --customise Configuration/DataProcessing/Utils.addMonitoring \
  --step NANO \
  --conditions 130X_mcRun3_2022_realistic_v5 \
  --datatier NANOAODSIM \
  --eventcontent NANOAODSIM \
  --python_filename "$TAG"__NanoAODv12__cfg.py \
  --fileout file:\$Output_filename \
  --filein file:\$Input_filename \
  -n -1 \
  --nThreads ${NTHREADS} \
  --no_exec \
  --mc

echo "Run cmssw with configuration file"
cmsRun -j FrameworkJobReport.xml "$TAG"__NanoAODv12__cfg.py

# Validate the final output before removing the MiniAOD input.
if [ ! -s \$Output_filename ]; then
  echo "ERROR: missing NanoAOD output: \$Output_filename" >&2
  exit 20
fi
root -l -b -q -e "TFile f(\"\$Output_filename\"); if(f.IsZombie()){gSystem->Exit(21);} auto t=(TTree*)f.Get(\"Events\"); if(!t || t->GetEntries()==0){gSystem->Exit(21);} std::cout << \"Events = \" << t->GetEntries() << std::endl; gSystem->Exit(0);"


#-------CleanUp---------

echo "Clean up files"
rm -f \$Input_filename
rm -f ${TAG}__LHE__cfg.py
rm -f ${TAG}__DIGIPREMIX__cfg.py
rm -f ${TAG}__AOD__cfg.py
rm -f ${TAG}__MiniAODv4__cfg.py
rm -f ${TAG}__NanoAODv12__cfg.py

echo "Production completed: \$Output_filename"

# End of "$TAG"_cmd.sh file
EndOfTestFile

echo "Made ${TAG}_cmd.sh"
chmod +x job_scripts/"$TAG"_cmd.sh

./job_scripts/${TAG}_cmd.sh
