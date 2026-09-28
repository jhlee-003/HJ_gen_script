import os
import string

from WMCore.Configuration import Configuration

config = Configuration()

# Submit on KNU from an initialized CMSSW/CRAB environment:
# crab submit -c config/HJ_M-125_TuneCP5_2022_crabConfig.py
# CRAB supplies the worker proxy; do not ship a separate voms_proxy.txt.
project_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
env_filename = "HJ_M-125_TuneCP5_2022_FullSim.env"
script_filename = "ProduceHJ_M-125_TuneCP5_2022.sh"

# Read the shared naming file without executing shell code.
names = {}
with open(os.path.join(project_dir, "config", env_filename)) as env_file:
    for line in env_file:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        key, value = line.split("=", 1)
        names[key.strip()] = string.Template(value.strip().strip('"')).substitute(names)

base = names["BASE_TAG"]
work_area = os.path.join(project_dir, "crab_projects")

# Preserve the existing production size; tune events_per_job after a pilot.
# The complete LHE-to-NanoAOD chain must fit maxJobRuntimeMin.
events_per_job = 10000
number_of_jobs = 10000
num_cores = 2

# Numbered tasks avoid collisions without adding a timestamp to names.
task_number = 1
while os.path.exists(os.path.join(work_area, "crab_{}_{}".format(base, task_number))):
    task_number += 1
request_name = "{}_{}".format(base, task_number)

config.section_("General")
config.General.requestName = request_name
config.General.workArea = work_area
config.General.transferOutputs = True
config.General.transferLogs = False

config.section_("JobType")
config.JobType.pluginName = "PrivateMC"
config.JobType.psetName = os.path.join(project_dir, "PSet.py")
config.JobType.scriptExe = os.path.join(project_dir, "crab_convert_wrapper.sh")
config.JobType.numCores = num_cores
config.JobType.maxMemoryMB = 5000
config.JobType.maxJobRuntimeMin = 1200
config.JobType.eventsPerLumi = 200
config.JobType.inputFiles = [
    os.path.join(project_dir, "config", env_filename),
    os.path.join(project_dir, "config", names["Fragment_filename"]),
    os.path.join(project_dir, script_filename),
]
# The final cmsRun report and this explicit name refer to the same file.
config.JobType.disableAutomaticOutputCollection = True
config.JobType.outputFiles = [names["NANOAOD_NAME"] + "__job.root"]
config.JobType.scriptArgs = [
    "script=" + script_filename,
    "events=" + str(events_per_job),
    "names=" + env_filename,
    "threads=" + str(num_cores),
]

config.section_("Data")
config.Data.splitting = "EventBased"
config.Data.unitsPerJob = events_per_job
config.Data.totalUnits = events_per_job * number_of_jobs
config.Data.outputPrimaryDataset = base.split("__", 1)[0]
config.Data.outputDatasetTag = "Run3Summer22_NanoAODv12_privateProduction_{}".format(task_number)
config.Data.publication = False
# Use CRAB's default /store/user/<username>/ output base at the storage site.

config.section_("Site")
config.Site.storageSite = "T3_KR_KNU"

config.section_("Debug")
config.Debug.extraJDL = ["request_disk = 8000000"]
