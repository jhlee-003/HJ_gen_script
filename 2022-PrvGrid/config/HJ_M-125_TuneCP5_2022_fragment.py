import os

import FWCore.ParameterSet.Config as cms


# Set by the producer in the worker sandbox, not a developer-machine path.
# An absolute path also works in ExternalLHEProducer's concurrent thread dirs.
_gridpack_path = os.environ.get("GRIDPACK_PATH", "")
if not os.path.isabs(_gridpack_path):
    raise RuntimeError("GRIDPACK_PATH must be absolute; run the production script to set it.")
if not os.path.isfile(_gridpack_path) or not os.access(_gridpack_path, os.R_OK):
    raise RuntimeError("Custom gridpack is missing or unreadable: " + _gridpack_path)


externalLHEProducer = cms.EDProducer(
    "ExternalLHEProducer",
    args=cms.vstring(_gridpack_path),
    nEvents=cms.untracked.uint32(5000),
    numberOfParameters=cms.uint32(1),
    outputFile=cms.string("cmsgrid_final.lhe"),
    generateConcurrently=cms.untracked.bool(True),
    scriptName=cms.FileInPath(
        "GeneratorInterface/LHEInterface/data/run_generic_tarball_cvmfs.sh"
    ),
)


from Configuration.Generator.Pythia8CommonSettings_cfi import *
from Configuration.Generator.MCTunesRun3ECM13p6TeV.PythiaCP5Settings_cfi import *
from Configuration.Generator.PSweightsPythia.PythiaPSweightsSettings_cfi import *


generator = cms.EDFilter(
    "Pythia8ConcurrentHadronizerFilter",
    maxEventsToPrint=cms.untracked.int32(1),
    pythiaPylistVerbosity=cms.untracked.int32(1),
    filterEfficiency=cms.untracked.double(1.0),
    pythiaHepMCVerbosity=cms.untracked.bool(False),
    comEnergy=cms.double(13600.0),
    PythiaParameters=cms.PSet(
        pythia8CommonSettingsBlock,
        pythia8CP5SettingsBlock,
        pythia8PSweightsSettingsBlock,
        processParameters=cms.vstring(
            # HJ MiNNLO shower-matching prescription
            "SpaceShower:pTmaxMatch = 1",
            "TimeShower:pTmaxMatch = 1",
            # mH = 125 GeV and H -> Z gamma
            "25:m0 = 125.0",
            "25:onMode = off",
            "25:onIfMatch = 23 22",
            # Z -> ee, mumu, or tautau with mZ > 50 GeV
            "23:mMin = 50.0",
            "23:onMode = off",
            "23:onIfAny = 11 13 15",
        ),
        parameterSets=cms.vstring(
            "pythia8CommonSettings",
            "pythia8CP5Settings",
            "pythia8PSweightsSettings",
            "processParameters",
        ),
    ),
)


ProductionFilterSequence = cms.Sequence(generator)
