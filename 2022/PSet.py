import FWCore.ParameterSet.Config as cms

# CRAB PrivateMC bootstrap. The wrapper runs the full production chain.
process = cms.Process("CRABBOOTSTRAP")
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(1))
process.source = cms.Source("EmptySource")
process.options = cms.untracked.PSet(
    numberOfThreads=cms.untracked.uint32(2),
    numberOfStreams=cms.untracked.uint32(0),
)
