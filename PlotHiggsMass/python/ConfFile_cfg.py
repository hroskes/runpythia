import FWCore.ParameterSet.Config as cms
import sys

argv = sys.argv[2:]
try:
    infile = argv[0]
    firstevent = int(argv[1])
    lastevent = int(argv[2])
    outfile = argv[3]
except:
    print sys.argv
    print "cmsRun", sys.argv[1], "infile firstevent lastevent outfile"
    raise

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("RecoJets.Configuration.GenJetParticles_cff")
process.load("RecoJets.JetProducers.ak5GenJets_cfi")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(lastevent - firstevent + 1)
)

process.TFileService = cms.Service("TFileService", 
    fileName = cms.string(outfile),
    closeFileFast = cms.untracked.bool(True)
)

process.demo = cms.EDAnalyzer('PlotHiggsMass'
)


process.p = cms.Path(process.genJetParticles + process.ak5GenJets + process.demo)


process.schedule = cms.Schedule(
    process.p
)


myfilelist = cms.untracked.vstring()
myfilelist.extend( [
        infile,
])

process.source = cms.Source("PoolSource",
    fileNames = myfilelist,
    duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
)
