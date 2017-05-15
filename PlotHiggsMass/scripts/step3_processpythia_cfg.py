#!/usr/bin/env cmsRun

import FWCore.ParameterSet.Config as cms
import os
import sys

argv = sys.argv[2:]
try:
    outfile = sys.argv[2]
    if not outfile.endswith(".root") or os.path.exists(outfile):
        raise ValueError("First argument {} should end with .root and not exist".format(sys.argv[2]))

    VBForVH = sys.argv[3]
    if VBForVH in ("ZH", "WH"):
        VBForVH = "VH"
    if VBForVH not in ("VBF", "VH", "ggH", "qqZZ"):
        raise ValueError("Second argument should be VBF, VH, ggH, or qqZZ")

    infiles = sys.argv[4:]
    if not infiles:
        raise ValueError("No input files!")
    for filename in infiles:
        if not filename.endswith(".root") or not os.path.exists(filename):
            raise ValueError("Third argument and further {} should end with .root and exist".format(filename))
except:
    print sys.argv
    print "cmsRun", sys.argv[1], "outputfile.root VBForVHorggHorqqZZ inputfile1.root inputfile2.root ..."
    raise

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.load("RecoJets.Configuration.GenJetParticles_cff")
process.load("RecoJets.JetProducers.ak5GenJets_cfi")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.TFileService = cms.Service("TFileService", 
    fileName = cms.string("file:"+outfile),
    closeFileFast = cms.untracked.bool(True)
)

process.demo = cms.EDAnalyzer('PlotHiggsMass',
    VBForVH = cms.string(VBForVH),
    smearelectronpt = cms.double(2.399/6),
    smearelectroneta = cms.double(0),
    smearelectronphi = cms.double(0),
    smearmuonpt = cms.double(2.169/6),
    smearmuoneta = cms.double(0),
    smearmuonphi = cms.double(0),
    smearjetpt = cms.double(18./6),
    smearjeteta = cms.double(0),
    smearjetphi = cms.double(0),
    randomseed = cms.uint32(hash(tuple(sorted(infiles))) % 0xFFFFFFFF),  #random, but deterministic for given infiles
    keepallevents = cms.bool(True),
)


process.p = cms.Path(process.genJetParticles + process.ak5GenJets + process.demo)


process.schedule = cms.Schedule(
    process.p
)


myfilelist = cms.untracked.vstring()
myfilelist.extend( [
        "file:"+infile for infile in infiles
])

process.source = cms.Source("PoolSource",
    fileNames = myfilelist,
    duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
)
