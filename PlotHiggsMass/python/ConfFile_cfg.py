import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("RecoJets.Configuration.GenJetParticles_cff")
process.load("RecoJets.JetProducers.ak5GenJets_cfi")

#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(15000) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100000) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.TFileService = cms.Service("TFileService", 
#    fileName = cms.string("Tree_WH_125_JHUgen.root"),
#    fileName = cms.string("Tree_WminusH_125_Powheg.root"),
#    fileName = cms.string("Tree_WplusH_125_Powheg.root"),
#    fileName = cms.string("Tree_VBF_125_phantom.root"),
#    fileName = cms.string("Tree_VBF_125_phantom_pTMaxMatch.root"),
#    fileName = cms.string("Tree_VBF_125_Powheg.root"),
#    fileName = cms.string("Tree_ttH_125_JHUgen.root"),
#    fileName = cms.string("Tree_ttH_125_JHUgen_pTmaxmatch.root"),
#    fileName = cms.string("Tree_ttH_125_Powheg.root"),
#    fileName = cms.string("Tree_ggH_750_Powheg.root"),
#    fileName = cms.string("Tree_ggH_125_Powheg.root"),
#    fileName = cms.string("Tree_ggH2bplus_750_JHU.root"),
#    fileName = cms.string("Tree_ggH2bplus_750_Fudge0p5_JHU.root"),
#    fileName = cms.string("Tree_ggH2bplus_750_Fudge0p5_Match1_JHU.root"),
#    fileName = cms.string("Tree_ggH2bplus_750_Fudge0p25_Match1_JHU.root"),
#    fileName = cms.string("Tree_ggH2bplus_125_JHU.root"),
    fileName = cms.string("Tree_ggH2bplus_125_JHU_fixscale_time.root"),
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
        #'file:WH_125_JHUgen_GEN.root'
        #'/store/mc/RunIISpring15DR74/WminusH_HToZZTo4L_M125_13TeV_powheg-minlo-HWJ_JHUgen_pythia8/AODSIM/Asympt25ns_MCRUN2_74_V9-v2/70000/E8B31FD4-900C-E511-BEF7-848F69FD46E8.root'
        #'/store/mc/RunIISpring15DR74/WplusH_HToZZTo4L_M125_13TeV_powheg-minlo-HWJ_JHUgen_pythia8/AODSIM/Asympt25ns_MCRUN2_74_V9-v1/70000/FA7F53FE-A506-E511-8078-000F530E4BD4.root'
        #'file:../../CMSSW_7_1_20_patch2/src/VBF_125_phantom_GEN.root'
        #'file:../../CMSSW_7_1_20_patch2/src/VBF_125_phantom_GEN_pTMaxMatch.root'
        #'/store/test/xrootd/T2_US_Purdue//store/mc/RunIISpring15DR74/VBF_HToZZTo4L_M125_13TeV_powheg_JHUgen_pythia8/AODSIM/Asympt25ns_MCRUN2_74_V9-v1/00000/16584C9C-2713-E511-A640-0CC47A13D216.root',
        #'/store/test/xrootd/T2_US_Purdue//store/mc/RunIISpring15DR74/VBF_HToZZTo4L_M125_13TeV_powheg_JHUgen_pythia8/AODSIM/Asympt25ns_MCRUN2_74_V9-v1/00000/88B85EDA-2713-E511-8795-0CC47A01035C.root',
        #'/store/test/xrootd/T2_US_Purdue//store/mc/RunIISpring15DR74/VBF_HToZZTo4L_M125_13TeV_powheg_JHUgen_pythia8/AODSIM/Asympt25ns_MCRUN2_74_V9-v1/00000/9A2EDBC1-2713-E511-BA83-001EC9B0AD32.root',
        #'/store/test/xrootd/T2_US_Purdue//store/mc/RunIISpring15DR74/VBF_HToZZTo4L_M125_13TeV_powheg_JHUgen_pythia8/AODSIM/Asympt25ns_MCRUN2_74_V9-v1/00000/F2B14B71-2713-E511-8A76-002590A37128.root',
        #'/store/test/xrootd/T2_US_Purdue//store/mc/RunIISpring15DR74/VBF_HToZZTo4L_M125_13TeV_powheg_JHUgen_pythia8/AODSIM/Asympt25ns_MCRUN2_74_V9-v1/10000/14DFD077-E312-E511-BB0E-6C3BE5B5F228.root',
        #'/store/test/xrootd/T2_US_Purdue//store/mc/RunIISpring15DR74/VBF_HToZZTo4L_M125_13TeV_powheg_JHUgen_pythia8/AODSIM/Asympt25ns_MCRUN2_74_V9-v1/10000/40206E0E-C012-E511-B1DE-0025905A6138.root',
        #'/store/test/xrootd/T2_US_Purdue//store/mc/RunIISpring15DR74/VBF_HToZZTo4L_M125_13TeV_powheg_JHUgen_pythia8/AODSIM/Asympt25ns_MCRUN2_74_V9-v1/10000/4A200BF4-E212-E511-B98B-0CC47A13CB36.root',
        #'/store/test/xrootd/T2_US_Purdue//store/mc/RunIISpring15DR74/VBF_HToZZTo4L_M125_13TeV_powheg_JHUgen_pythia8/AODSIM/Asympt25ns_MCRUN2_74_V9-v1/10000/54277375-E312-E511-920E-0025905C6448.root',
        #'/store/test/xrootd/T2_US_Purdue//store/mc/RunIISpring15DR74/VBF_HToZZTo4L_M125_13TeV_powheg_JHUgen_pythia8/AODSIM/Asympt25ns_MCRUN2_74_V9-v1/10000/A853A1C6-CF12-E511-AC09-0025905A610A.root',
        #'/store/test/xrootd/T2_US_Purdue//store/mc/RunIISpring15DR74/VBF_HToZZTo4L_M125_13TeV_powheg_JHUgen_pythia8/AODSIM/Asympt25ns_MCRUN2_74_V9-v1/10000/E67AB1D6-E012-E511-8A98-003048344A94.root',       
        #'/store/test/xrootd/T2_US_Purdue//store/mc/RunIISpring15DR74/VBF_HToZZTo4L_M125_13TeV_powheg_JHUgen_pythia8/AODSIM/Asympt25ns_MCRUN2_74_V9-v1/10000/F4628D42-DA12-E511-9D8E-0CC47A13D09C.root',
        #'/store/test/xrootd/T2_US_Purdue//store/mc/RunIISpring15DR74/VBF_HToZZTo4L_M125_13TeV_powheg_JHUgen_pythia8/AODSIM/Asympt25ns_MCRUN2_74_V9-v1/50000/34F9B89D-AA12-E511-B8A7-0CC47A13D110.root',
        #'/store/test/xrootd/T2_US_Purdue//store/mc/RunIISpring15DR74/VBF_HToZZTo4L_M125_13TeV_powheg_JHUgen_pythia8/AODSIM/Asympt25ns_MCRUN2_74_V9-v1/50000/3C69E7BC-4912-E511-8FB1-0025904A96BC.root',
        #'/store/test/xrootd/T2_US_Purdue//store/mc/RunIISpring15DR74/VBF_HToZZTo4L_M125_13TeV_powheg_JHUgen_pythia8/AODSIM/Asympt25ns_MCRUN2_74_V9-v1/50000/60E9C84F-4012-E511-A405-0CC47A13CDA0.root',
        #'/store/test/xrootd/T2_US_Purdue//store/mc/RunIISpring15DR74/VBF_HToZZTo4L_M125_13TeV_powheg_JHUgen_pythia8/AODSIM/Asympt25ns_MCRUN2_74_V9-v1/50000/90CAF271-AB12-E511-9061-00259055C92C.root',
        #'/store/test/xrootd/T2_US_Purdue//store/mc/RunIISpring15DR74/VBF_HToZZTo4L_M125_13TeV_powheg_JHUgen_pythia8/AODSIM/Asympt25ns_MCRUN2_74_V9-v1/50000/A08FDAE7-5012-E511-BE65-0025904A96BC.root',
        #'/store/test/xrootd/T2_US_Purdue//store/mc/RunIISpring15DR74/VBF_HToZZTo4L_M125_13TeV_powheg_JHUgen_pythia8/AODSIM/Asympt25ns_MCRUN2_74_V9-v1/50000/C48DDE73-3F12-E511-AD79-00304865C45A.root',
        #'/store/test/xrootd/T2_US_Purdue//store/mc/RunIISpring15DR74/VBF_HToZZTo4L_M125_13TeV_powheg_JHUgen_pythia8/AODSIM/Asympt25ns_MCRUN2_74_V9-v1/50000/EABC13FF-5112-E511-9B68-0CC47A13CEAC.root',
        #'/store/test/xrootd/T2_US_Purdue//store/mc/RunIISpring15DR74/VBF_HToZZTo4L_M125_13TeV_powheg_JHUgen_pythia8/AODSIM/Asympt25ns_MCRUN2_74_V9-v1/60000/3A9291DE-5A13-E511-8354-0026B92785E9.root'
        #'/store/mc/RunIISpring15DR74/VBF_HToZZTo4L_M125_13TeV_powheg_JHUgen_pythia8/AODSIM/Asympt25ns_MCRUN2_74_V9-v1/00000/9A2EDBC1-2713-E511-BA83-001EC9B0AD32.root',
        #'file:ttH_125_JHUgen_GEN.root'
        #'file:ttH_125_JHUgen_pTmaxmatch_GEN.root'
        #'/store/mc/RunIISpring15DR74/ttH_HToZZ_4LFilter_M125_13TeV_powheg_JHUgen_pythia8/AODSIM/Asympt25ns_MCRUN2_74_V9-v1/20000/4410D124-4034-E511-AF67-003048FFCB9E.root',
        #'/store/mc/RunIISpring15DR74/ttH_HToZZ_4LFilter_M125_13TeV_powheg_JHUgen_pythia8/AODSIM/Asympt25ns_MCRUN2_74_V9-v1/20000/C642E4A0-F530-E511-AFB5-0002C90C5A08.root',
        #'/store/mc/RunIISpring15DR74/ttH_HToZZ_4LFilter_M125_13TeV_powheg_JHUgen_pythia8/AODSIM/Asympt25ns_MCRUN2_74_V9-v1/50000/02B0EF99-162B-E511-9789-0026189437F0.root',
        #'/store/mc/RunIISpring15DR74/ttH_HToZZ_4LFilter_M125_13TeV_powheg_JHUgen_pythia8/AODSIM/Asympt25ns_MCRUN2_74_V9-v1/50000/0CA1D3E8-1A2B-E511-868E-0025905A609E.root',
        #'/store/mc/RunIISpring15DR74/ttH_HToZZ_4LFilter_M125_13TeV_powheg_JHUgen_pythia8/AODSIM/Asympt25ns_MCRUN2_74_V9-v1/50000/0E2ECC84-1A2B-E511-8063-842B2B75FFFD.root'
        #'file:/afs/cern.ch/work/d/dsperka/Run2MC/4B/CMSSW_7_1_20_patch2/src/ggH_750_GEN.root'
        #'file:/afs/cern.ch/work/d/dsperka/Run2MC/4B/CMSSW_7_1_20_patch2/src/ggH_125_GEN.root'
        #'file:/afs/cern.ch/work/d/dsperka/Run2MC/4B/CMSSW_7_1_20_patch2/src/ggH2bplus_750_GEN.root'
        #'file:/afs/cern.ch/work/d/dsperka/Run2MC/4B/CMSSW_7_1_20_patch2/src/ggH2bplus_750_Fudge0p5_GEN.root'
        #'file:/afs/cern.ch/work/d/dsperka/Run2MC/4B/CMSSW_7_1_20_patch2/src/ggH2bplus_750_Fudge0p5_Match1_GEN.root'
        #'file:/afs/cern.ch/work/d/dsperka/Run2MC/4B/CMSSW_7_1_20_patch2/src/ggH2bplus_750_Fudge0p25_Match1_GEN.root'
        #'file:/afs/cern.ch/work/d/dsperka/Run2MC/4B/CMSSW_7_1_20_patch2/src/ggH2bplus_125_GEN.root'
        #'file:/afs/cern.ch/work/d/dsperka/Run2MC/4B/CMSSW_7_1_20_patch2/src/ggH2bplus_125_GEN_fixscale.root'
        'file:/afs/cern.ch/work/d/dsperka/Run2MC/4B/CMSSW_7_1_20_patch2/src/ggH2bplus_125_GEN_fixscale_time.root'
])

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = myfilelist,
    duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
#    inputCommands = cms.untracked.vstring('drop LHEEventProduct_*_*_*',
#                                          'drop LHERunInfoProduct_*_*_*')
                            
)
