# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: step1 --filein /store/user/dsperka/Run2MC/TunePowheg/submission_LHE_ggH_JHU_125_Apr22_hfact0p5/cmsgrid_final_1.lhe --fileout file:ggH_125_hfact0p5_LHE.root --mc --eventcontent LHE --datatier GEN --conditions MCRUN2_71_V1::All --step NONE --python_filename ggH_125_hfact0p5_LHEtoEDM_cfg.py --no_exec --customise Configuration/DataProcessing/Utils.addMonitoring -n -1
import FWCore.ParameterSet.Config as cms
import sys
import os

try:
    if not sys.argv[2].endswith(".root") or os.path.exists(sys.argv[2]):
        raise ValueError("First argument {} should end with .root and not exist".format(sys.argv[2]))

    for filename in sys.argv[3:]:
        if not filename.endswith(".lhe") or not os.path.exists(filename):
            raise ValueError("Second argument and further {} should end with .lhe and exist".format(filename))
except:
    print sys.argv
    print "cmsRun", sys.argv[1], "outputfile.root inputfile1.lhe inputfile2.lhe ..."
    raise

process = cms.Process('LHE')

# import of standard configurations
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.MessageLogger.cerr.FwkReport.reportEvery = 5000

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Input source
process.source = cms.Source("LHESource",
    fileNames = cms.untracked.vstring(
        *('file:'+_ for _ in sys.argv[3:])
    )
)   

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.19 $'),
    annotation = cms.untracked.string('step1 nevts:-1'),
    name = cms.untracked.string('Applications')
)

# Output definition

process.LHEoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.LHEEventContent.outputCommands,
    fileName = cms.untracked.string('file:{}'.format(sys.argv[2])),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN')
    )
)

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'MCRUN2_71_V1', '')

# Path and EndPath definitions
process.LHEoutput_step = cms.EndPath(process.LHEoutput)

# Schedule definition
process.schedule = cms.Schedule(process.LHEoutput_step)

# customisation of the process.

# Automatic addition of the customisation function from Configuration.DataProcessing.Utils
from Configuration.DataProcessing.Utils import addMonitoring 

#call to customisation function addMonitoring imported from Configuration.DataProcessing.Utils
process = addMonitoring(process)

# End of customisation functions
