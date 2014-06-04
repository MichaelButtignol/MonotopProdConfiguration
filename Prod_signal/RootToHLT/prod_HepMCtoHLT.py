# Auto generated configuration file
# using: 
# Revision: 1.381.2.27 
# Source: /local/reps/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: step0 --filein file:Monotop_S1_mres1000_mchi100_HepMC.root --fileout file:Monotop_S1_mres1000_mchi100_HLT.root --mc --eventcontent RECOSIM --datatier GEN-SIM-RAW --conditions START53_V7C::All --pileup 2012_Summer_50ns_PoissonOOTPU --step GEN,SIM,DIGI,L1,DIGI2RAW,HLT:7E33v2 --python_filename prod_HepMCtoHLT.py --no_exec -n 10000
import FWCore.ParameterSet.Config as cms

process = cms.Process('HLT')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mix_2012_Summer_50ns_PoissonOOTPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.GeometrySimDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedRealistic8TeVCollision_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.DigiToRaw_cff')
process.load('HLTrigger.Configuration.HLT_7E33v2_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

#process.MessageLogger = cms.Service("MessageLogger",
#     statistics = cms.untracked.vstring('cout'),
#     cout = cms.untracked.PSet(
#         threshold = cms.untracked.string('DEBUG'),
#         noLineBreaks = cms.untracked.bool(True),
#         INFO = cms.untracked.PSet(
#             limit = cms.untracked.int32(0)
#         ),
#         DEBUG = cms.untracked.PSet(
#             limit = cms.untracked.int32(0)
#         ),
#         WARNING = cms.untracked.PSet(
#             limit = cms.untracked.int32(0)
#         ),
#         ERROR = cms.untracked.PSet(
#             limit = cms.untracked.int32(0)
#         ),
#     ),
#     
#     destinations = cms.untracked.vstring('cout',  ## .log automatically
#         )
# )

# Input source
process.source = cms.Source("PoolSource",
    secondaryFileNames= cms.untracked.vstring(),
     fileNames = cms.untracked.vstring('/store/user/mbuttign/Prod_S1/prod_S1_mres1000p0_mchi825p0/prod_S1_mres1000p0_mchi825p0_1_HepMC.root')

)


process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.381.2.27 $'),
    annotation = cms.untracked.string('step0 nevts:-1'),
    name = cms.untracked.string('PyReleaseValidation')
)

# Output definition

process.RECOSIMoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.RECOSIMEventContent.outputCommands,
    fileName = cms.untracked.string('file:prod_S1_mres1000p0_mchi825p0_1_HLT_n1000.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM-RAW')
    ),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    )
)

# Additional output definition

# Other statements
process.genstepfilter.triggerConditions=cms.vstring("generation_step")
from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'START53_V7C::All', '')
process.GlobalTag = GlobalTag(process.GlobalTag, 'START53_V26::All', '')

# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.simulation_step = cms.Path(process.psim)
process.digitisation_step = cms.Path(process.pdigi)
process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.digi2raw_step = cms.Path(process.DigiToRaw)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RECOSIMoutput_step = cms.EndPath(process.RECOSIMoutput)

process.g4SimHits.Generator.HepMCProductLabel = 'source'
process.VtxSmearedCommon.src = 'source'
process.genParticles.src = 'source'
#process.BetafuncEvtVtxGenerator/'VtxSmeared'
process.VtxSmeared.src='source'

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.simulation_step,process.digitisation_step,process.L1simulation_step,process.digi2raw_step)
process.schedule.extend(process.HLTSchedule)
process.schedule.extend([process.endjob_step,process.RECOSIMoutput_step])

# customisation of the process.

# Automatic addition of the customisation function from HLTrigger.Configuration.customizeHLTforMC
from HLTrigger.Configuration.customizeHLTforMC import customizeHLTforMC 

#call to customisation function customizeHLTforMC imported from HLTrigger.Configuration.customizeHLTforMC
process = customizeHLTforMC(process)

# End of customisation functions
