# Starting with a skeleton process which gets imported with the following line
from PhysicsTools.PatAlgos.patTemplate_cfg import *

from PhysicsTools.PatAlgos.tools.coreTools import *

###############################
####### Parameters ############
###############################
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')

options.register ('useData',
                  False,
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.int,
                  'Run this on real data')

options.register ('globalTag',
                  '',
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.string,
                  'Overwrite defaul globalTag')


options.register ('forceCheckClosestZVertex',
                  False,
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.int,
                  "Force the check of the closest z vertex")


options.register ('runOnFastSim',
                  False,
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.int,
                  "Option needed to run on fastsim.")



print options.useData 
print options.useData 
print options.useData 
print options.useData 
print options.useData 
print options.useData 
print options.useData 
print options.useData 


if not options.useData :
    inputJetCorrLabel = ('AK5PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'])
    process.source.fileNames = [
        #'/store/mc/Summer12_DR53X/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V7C-v1/00000/FE7C71D8-DB25-E211-A93B-0025901D4C74.root' 
    #'file:000C5D15-AB1A-E211-8BDE-00215E22053A.root',
    #'file:0010005A-421A-E211-9E3C-E41F13181DA4.root'
#	'file:/opt/sbg/data/data4/cms/aaubin/synchronizationSample/64542608-69C7-E111-AFC8-002618943953.root',
#	'file:/opt/sbg/data/data4/cms/aaubin/synchronizationSample/8C41444F-62C7-E111-9FE7-002618943833.root',
#	'file:/opt/sbg/data/data4/cms/aaubin/synchronizationSample/D80E86FF-7EC7-E111-BB09-003048678DD6.root'
    #'file:ZjetsSynchro/70300E2E-27D2-E111-92BD-001E67397AE4.root',
    #'file:ZjetsSynchro/7041870B-D3D2-E111-8CFE-001E67397008.root',
    #'file:ZjetsSynchro/70B638CE-C6D2-E111-8430-003048673F0A.root',
    #'file:ZjetsSynchro/70CC5B25-C3D2-E111-85A7-001E6739751C.root',
    #'file:ZjetsSynchro/70EA5873-5AD2-E111-AE4F-003048D462C4.root'
    
    
    ]

else :
    inputJetCorrLabel = ('AK5PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual'])
    process.source.fileNames = [
        ''
    ]
#process.source.eventsToProcess = cms.untracked.VEventRange( ['1:100000000'] )

#process.source.skipEvents = cms.untracked.uint32(17268) 
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(500) )


print options

print 'Running jet corrections: '
print inputJetCorrLabel

import sys

#sys.path.append('s4_m500')
#import INPUT
#process.source.fileNames = INPUT.fileNames

##################################################################3



###############################
####### Global Setup ##########
###############################

if options.useData :
    if options.globalTag is '':
        process.GlobalTag.globaltag = cms.string( 'GR_P_V40_AN4::All' )
    else:
        process.GlobalTag.globaltag = cms.string( options.globalTag )
else :
    if options.globalTag is '':
        process.GlobalTag.globaltag = cms.string( 'START53_V26::All' ) #V20 should be used
        #process.GlobalTag.globaltag = cms.string( 'START53_V7G::All' ) #V20 should be used
    else:
        process.GlobalTag.globaltag = cms.string( options.globalTag )


from PhysicsTools.PatAlgos.patTemplate_cfg import *




###############################
########## PF Setup ###########
###############################

# Default PF2PAT with AK5 jets. Make sure to turn ON the L1fastjet stuff. 
from PhysicsTools.PatAlgos.tools.pfTools import *
postfix = "PFlow"
usePF2PAT(process,runPF2PAT=True, jetAlgo='AK5', runOnMC=not options.useData, postfix=postfix,
	  jetCorrections=inputJetCorrLabel, pvCollection=cms.InputTag('goodOfflinePrimaryVertices'), typeIMetCorrections=True)

#useGsfElectrons(process,postfix,dR="03")



####################################
#  top projections in PF2PAT:
####################################
getattr(process,"pfNoPileUp"  +postfix).enable = True
getattr(process,"pfNoMuon"    +postfix).enable = False 
getattr(process,"pfNoElectron"+postfix).enable = False
getattr(process,"pfNoTau"     +postfix).enable = False
getattr(process,"pfNoJet"     +postfix).enable = False





if not options.forceCheckClosestZVertex :
    process.pfPileUpPFlow.checkClosestZVertex = False

# change the cone size of electron isolation to 0.3 as default.
process.pfIsolatedElectronsPFlow.isolationValueMapsCharged = cms.VInputTag(cms.InputTag("elPFIsoValueCharged03PFIdPFlow"))
process.pfIsolatedElectronsPFlow.deltaBetaIsolationValueMap = cms.InputTag("elPFIsoValuePU03PFIdPFlow")
process.pfIsolatedElectronsPFlow.isolationValueMapsNeutral = cms.VInputTag(cms.InputTag("elPFIsoValueNeutral03PFIdPFlow"), cms.InputTag("elPFIsoValueGamma03PFIdPFlow"))

process.pfElectronsPFlow.isolationValueMapsCharged  = cms.VInputTag(cms.InputTag("elPFIsoValueCharged03PFIdPFlow"))
process.pfElectronsPFlow.deltaBetaIsolationValueMap = cms.InputTag("elPFIsoValuePU03PFIdPFlow" )
process.pfElectronsPFlow.isolationValueMapsNeutral  = cms.VInputTag(cms.InputTag( "elPFIsoValueNeutral03PFIdPFlow"), cms.InputTag("elPFIsoValueGamma03PFIdPFlow"))

process.patElectronsPFlow.isolationValues = cms.PSet(
        pfChargedHadrons = cms.InputTag("elPFIsoValueCharged03PFIdPFlow"),
        pfChargedAll = cms.InputTag("elPFIsoValueChargedAll03PFIdPFlow"),
        pfPUChargedHadrons = cms.InputTag("elPFIsoValuePU03PFIdPFlow"),
        pfNeutralHadrons = cms.InputTag("elPFIsoValueNeutral03PFIdPFlow"),
        pfPhotons = cms.InputTag("elPFIsoValueGamma03PFIdPFlow")
        )


applyPostfix(process,"pfIsolatedElectrons",postfix).combinedIsolationCut = cms.double(9999.)
 
 


# change the cone size of muons isolation to 0.3 as default.
applyPostfix(process,"pfIsolatedMuons",postfix).isolationValueMapsCharged = cms.VInputTag( cms.InputTag( 'muPFIsoValueCharged03PFlow' ) )
applyPostfix(process,"pfIsolatedMuons",postfix).isolationValueMapsNeutral = cms.VInputTag( cms.InputTag( 'muPFIsoValueNeutral03PFlow' ), cms.InputTag( 'muPFIsoValueGamma03PFlow' ) )
applyPostfix(process,"pfIsolatedMuons",postfix).deltaBetaIsolationValueMap = cms.InputTag( 'muPFIsoValuePU03PFlow' )
applyPostfix(process,"patMuons",postfix).isolationValues.pfNeutralHadrons = cms.InputTag( 'muPFIsoValueNeutral03PFlow' )
applyPostfix(process,"patMuons",postfix).isolationValues.pfPhotons = cms.InputTag( 'muPFIsoValueGamma03PFlow' )
applyPostfix(process,"patMuons",postfix).isolationValues.pfChargedHadrons = cms.InputTag( 'muPFIsoValueCharged03PFlow' )
applyPostfix(process,"patMuons",postfix).isolationValues.pfPUChargedHadrons = cms.InputTag( 'muPFIsoValuePU03PFlow' ) 
applyPostfix(process,"pfIsolatedMuons",postfix).combinedIsolationCut = cms.double(9999.)




applyPostfix(process,"pfIsolatedMuons",postfix).isolationCut = cms.double(9999.)
applyPostfix(process,"pfIsolatedElectrons",postfix).isolationCut = cms.double(9999.)   



# turn to false when running on data
if options.useData :
    removeMCMatching( process, ['All'] )

###############################
###### Electron ID ############
###############################

process.load('EgammaAnalysis.ElectronTools.electronIdMVAProducer_cfi') 
process.eidMVASequence = cms.Sequence(  process.mvaTrigV0 + process.mvaNonTrigV0 )
#Electron ID
process.patElectronsPFlow.electronIDSources.mvaTrigV0	 = cms.InputTag("mvaTrigV0")
process.patElectronsPFlow.electronIDSources.mvaNonTrigV0 = cms.InputTag("mvaNonTrigV0") 
process.patPF2PATSequencePFlow.replace( process.patElectronsPFlow, process.eidMVASequence * process.patElectronsPFlow )


#Convesion Rejection
# this should be your last selected electron collection name since currently index is used to match with electron later. We can fix this using reference pointer.
process.patConversionsPFlow = cms.EDProducer("PATConversionProducer",
                                             electronSource = cms.InputTag("selectedPatElectronsPFlow")      
                                             )
					     
process.patPF2PATSequencePFlow += process.patConversionsPFlow


switchJetCollection(process,cms.InputTag('ak5PFJets'),
		    doJTA        = False,
		    doBTagging   = True,
		    jetCorrLabel = inputJetCorrLabel,
		    doType1MET   = True,
		    genJetCollection=cms.InputTag("ak5GenJetsNoNu"),
		    doJetID      = True
		    )
		    
process.load("JetMETCorrections.Type1MET.pfMETCorrectionType0_cfi")
process.pfType1CorrectedMet.applyType0Corrections = cms.bool(True)
process.pfType1CorrectedMet.srcType1Corrections = cms.VInputTag(
    cms.InputTag('pfMETcorrType0'),
    cms.InputTag('pfJetMETcorr', 'type1')        
)		    



#################################################
## Produce trigger infos #
##################################################

process.load("PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cff")
from PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cfi import *
from PhysicsTools.PatAlgos.triggerLayer1.triggerEventProducer_cfi import *
process.patTrigger = patTrigger.clone()
process.patTriggerEvent = patTriggerEvent.clone()




#################################################
## Type I and Type 0 correction#
##################################################

#process.load("JetMETCorrections.Type1MET.pfMETCorrections_cff")




#process.load("JetMETCorrections.Type1MET.pfMETCorrections_cff")
#process.load("JetMETCorrections.Type1MET.pfMETCorrectionType0_cfi")
#process.pfType1CorrectedMet.applyType0Corrections = cms.bool(False)
#process.pfType1CorrectedMet.srcType1Corrections = cms.VInputTag(
#    cms.InputTag('pfMETcorrType0'),
#    cms.InputTag('pfJetMETcorr', 'type1')        
#)




###############################
####### Global Setup ##########
###############################

### The Strasbourg's ntupler

process.MiniTreeProduction = cms.EDProducer('MiniTreeProducer',
# ---------------------- General info -------------------------------
         verbose             = cms.uint32(0),   #0: nothing - >1 
         isAOD               = cms.bool(True),  # true if processing AOD 
         isData              = cms.bool(False), # true if processing AOD data
# ----------------------   Trigger -------------------------------
         doTrigger           = cms.bool(True),
         saveAllTriggers     = cms.bool(False), #should be True by default !!
         triggerList         = cms.vstring("HLT_Mu17_Mu8_v17", "HLT_Mu17_TkMu8_v10", "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v17", "HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7", "HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7"),
# ----------------------  Electrons -------------------------------
         doElectrons         = cms.bool(True),
         electron_cut_pt     = cms.double(10),
         electron_cut_eta    = cms.double(2.5),
         electron_saveAllID  = cms.bool(False),
         electron_IDlist     = cms.vstring("mvaTrigV0","mvaNonTrigV0", "eidLoose","eidRobustLoose", "eidTight"),
         electronHLTmatching = cms.vstring(""),
	 electron_Isolist = cms.vstring("PAT","PF03"),
         electron_rhoCorrSrc = cms.vstring("kt6PFJetsForIsolation2011","GammaAndNeutralHadronIso03","Data2012"),
##         electronProducer  = cms.VInputTag(cms.InputTag("selectedPatElectronsPF2PAT")),
         electronProducer    = cms.VInputTag(cms.InputTag("selectedPatElectronsPFlow")),        
	 doElectronMatch     = cms.InputTag(""),
# ----------------------   Photons -------------------------------
         doPhotons           = cms.bool(False),
         photon_cut_pt       = cms.double(10),
         photon_cut_eta      = cms.double(2.5),
         photonHLTmatching   = cms.vstring(""),
##         photonProducer      = cms.VInputTag(cms.InputTag("selectedPatPhotonsPF2PAT")),
         photonProducer      = cms.VInputTag(cms.InputTag("selectedPatPhotonsPFlow")),
# -----------------------   Muons -------------------------------
         doMuons             = cms.bool(True),
         muon_cut_pt         = cms.double(10),
         muon_cut_eta        = cms.double(2.5),
         muon_cut_keepStandaloneMu  = cms.bool(False),
         muon_cut_keepTrackerMu     = cms.bool(True),
         muon_cut_keepCaloMu        = cms.bool(False),
         muon_cut_keepGlobalMu      = cms.bool(True),
         muon_IDlist                = cms.vstring("GlobalMuonPromptTight"),
         muonHLTmatching            = cms.vstring(""),
##         muonProducer        = cms.VInputTag(cms.InputTag("selectedPatMuonsPF2PAT")),
         muonProducer               = cms.VInputTag(cms.InputTag("selectedPatMuonsPFlow")),        
	 muon_Isolist               = cms.vstring("PF03"),
         doMuonMatch                = cms.InputTag(""),
# -----------------------   Taus ----------------------------------
         doTaus              = cms.bool(False),
         tau_cut_pt          = cms.double(10),
         tau_cut_eta         = cms.double(2.4),
         tau_saveAllID       = cms.bool(True),
         tau_IDlist          = cms.vstring(""),
         tauHLTmatching      = cms.vstring(""),
##         tauProducer         = cms.VInputTag(cms.InputTag("selectedPatTausPF2PAT")),
         tauProducer         = cms.VInputTag(cms.InputTag("selectedPatTausPFlow")),
# -----------------------   Tracks ---------------------------------
         doTracks            = cms.bool(False),
         track_cut_pt        = cms.double(0.5),
         track_cut_eta       = cms.double(2.4),
         trackProducer       = cms.VInputTag(cms.InputTag("generalTracks")),
# ----------------------- PFCandidates ---------------------------------
        doPFCandidates = cms.bool(False), # Enabled
        pfcandidate_cut_dR = cms.double(0.3),
        pfcandidate_cut_dz = cms.double(0.05),
        pfcandidate_cut_minPt = cms.double(4.0),
        pfcandidate_VertexTag = cms.VInputTag(cms.InputTag("offlinePrimaryVertices")),
        pfcandidate_InputTag = cms.VInputTag(cms.InputTag("particleFlow")),
# -----------------------  Vertices ---------------------------------
         doVertices          = cms.bool(True),
         saveAllVertex       = cms.bool(True),
         vertexProducer      = cms.VInputTag(cms.InputTag("offlinePrimaryVertices")),
# -----------------------  BeamSpot ---------------------------------
         doBeamSpot          = cms.bool(True),
         beamSpotProducer    = cms.InputTag("offlineBeamSpot"),
# -----------------------   JetMet ----------------------------------
         doJetMet            = cms.bool(True),
         doMuonCorrection    = cms.bool(True),
         jet_cut_pt          = cms.double(10),
         jet_cut_eta         = cms.double(2.5),
         jetIDList           = cms.vstring("LOOSE","TIGHT"),
         jetBTagList         = cms.vstring("trackCountingHighEffBJetTags","jetProbabilityBJetTags","jetBProbabilityBJetTags","combinedSecondaryVertexBJetTags"),
##         jetHLTmatching      = cms.vstring("jetMatchHLTJets"),
         jetHLTmatching      = cms.vstring(""),
         jetmetProducer      = cms.VPSet(
                    cms.PSet(
                 jet    = cms.untracked.string("selectedPatJetsPFlow"),
                 met    = cms.untracked.string("patMETsPFlow"),
                 ##met    = cms.untracked.string("pfType1CorrectedMet"),
                 algo   = cms.untracked.string("pf"),                              
		 fillJetConstituents = cms.untracked.bool(False),
		 fillSubJetConstituents = cms.untracked.bool(False)
               )
         ),
# -----------------------  Pile-Up ----------------------------------
         doPileUp            = cms.bool(True),
         rho_PUUE_dens       = cms.InputTag("kt6PFJets", "rho"),
##         neutralRho_PUUE_dens= cms.InputTag("kt6NeutralPFJets", "rho"),
         neutralRho_PUUE_dens= cms.InputTag("kt6PFJetsCentralNeutral", "rho"),
# -----------------------  MonteCarlo -------------------------------
         doGenParticleCollection = cms.bool(True),
         mcDescentMax = cms.uint32(4),
         mcNGenPartMax = cms.uint32(100),       
         mcTauDecayMode = cms.uint32(2),
         mcHeavyQuarkDecayMode = cms.uint32(0)
)

if options.useData:
    process.MiniTreeProduction.isData                  = cms.bool(True)
    process.MiniTreeProduction.doPileUp                = cms.bool(False)
    process.MiniTreeProduction.doGenParticleCollection = cms.bool(False)
else:
    process.MiniTreeProduction.isData                  = cms.bool(False)
    process.MiniTreeProduction.doGenParticleCollection = cms.bool(True)
    process.MiniTreeProduction.doPileUp                = cms.bool(True)

# loads your analyzer
process.MyModule = cms.EDAnalyzer('NTupleProducer',
         verbose             = cms.uint32(0),   #0: nothing - >1 
# -------------------------------------------------------------------
#                         GENERAL SKIM
# -------------------------------------------------------------------
    general_skim = cms.PSet(
         verbose             = cms.uint32(0),   #0: nothing - >1 
# ----------------------   Trigger -------------------------------
         skimTrigger      = cms.bool(False),
         skimGenParticles = cms.bool(True),
         skimGenTaus      = cms.bool(True),
         triggerList      = cms.vstring(""),
# ----------------------  Electrons -------------------------------
         skimElectrons = cms.bool(True),
         electron_keepAllCollections = cms.bool(True),
         electron_collectionList     = cms.vstring(""),
         electron_pt                 = cms.double(15),
         electron_eta                = cms.double(2.5),
# ----------------------  Muons   -------------------------------
         skimMuons = cms.bool(True),
         muon_keepAllCollections = cms.bool(True),
         muon_collectionList     = cms.vstring(""),
         muon_pt                 = cms.double(15),
         muon_eta                = cms.double(2.5),
# ----------------------  Photons   -------------------------------
         skimPhotons = cms.bool(False),
         photon_keepAllCollections = cms.bool(True),
         photon_collectionList     = cms.vstring(""),
         photon_pt                 = cms.double(7),
         photon_eta                = cms.double(2.5),
# ----------------------  Taus   -------------------------------
         skimTaus = cms.bool(False),
         tau_keepAllCollections = cms.bool(True),
         tau_collectionList     = cms.vstring(""),
         tau_pt                 = cms.double(7),
         tau_eta                = cms.double(2.5),
# ----------------------  Jets   -------------------------------
         skimJets = cms.bool(True),
         jet_keepAllCollections = cms.bool(True),
         jet_collectionList     = cms.vstring("pf"),
         jet_pt                 = cms.double(20),
         jet_eta                = cms.double(2.5),
# ----------------------  Tracks   -------------------------------
         skimTracks = cms.bool(False),
         track_keepAllCollections = cms.bool(True),
         track_collectionList     = cms.vstring(""),
         track_pt                 = cms.double(7),
         track_eta                = cms.double(2.5),
# ----------------------  Vertices   -------------------------------
         skimVertices = cms.bool(False),
         vertex_keepAllCollections = cms.bool(True),
         vertex_collectionList     = cms.vstring("")
     ),
# -------------------------------------------------------------------
#                         TOPDILEPTON SKIM
# -------------------------------------------------------------------
    topdilepton_skim = cms.PSet(
# ----------------------   Trigger -------------------------------
     doTriggerSkimming     = cms.bool(False), # skim on trigger decisions
     triggerSkimList       = cms.vstring("HLT_QuadJet15U"),
# ----------------------  Muons   -------------------------------
     numberOfLept     = cms.int32(1),# for skims ! Total number of 
     numberOfMuon     = cms.int32(-1),# number of sel muon
     muon_cut_pt      = cms.double(10),
     muon_cut_iso     = cms.double(-1),  # PLEASE NO ISO FOR SKIMMING!!!
     useMuonIdSkim     = cms.bool(False),
# ----------------------  Electrons -------------------------------
     numberOfElec       = cms.int32(-1),# number of sel electron
     useElectronIdSkim  = cms.bool(False),
     electron_cut_pt    = cms.double(7),
     electron_cut_iso   = cms.double(-1), # PLEASE NO ISO FOR SKIMMING!!!
# ----------------------  MonteCarlo -------------------------------
     doTMEMESkimming       = cms.bool(False), # skim on the TMEME
     TMEMESkimList         = cms.vint32(),
     doMCDiLepSkimming     = cms.bool(False),
     MCDiLepList           = cms.vstring(""),
# ----------------------  Taus   -------------------------------
     doTauSkimming    = cms.bool(False), # skim on the number of reco 
     numberOfTau      = cms.int32(1),
     tau_cut_pt       = cms.double(5),
     tau_algo         = cms.string("selectedPatTaus"),
# ----------------------  Jets   -------------------------------
     doJetSkimming         = cms.bool(False), # skim on the number of jets
     numberOfJet      = cms.int32(1),
     jet_cut_pt       = cms.double(20),
     jet_cut_eta      = cms.double(2.5),
     jet_algo         = cms.string("pf")
    )
)


process.out.outputCommands = cms.untracked.vstring('drop *',
                            'keep edmTriggerResults_*_*_*',
                            'keep IPHCTreeMTEvent_*_*_*'
                            )

process.TFileService = cms.Service( "TFileService",
                           fileName = cms.string( 'NTuple.root' )
#                           fileName = cms.string( 'MET_Run2012A-13Jul2012-v1_AOD_CMSSW_5_3_3_Ov13.0.root' )
                           )


## The iso-based HBHE noise filter ___________________________________________||
process.load('CommonTools.RecoAlgos.HBHENoiseFilter_cfi')

## The CSC beam halo tight filter ____________________________________________||
process.load('RecoMET.METAnalyzers.CSCHaloFilter_cfi')

## The HCAL laser filter _____________________________________________________||
process.load('RecoMET.METFilters.hcalLaserEventFilter_cfi')

## The ECAL dead cell trigger primitive filter _______________________________||
process.load('RecoMET.METFilters.EcalDeadCellTriggerPrimitiveFilter_cfi') 

## The EE bad SuperCrystal filter ____________________________________________||
process.load('RecoMET.METFilters.eeBadScFilter_cfi') 

## The ECAL laser correction filter
process.load('RecoMET.METFilters.ecalLaserCorrFilter_cfi')

## The tracking failure filter _______________________________________________||
process.load('RecoMET.METFilters.trackingFailureFilter_cfi')

## The tracking POG filters __________________________________________________||
process.load('RecoMET.METFilters.trackingPOGFilters_cff')


## The beam scraping filter __________________________________________________||
process.noscraping = cms.EDFilter(
    "FilterOutScraping",
    applyfilter = cms.untracked.bool(True),
    debugOn = cms.untracked.bool(False),
    numtrack = cms.untracked.uint32(10),
    thresh = cms.untracked.double(0.25)
    )

## The iso-based HBHE noise filter ___________________________________________||
process.load('CommonTools.RecoAlgos.HBHENoiseFilter_cfi')

## The Good vertices collection needed by the tracking failure filter ________||
process.goodVertices = cms.EDFilter(
  "VertexSelector",
  filter = cms.bool(False),
  src = cms.InputTag("offlinePrimaryVertices"),
  cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.rho < 2")
)



###############################
####### DAF PV's     ##########
###############################

pvSrc = 'offlinePrimaryVertices'

## The good primary vertex filter ____________________________________________||
process.primaryVertexFilter = cms.EDFilter(
    "VertexSelector",
    src = cms.InputTag("offlinePrimaryVertices"),
    cut = cms.string("!isFake & ndof > 4 & abs(z) <= 24 & position.Rho < 2"),
    filter = cms.bool(True)
    )


from PhysicsTools.SelectorUtils.pvSelector_cfi import pvSelector

process.goodOfflinePrimaryVertices = cms.EDFilter(
    "PrimaryVertexObjectFilter",
    filterParams = pvSelector.clone( maxZ = cms.double(24.0),
                                     minNdof = cms.double(4.0) # this is >= 4
                                     ),
    src=cms.InputTag(pvSrc)
    )
    
process.trackingFailureFilter.VertexSource = 'primaryVertexFilter'
  
process.filtersSeq = cms.Sequence(
   process.primaryVertexFilter *
   #process.noscraping *
   process.HBHENoiseFilter* 
   process.CSCTightHaloFilter *
   process.hcalLaserEventFilter *
   process.EcalDeadCellTriggerPrimitiveFilter *
   process.trackingFailureFilter *
   process.eeBadScFilter *
   process.ecalLaserCorrFilter *
   process.trkPOGFilters
)

process.patseq = cms.Sequence(
    process.filtersSeq* 
    process.type0PFMEtCorrection*
    process.producePFMETCorrections*
    process.goodOfflinePrimaryVertices*
    getattr(process,"patPF2PATSequence"+postfix)
    )


if options.runOnFastSim:
    process.patseq.remove( process.HBHENoiseFilter )


process.p = cms.Path(
       process.patTrigger                      
       * process.patTriggerEvent 
       * process.patseq
       * process.MiniTreeProduction
       * process.MyModule
       #* process.rootNTuples
)


process.out.outputCommands = cms.untracked.vstring('drop *')#,
                                                   #'keep IPHCTreeMTEvent_*_*_*',
                                                   #'keep *_MiniTreeSkimming_*_*')


process.MessageLogger.cerr.FwkReport.reportEvery = 100

