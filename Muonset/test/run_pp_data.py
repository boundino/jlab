import FWCore.ParameterSet.Config as cms
process = cms.Process('HiForest')
process.options = cms.untracked.PSet()

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
                            'file:/afs/cern.ch/work/w/wangj/public/DoubleMuon/Run2015E-PromptReco-v1/6E4E5876-9291-E511-9E3C-02163E01417C.root'
					)
)


# Number of events we want to process, -1 = all events
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(200))

process.load('Configuration.StandardSequences.Services_cff')
# process.load('Configuration.Geometry.GeometryRecoDB_cff')
# process.load('Configuration.StandardSequences.MagneticField_38T_cff')
# process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
#process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load('Configuration.StandardSequences.ReconstructionHeavyIons_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")

from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data', '')
# process.GlobalTag = GlobalTag(process.GlobalTag, '75X_dataRun2_v13', '')  #for now track GT manually, since centrality tables updated ex post facto
# process.HiForest.GlobalTagLabel = process.GlobalTag.globaltag

process.load("RecoHI.HiCentralityAlgos.CentralityBin_cfi")
process.centralityBin.Centrality = cms.InputTag("hiCentrality")
process.centralityBin.centralityVariable = cms.string("HFtowers")

process.load('HeavyIonsAnalysis.EventAnalysis.hltanalysis_cff')
from HeavyIonsAnalysis.EventAnalysis.dummybranches_cff import addHLTdummybranchesForPP
addHLTdummybranchesForPP(process)
process.ana_step = cms.Path(process.hltanalysis)


process.TFileService = cms.Service("TFileService",
                                   fileName=cms.string("HiForestAOD_DATA.root"))

runOnMC = False
VtxLabel = "offlinePrimaryVerticesWithBS"
# TrkLabel = "generalTracks"
from Muonset.muonfinder.muonfinder_cff import finderMaker_75X
finderMaker_75X(process, runOnMC, VtxLabel)

process.p = cms.Path(process.MuonsetSequence)
