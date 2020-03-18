import FWCore.ParameterSet.Config as cms

def finderMaker_75X(process, runOnMC = True, VtxLabel = "hiSelectedVertex"):
	process.load("PhysicsTools.PatAlgos.patSequences_cff")
	process.particleFlowPtrs.src = "particleFlowTmp"
	process.pfPileUpIsoPFBRECO.Vertices = cms.InputTag(VtxLabel)
	process.pfPileUpPFBRECO.Vertices = cms.InputTag(VtxLabel)
	## patMuonsWithTrigger
	process.load("MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff")
	from MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff import addMCinfo, useL1MatchingWindowForSinglets, changeTriggerProcessName, switchOffAmbiguityResolution, addHLTL1Passthrough

	process.patMuonsWithTriggerSequence = cms.Sequence(process.patMuonsWithTriggerSequence)
	process.patMuonsWithoutTrigger.isoDeposits = cms.PSet()
	process.patMuonsWithoutTrigger.pvSrc = cms.InputTag(VtxLabel)
	process.patMuonsWithoutTrigger.isolationValues = cms.PSet()
	changeTriggerProcessName(process, "HLT")
	switchOffAmbiguityResolution(process) # Switch off ambiguity resolution: allow multiple reco muons to match to the same trigger muon
	addHLTL1Passthrough(process)

	process.patTrigger.collections.remove("hltL3MuonCandidates")
	process.patTrigger.collections.append("hltHIL3MuonCandidates")

	process.muonL1Info.maxDeltaR = 0.3
	process.muonL1Info.fallbackToME1 = True
	process.muonMatchHLTL1.maxDeltaR = 0.3
	process.muonMatchHLTL1.fallbackToME1 = True
	process.muonMatchHLTL2.maxDeltaR = 0.3
	process.muonMatchHLTL2.maxDPtRel = 10.0
	process.muonMatchHLTL3.maxDeltaR = 0.1
	process.muonMatchHLTL3.maxDPtRel = 10.0
	process.muonMatchHLTCtfTrack.maxDeltaR = 0.1
	process.muonMatchHLTCtfTrack.maxDPtRel = 10.0
	process.muonMatchHLTTrackMu.maxDeltaR = 0.1
	process.muonMatchHLTTrackMu.maxDPtRel = 10.0
	process.muonMatchHLTL3.matchedCuts = cms.string('coll("hltHIL3MuonCandidates")')

	process.Muonset = cms.EDAnalyzer('Muonset',
        MuonLabel = cms.InputTag('patMuonsWithTrigger'),
        BSLabel = cms.InputTag("offlineBeamSpot"),
        PVLabel = cms.InputTag(VtxLabel),
	)

	process.MuonsetSequence = cms.Sequence(process.patMuonsWithTriggerSequence*process.Muonset)
