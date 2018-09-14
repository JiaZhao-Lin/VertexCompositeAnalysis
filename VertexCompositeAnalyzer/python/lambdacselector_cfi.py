import FWCore.ParameterSet.Config as cms

lambdacselector = cms.EDProducer('VertexCompositeSelector',
  doGenMatching = cms.untracked.bool(False),
  hasSwap = cms.untracked.bool(False),
  decayInGen = cms.untracked.bool(False),
  twoLayerDecay = cms.untracked.bool(True),
  #PID used only for GEN and/or GEN match
  PID = cms.untracked.int32(4122),
  PID_dau1 = cms.untracked.int32(310),
  PID_dau2 = cms.untracked.int32(2212),
  VertexCollection = cms.untracked.InputTag("offlinePrimaryVertices"),
  TrackCollection = cms.untracked.InputTag("generalTracks"),
  VertexCompositeCollection = cms.untracked.InputTag("generalLambdaCCandidatesNew:LambdaCToKsP"),
  GenParticleCollection = cms.untracked.InputTag("genParticles"),
  MuonCollection = cms.untracked.InputTag("null"),
  doMuon = cms.untracked.bool(False),

  selectGenMatch = cms.untracked.bool(False),
  selectGenUnMatch = cms.untracked.bool(False),
  selectGenMatchSwap = cms.untracked.bool(False),
  selectGenMatchUnSwap = cms.untracked.bool(False),

  useAnyMVA = cms.bool(False),
  useExistingMVA = cms.bool(False),
  mvaType = cms.string('BDT'),
  GBRForestLabel = cms.string(''),
  GBRForestFileName = cms.string(''),
  MVACollection = cms.InputTag("generalLambdaCCandidatesNew:MVAValues"),
  mvaMax = cms.untracked.double(999.9),
  mvaMin = cms.untracked.double(-999.9),

  trkPMin = cms.untracked.double(0.),
  trkPtMin = cms.untracked.double(0.),
  trkEtaMax = cms.untracked.double(999.),
  trkPSumMin = cms.untracked.double(0.),
  trkPtSumMin = cms.untracked.double(0.),
  trkPtAsymMin = cms.untracked.double(0.),
  trkEtaDiffMax = cms.untracked.double(999.),
  trkPtErrMax = cms.untracked.double(999.),
  trkNHitMin = cms.untracked.int32(0),
  candpTMin = cms.untracked.double(-999.),
  candpTMax = cms.untracked.double(999.),
  candYMin = cms.untracked.double(-999.),
  candYMax = cms.untracked.double(999.),
  cand3DDecayLengthSigMin = cms.untracked.double(0.),
  cand3DPointingAngleMax = cms.untracked.double(999.),
  candVtxProbMin = cms.untracked.double(0.)
                              )

lambdacselectorMC = cms.EDProducer('VertexCompositeSelector',

  doGenMatching = cms.untracked.bool(True),
  hasSwap = cms.untracked.bool(False),
  decayInGen = cms.untracked.bool(True),
  twoLayerDecay = cms.untracked.bool(True),
  #PID used only for GEN and/or GEN match
  PID = cms.untracked.int32(4122),
  PID_dau1 = cms.untracked.int32(310),
  PID_dau2 = cms.untracked.int32(2212),
  deltaR = cms.untracked.double(0.03),

  VertexCollection = cms.untracked.InputTag("offlinePrimaryVertices"),
  TrackCollection = cms.untracked.InputTag("generalTracks"),
  VertexCompositeCollection = cms.untracked.InputTag("generalLambdaCCandidatesNew:LambdaCToKsP"),
  GenParticleCollection = cms.untracked.InputTag("genParticles"),
  MuonCollection = cms.untracked.InputTag("null"),
  doMuon = cms.untracked.bool(False),

  selectGenMatch = cms.untracked.bool(False),
  selectGenUnMatch = cms.untracked.bool(False),
  selectGenMatchSwap = cms.untracked.bool(False),
  selectGenMatchUnSwap = cms.untracked.bool(False),

  useAnyMVA = cms.bool(False),
  useExistingMVA = cms.bool(False),
  mvaType = cms.string('BDT'),
  GBRForestLabel = cms.string(''),
  GBRForestFileName = cms.string(''),
  MVACollection = cms.InputTag("generalLambdaCCandidatesNew:MVAValues"),
  mvaMax = cms.untracked.double(999.9),
  mvaMin = cms.untracked.double(-999.9),

  trkPMin = cms.untracked.double(0.),
  trkPtMin = cms.untracked.double(0.),
  trkEtaMax = cms.untracked.double(999.),
  trkPSumMin = cms.untracked.double(0.),
  trkPtSumMin = cms.untracked.double(0.),
  trkPtAsymMin = cms.untracked.double(0.),
  trkEtaDiffMax = cms.untracked.double(999.),
  trkPtErrMax = cms.untracked.double(999.),
  trkNHitMin = cms.untracked.int32(0),
  candpTMin = cms.untracked.double(-999.),
  candpTMax = cms.untracked.double(999.),
  candYMin = cms.untracked.double(-999.),
  candYMax = cms.untracked.double(999.),
  cand3DDecayLengthSigMin = cms.untracked.double(0.),
  cand3DPointingAngleMax = cms.untracked.double(999.),
  candVtxProbMin = cms.untracked.double(0.)
                              )
