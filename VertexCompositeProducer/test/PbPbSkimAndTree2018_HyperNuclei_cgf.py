import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras
process = cms.Process('ANASKIM',eras.Phase2C4_timing_layer_bar)

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.Geometry.GeometryExtended2023D35Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')

# Limit the output messages
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 200
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

# Define the input source
process.source = cms.Source("PoolSource",
   fileNames = cms.untracked.vstring('/store/user/anstahll/MTD/MC/20201103/PY8EGun_TuneCP5_4LambdaH_p12_eta3p2_To4HePi_Pythia8_pp_14TeV_RECO_20201103/MTD/PY8EGun_TuneCP5_4LambdaH_p12_eta3p2_To4HePi_Pythia8_pp_14TeV_RECO_20201103/201103_080620/0000/step3_PY8EGun_TuneCP5_4LambdaH_p12_eta3p2_To4HePi_Pythia8_pp_14TeV_RECO_20201103_1.root'),
)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1))

# Set the global tag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.GlobalTag.globaltag = cms.string('103X_upgrade2023_realistic_v2')

# Add PbPb centrality
process.GlobalTag.snapshotTime = cms.string("9999-12-31 23:59:59.000")
process.GlobalTag.toGet.extend([
    cms.PSet(record = cms.string("HeavyIonRcd"),
        tag = cms.string("CentralityTable_HFtowers200_HydjetTuneCP5MTD_v1040mtd4x1_mc"),
        connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"),
        label = cms.untracked.string("HFtowers")
        ),
    ])
process.load('RecoHI.HiCentralityAlgos.HiCentrality_cfi')
process.hiCentrality.produceHFhits = False
process.hiCentrality.produceHFtowers = True
process.hiCentrality.produceEcalhits = False
process.hiCentrality.produceZDChits = False
process.hiCentrality.produceETmidRapidity = False
process.hiCentrality.producePixelhits = False
process.hiCentrality.produceTracks = False
process.hiCentrality.producePixelTracks = False
process.hiCentrality.reUseCentrality = False
process.hiCentrality.srcReUse = cms.InputTag("hiCentrality","","RECO")
process.hiCentrality.srcTracks = cms.InputTag("generalTracks")
process.hiCentrality.srcVertex = cms.InputTag("offlinePrimaryVertices4D")
process.load("RecoHI.HiCentralityAlgos.CentralityBin_cfi")
process.centralityBin.Centrality = cms.InputTag("hiCentrality")
process.centralityBin.centralityVariable = cms.string("HFtowers")
process.centralityBin.nonDefaultGlauberModel = cms.string("")
process.hiCentrality.srcEBhits = cms.InputTag("HGCalRecHit","HGCHEBRecHits")
process.hiCentrality.srcEEhits = cms.InputTag("HGCalRecHit","HGCEERecHits")
process.cent_seq = cms.Sequence(process.hiCentrality * process.centralityBin)

# Add the VertexComposite producer
from VertexCompositeAnalysis.VertexCompositeProducer.generalParticles_cff import generalParticles
process.generalCandidates = generalParticles.clone(
    pdgId = cms.int32(1010010040),
    mass = cms.double(3.9226),
    doSwap = cms.bool(True),
    width = cms.double(0.3),

    dedxHarmonic2 = cms.InputTag('dedxPixelHarmonic2'),
    mtdValueNames = cms.vstring(["tMTD", "tMTDErr", "pathLength"]),
    mtdValueLabels = cms.vstring(["trackExtenderWithMTD:generalTracktmtd:RECO", "trackExtenderWithMTD:generalTracksigmatmtd:RECO", "trackExtenderWithMTD:generalTrackPathLength:RECO"]),

    preSelection = cms.string(""
       "charge==0"
       ),
    pocaSelection = cms.string(""
       "userFloat('bestMass') >= 3.51 && userFloat('bestMass') <= 4.32 && pt >= 0.0"
       "&& userFloat('dca') >= 0 && userFloat('dca') <= 9999."
       ),
    postSelection = cms.string(""
       "userFloat('vertexProb') >= 0.02"
       "&& userFloat('normChi2') <= 9999.0"
       ),
    finalSelection = cms.string(""
       "abs(rapidity) < 3.0"
       ),

    # daughter information
    daughterInfo = cms.VPSet([
        cms.PSet(pdgId = cms.int32(211), charge = cms.int32(-1),
           selection = cms.string(
              "p>=0.7 && abs(eta)<=3.0"
              "&& quality('highPurity') && ptError/pt<0.1"
              "&& normalizedChi2<7."
              "&& numberOfValidHits >=11"
              ),
           finalSelection = cms.string(''
              'abs(userFloat("tMTDErr")) > 0.0'
              )
        ),
        cms.PSet(pdgId = cms.int32(1000020040), mass = cms.double(3.72742), charge = cms.int32(+1),
           selection = cms.string(
              "p>=0.7 && abs(eta)<=3.0"
              "&& quality('highPurity') && ptError/pt<0.1"
              "&& normalizedChi2<7."
              "&& numberOfValidHits >=11"
              ),
           finalSelection = cms.string(''
              'abs(userFloat("tMTDErr")) > 0.0'
              )
        )
    ]),

)

# Define the analysis steps
process.pcentandep_step = cms.Path(process.cent_seq)
process.rereco_step = cms.Path(process.generalCandidates)

# Add the VertexComposite tree
from VertexCompositeAnalysis.VertexCompositeAnalyzer.particle_tree_cff import particleAna_mc
process.ana = particleAna_mc.clone(
  recoParticles = cms.InputTag("generalCandidates"),
)

# Define the output
process.TFileService = cms.Service("TFileService", fileName = cms.string('ana.root'))
process.p = cms.EndPath(process.ana)

# Define the process schedule
process.schedule = cms.Schedule(
    process.pcentandep_step,
    process.rereco_step,
    process.p
)
