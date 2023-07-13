import FWCore.ParameterSet.Config as cms

# load era modifier to run on 2018 data
from Configuration.Eras.Modifier_ctpps_2022_cff import ctpps_2022

process = cms.Process( 'TEST',ctpps_2022)

# LHCInfo plotter
process.load("Validation.CTPPS.ctppsLHCInfoPlotter_cfi")
process.ctppsLHCInfoPlotter.outputFile = "outALCAPPS_RECO_EXPRESS_Fill9019/alcareco_lhc_info_express.root"

# Load geometry from DB
process.load("Geometry.VeryForwardGeometry.geometryRPFromDB_cfi")

# track distribution plotter
process.ctppsTrackDistributionPlotter = cms.EDAnalyzer("CTPPSTrackDistributionPlotter",
  tagTracks = cms.InputTag("ctppsLocalTrackLiteProducerAlCaRecoProducer"),

  rpId_45_F = cms.uint32(23),
  rpId_45_N = cms.uint32(3),
  rpId_56_N = cms.uint32(103),
  rpId_56_F = cms.uint32(123),

  outputFile = cms.string("outALCAPPS_RECO_EXPRESS_Fill9019/alcareco_tracks_express.root")
)

# reconstruction plotter
process.ctppsProtonReconstructionPlotter = cms.EDAnalyzer("CTPPSProtonReconstructionPlotter",
  tagTracks = cms.InputTag("ctppsLocalTrackLiteProducerAlCaRecoProducer"),
  tagRecoProtonsSingleRP = cms.InputTag("ctppsProtonsAlCaRecoProducer", "singleRP"),
  tagRecoProtonsMultiRP = cms.InputTag("ctppsProtonsAlCaRecoProducer", "multiRP"),

  rpId_45_F = cms.uint32(23),
  rpId_45_N = cms.uint32(3),
  rpId_56_N = cms.uint32(103),
  rpId_56_F = cms.uint32(123),

  outputFile = cms.string("outALCAPPS_RECO_EXPRESS_Fill9019/alcareco_protons_express.root")
)

process.p = cms.Path(
  process.ctppsLHCInfoPlotter
  * process.ctppsTrackDistributionPlotter
  * process.ctppsProtonReconstructionPlotter
)

process.options = cms.PSet(
    wantSummary = cms.untracked.bool( True ),
)

# load GT
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run3_data_express')

process.source = cms.Source( 'PoolSource',
    fileNames = cms.untracked.vstring(
        'file:outALCAPPS_RECO_EXPRESS_Fill9019/outputALCAPPS_RECO_express_new.root',
    ),
    inputCommands = cms.untracked.vstring(
        'keep *'
    )
)

# limit the number of events to be processed
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32( -1 )
)
