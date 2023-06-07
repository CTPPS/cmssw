import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Run3_cff import Run3

# needed to have ctppsDQMCalibrationSource properly working
#from Configuration.Eras.Modifier_ctpps_cff import ctpps

process = cms.Process('RECO', Run3)

# minimum of logs
process.MessageLogger = cms.Service("MessageLogger",
  statistics = cms.untracked.vstring(),
  destinations = cms.untracked.vstring('cerr'),
  cerr = cms.untracked.PSet(
      threshold = cms.untracked.string('WARNING')
  )
)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

# load DQM framework
#process.load("DQM.Integration.config.environment_cfi")
#process.dqmEnv.subSystemFolder = "CTPPS"
#process.dqmEnv.eventInfoFolder = "EventInfo"
#process.dqmSaver.path = ""
#process.dqmSaver.tag = "CTPPS"

# raw data source, T2 + TrackingStrips either in run 368593 (correct) or 368594 (no)
process.source = cms.Source("NewEventStreamFileReader",
    fileNames = cms.untracked.vstring(
        '/store/t0streamer/Minidaq/A/000/368/593/run368593_ls0001_streamA_StorageManager.dat',
#        '/store/t0streamer/Minidaq/A/000/368/594/run368594_ls0001_streamA_StorageManager.dat',
    )
)
# raw data source
#alignment run without T2, TimingDiamond and TrackingStrip not affected by PR changes


process.maxEvents = cms.untracked.PSet(
  input = cms.untracked.int32(1000)
)

# raw-to-digi conversion
process.load("EventFilter.CTPPSRawToDigi.ctppsRawToDigi_cff")

# global tag - conditions for P5 cluster
#process.load("DQM.Integration.config.FrontierCondition_GT_cfi")

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '130X_dataRun3_Prompt_Candidate_2023_06_06_21_34_08', '')

# local RP reconstruction chain with standard settings
process.load("RecoPPS.Configuration.recoCTPPS_cff")

# CTPPS DQM modules
#process.load("DQM.CTPPS.ctppsDQM_cff")

process.path = cms.Path(
  process.ctppsRawToDigi *
  process.recoCTPPS 
)

process.output = cms.OutputModule("PoolOutputModule",
        fileName = cms.untracked.string("file:output-commonReco-368593-newPCGT-1000ev.root"),
    outputCommands = cms.untracked.vstring(
        'keep *',
    ),
)

process.outpath = cms.EndPath(process.output)

process.schedule = cms.Schedule(
  process.path,
  process.outpath
)
