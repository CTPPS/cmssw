import FWCore.ParameterSet.Config as cms
from Configuration.Eras.Era_Run3_cff import Run3

process = cms.Process('RECODQM', Run3)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(90000) )
process.verbosity = cms.untracked.PSet( input = cms.untracked.int32(10) ) 

# minimum of logs
process.MessageLogger = cms.Service("MessageLogger",
    cerr = cms.untracked.PSet(
        threshold = cms.untracked.string('INFO')
    ) 
)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

# load DQM framework
process.load("DQM.Integration.config.environment_cfi")
process.dqmEnv.subSystemFolder = "CTPPS"
process.dqmEnv.eventInfoFolder = "EventInfo"
process.dqmSaver.path = ""
process.dqmSaver.tag = "CTPPS"
#process.dqmSaver.runNumber = 999999


# raw data source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#'file:/eos/cms/store/data/Run2022B/AlCaPPS/RAW/v1/000/355/207/00000/c23440f4-49c0-44aa-b8f6-f40598fb4705.root',
#'file:/eos/cms/store/group/dpg_ctpps/comm_ctpps/ALCAPPS-RAW/000/355/988/00000/7728d214-b404-4578-b781-b0c207cfb875.root',
#'/store/data/Run2022D/AlCaPPSPrompt/ALCARECO/PPSCalMaxTracks-PromptReco-v2/000/357/900/00000/0434a2eb-2cea-4d83-8b21-fe91b755e62e.root',
#  RUN 440 - "357440": [[1,354]],
"/store/data/Run2022C/AlCaPPSPrompt/ALCARECO/PPSCalMaxTracks-PromptReco-v1/000/357/440/00000/53746080-d21b-4fdf-9d19-67c6dae347ac.root",
"/store/data/Run2022C/AlCaPPSPrompt/ALCARECO/PPSCalMaxTracks-PromptReco-v1/000/357/440/00000/b8d3f012-0b59-43b4-adbc-811bcb59c9c4.root",
"/store/data/Run2022C/AlCaPPSPrompt/ALCARECO/PPSCalMaxTracks-PromptReco-v1/000/357/440/00000/cfa9f82a-5296-49cf-8353-11f8761b675c.root",
"/store/data/Run2022C/AlCaPPSPrompt/ALCARECO/PPSCalMaxTracks-PromptReco-v1/000/357/440/00000/db3c6c33-9ad0-4d61-9a64-a1aa9578550e.root"    ),
)



from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_hlt_relval', '')
process.GlobalTag = GlobalTag(process.GlobalTag, '124X_dataRun3_Prompt_frozen_v4', '')
process.GlobalTag.toGet = cms.VPSet()
process.GlobalTag.toGet.append(
  cms.PSet(record = cms.string("PPSTimingCalibrationRcd"),
           tag =  cms.string("PPSDiamondTimingCalibration_Run3_recovered_v1"),
           label = cms.untracked.string('PPSTestCalibration'),
           connect = cms.string("frontier://FrontierPrep/CMS_CONDITIONS")
	)
)


# raw-to-digi conversion
process.load("EventFilter.CTPPSRawToDigi.ctppsRawToDigi_cff")

process.ctppsPixelDigis.inputLabel = cms.InputTag("hltPPSCalibrationRaw")
process.ctppsDiamondRawToDigi.rawDataTag = cms.InputTag("hltPPSCalibrationRaw")
process.totemRPRawToDigi.rawDataTag = cms.InputTag("hltPPSCalibrationRaw")
process.totemTimingRawToDigi.rawDataTag = cms.InputTag("hltPPSCalibrationRaw")


# local RP reconstruction chain with standard settings
process.load("RecoPPS.Configuration.recoCTPPS_cff")
process.ctppsDiamondRecHits.digiTag=cms.InputTag('ctppsDiamondRawToDigiAlCaRecoProducer:TimingDiamond')
process.ctppsDiamondRecHits.timingCalibrationTag=cms.string("GlobalTag:PPSTestCalibration")

#process.load('CondCore.CondDB.CondDB_cfi')
#process.CondDB.connect = 'sqlite_file:ppsDiamondTiming_calibration.sqlite' # SQLite input
#process.PoolDBESSource = cms.ESSource('PoolDBESSource',
#        process.CondDB,
#        DumpStats = cms.untracked.bool(True),
#        toGet = cms.VPSet(
#            cms.PSet(
#                record = cms.string('PPSTimingCalibrationRcd'),
#                tag = cms.string('DiamondTimingCalibration') 
#            )
#        )
#    )

#process.load('Geometry.VeryForwardGeometry.geometryRPFromDD_2021_cfi')
# CTPPS DQM modules
process.load("DQM.CTPPS.ctppsDQM_cff")
process.ctppsDiamondDQMSource.tagDigi=cms.InputTag('ctppsDiamondRawToDigiAlCaRecoProducer:TimingDiamond')
process.ctppsDiamondDQMSource.tagFEDInfo=cms.InputTag('ctppsDiamondRawToDigiAlCaRecoProducer:TimingDiamond')
process.ctppsDiamondDQMSource.tagStatus=cms.InputTag('ctppsDiamondRawToDigiAlCaRecoProducer:TimingDiamond')
process.ctppsDiamondDQMSource.tagPixelLocalTracks=cms.InputTag('ctppsPixelLocalTracksAlCaRecoProducer')
#process.ctppsDiamondDQMSource.excludeMultipleHits = cms.bool(True)
process.ctppsDiamondDQMSource.plotOnline = cms.untracked.bool(True)
process.ctppsDiamondDQMSource.plotOffline = cms.untracked.bool(False)

process.content = cms.EDAnalyzer("EventContentAnalyzer")

from DQMServices.Core.DQMEDAnalyzer import DQMEDAnalyzer
process.tdcCalibShiftDQM = DQMEDAnalyzer('TdcCalibShiftDQM')
process.path = cms.Path(
    #process.ctppsRawToDigi*
    #process.recoCTPPS*
    #process.ctppsDQMOnlineSource*
    ##process.ctppsDiamondDQMSource*
    #process.ctppsDQMOnlineHarvest
    #process.tdcCalibShiftDQM
    
    #process.content*
    process.ctppsDiamondRecHits *
    process.ctppsDiamondLocalTracks *
    #process.ctppsLocalTrackLiteProducer *
    #process.ctppsProtons *
    process.ctppsDiamondDQMSource*
    process.ctppsDQMOnlineHarvest
)

process.output = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string("file:test.root"),
    outputCommands = cms.untracked.vstring(
        'drop *',
        'keep *_ctpps*_*_*',
    ),
)


process.end_path = cms.EndPath(
    #process.output
    process.dqmEnv +
    process.dqmSaver
)

process.schedule = cms.Schedule(
    process.path,
    process.end_path
) 

#process.output = cms.OutputModule("PoolOutputModule",
#    fileName = cms.untracked.string("AOD.root"),
#    outputCommands = cms.untracked.vstring(
#        'drop *',
#        'keep *_ctpps*_*_*',
#    ),
#)  
#process.outpath = cms.EndPath(process.output)