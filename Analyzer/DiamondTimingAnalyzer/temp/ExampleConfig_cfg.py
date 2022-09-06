
import FWCore.ParameterSet.Config as cms
import copy

process = cms.Process('PPSTiming2')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(30000)
)

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = ''
process.MessageLogger.cerr.FwkReport.reportEvery = 1000




process.source = cms.Source ("PoolSource",
                             fileNames = cms.untracked.vstring(
    "file:/eos/home-m/mobrzut/357/734/00000/73c5eeec-578d-438f-92f1-3e1fff2ed390.root",   
    "file:/eos/home-m/mobrzut/357/734/00000/08e5ee70-ebff-4571-b14c-806961a3d9dd.root",
    "file:/eos/home-m/mobrzut/357/734/00000/d42eaf2c-57bb-48fe-9a22-317053f89115.root"
# "file:/eos/home-m/mobrzut/357/479/00000/68a5a64a-b756-4cf5-9a2a-00181e34f501.root"
# RUN 442 - "357442": [[1,1392]],
# "/store/data/Run2022C/AlCaPPSPrompt/ALCARECO/PPSCalMaxTracks-PromptReco-v1/000/357/442/00000/19b41c28-cc6e-4df9-b73a-ae8b187021c7.root",
# "/store/data/Run2022C/AlCaPPSPrompt/ALCARECO/PPSCalMaxTracks-PromptReco-v1/000/357/442/00000/4d2b38b9-f03f-4887-ad72-e53e448cc52d.root"
#"/store/data/Run2022C/AlCaPPSPrompt/ALCARECO/PPSCalMaxTracks-PromptReco-v1/000/357/442/00000/5987c37e-850f-4c00-9e69-d6a3ed8df910.root",
#"/store/data/Run2022C/AlCaPPSPrompt/ALCARECO/PPSCalMaxTracks-PromptReco-v1/000/357/442/00000/9450abb5-9478-4209-b1dd-7ff06beab620.root",
#"/store/data/Run2022C/AlCaPPSPrompt/ALCARECO/PPSCalMaxTracks-PromptReco-v1/000/357/442/00000/a2536da7-b15f-4d3c-be30-9feefb303f41.root",
#"/store/data/Run2022C/AlCaPPSPrompt/ALCARECO/PPSCalMaxTracks-PromptReco-v1/000/357/442/00000/a2f71e31-b36b-4a05-bb14-17e6b404ceec.root",
#"/store/data/Run2022C/AlCaPPSPrompt/ALCARECO/PPSCalMaxTracks-PromptReco-v1/000/357/442/00000/ca8d0e9c-5ac7-4987-ba7b-b771cd933485.root",
#"/store/data/Run2022C/AlCaPPSPrompt/ALCARECO/PPSCalMaxTracks-PromptReco-v1/000/357/442/00000/e6a28ca6-f978-430c-bacf-7b0d59f8c069.root",
#"/store/data/Run2022C/AlCaPPSPrompt/ALCARECO/PPSCalMaxTracks-PromptReco-v1/000/357/442/00000/eac6865b-3265-4d9d-aee9-579dc95e7f63.root",
#"/store/data/Run2022C/AlCaPPSPrompt/ALCARECO/PPSCalMaxTracks-PromptReco-v1/000/357/442/00000/ee2f7caa-dcb5-4ab3-b6cf-b302031dd105.root"
)
)

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
# 2022 prompt: to be updated
process.GlobalTag = GlobalTag(process.GlobalTag, "124X_dataRun3_v9") 


#JH - use new tag for timing calibrations



                                        
process.load('CondCore.CondDB.CondDB_cfi')
process.CondDB.connect = 'sqlite_file:/afs/cern.ch/user/m/mobrzut/public/timing_calibration_env/CMSSW_12_4_6/src/Analyzer/DiamondTimingAnalyzer/temp/ppsDiamondTiming_calibration357734.sqlite' # SQLite input TODO: migrate to using tag
process.PoolDBESSource = cms.ESSource('PoolDBESSource',
        process.CondDB,
        DumpStats = cms.untracked.bool(True),
        toGet = cms.VPSet(
            cms.PSet(
                record = cms.string('PPSTimingCalibrationRcd'),
                tag = cms.string('DiamondTimingCalibration'),
                           label = cms.untracked.string('PPSTestCalibration'), 
            )
        )
    )

apply_shift = True                                        

# JH - rerun reco sequence with new timing conditions
process.load("RecoPPS.Configuration.recoCTPPS_cff")
process.ctppsDiamondLocalTracks.recHitsTag = cms.InputTag("ctppsDiamondRecHits","","PPSTiming2")
process.ctppsLocalTrackLiteProducer.tagDiamondTrack = cms.InputTag("ctppsDiamondLocalTracks","","PPSTiming2")
process.ctppsProtons.tagLocalTrackLite = cms.InputTag("ctppsLocalTrackLiteProducer","","PPSTiming2")
process.ctppsLocalTrackLiteProducer.includeDiamonds = cms.bool(True)
process.ctppsLocalTrackLiteProducer.tagPixelTrack = cms.InputTag("ctppsPixelLocalTracksAlCaRecoProducer")
process.ctppsLocalTrackLiteProducer.includePixels = cms.bool(True)
process.ctppsDiamondRecHits.timingCalibrationTag=cms.string("PoolDBESSource:PPSTestCalibration")
process.ctppsDiamondRecHits.digiTag= cms.InputTag("ctppsDiamondRawToDigiAlCaRecoProducer:TimingDiamond")

from DQMServices.Core.DQMEDAnalyzer import DQMEDAnalyzer
if(apply_shift):
    process.diamondTimingWorker = DQMEDAnalyzer("DiamondTimingWorker",
        tagDigi = cms.InputTag("ctppsDiamondRawToDigiAlCaRecoProducer", "TimingDiamond"),
        tagRecHit = cms.InputTag("ctppsDiamondRecHits","","PPSTiming2"),
        tagPixelLocalTrack = cms.InputTag("ctppsPixelLocalTracksAlCaRecoProducer"),
        timingCalibrationTag=cms.string("PoolDBESSource:PPSTestCalibration"),
        tagLocalTrack = cms.InputTag("ctppsDiamondLocalTracks","","PPSTiming2"),
        tagValidOOT = cms.int32(-1),
        planesConfig = cms.string("planes.json"),
        Ntracks_Lcuts = cms.vint32([-1,1,-1,1]), # minimum number of tracks in pots [45-210, 45-220, 56-210, 56-220]
        Ntracks_Ucuts = cms.vint32([-1,6,-1,6]), # maximum number of tracks in pots [45-210, 45-220, 56-210, 56-220]
    )
else:
    process.diamondTimingWorker = DQMEDAnalyzer("DiamondTimingWorker",
    tagDigi = cms.InputTag("ctppsDiamondRawToDigiAlCaRecoProducer", "TimingDiamond"),
    tagRecHit = cms.InputTag("ctppsDiamondRecHitsAlCaRecoProducer"),
    tagPixelLocalTrack = cms.InputTag("ctppsPixelLocalTracksAlCaRecoProducer"),
    tagLocalTrack = cms.InputTag("ctppsDiamondLocalTracksAlCaRecoProducer"),
    timingCalibrationTag=cms.string("PoolDBESSource:PPSTestCalibration"),
    tagValidOOT = cms.int32(-1),
    planesConfig = cms.string("planes.json"),
    Ntracks_Lcuts = cms.vint32([-1,1,-1,1]), # minimum number of tracks in pots [45-210, 45-220, 56-210, 56-220]
    Ntracks_Ucuts = cms.vint32([-1,6,-1,6]), # maximum number of tracks in pots [45-210, 45-220, 56-210, 56-220]
) 

if(apply_shift):
    process.ALL = cms.Path(
        # Re-run the PPS local+proton reconstruction starting from AOD with new timing calibrations
        process.ctppsDiamondRecHits *
        process.ctppsDiamondLocalTracks *
        process.diamondTimingWorker
                        )
else: 
    process.ALL = cms.Path(process.diamondTimingWorker)

saveToDQM = False
if(saveToDQM):
    process.load("DQMServices.Core.DQMStore_cfi")
    process.load("DQMServices.Components.DQMEnvironment_cfi")
    process.dqmEnv.subSystemFolder = "CalibPPS"
    process.dqmSaver.convention = 'Offline'
    process.dqmSaver.workflow = "/CalibPPS/TimingCalibration/CMSSW_12_0_0_pre2"
    process.dqmSaver.saveByRun = -1
    process.dqmSaver.saveAtJobEnd = True
    process.dqmSaver.forceRunNumber = 999999


    process.dqmOutput = cms.OutputModule("DQMRootOutputModule",
        fileName = cms.untracked.string("worker_output.root")
    )

    process.load("DQMServices.Components.EDMtoMEConverter_cff")
    process.end_path = cms.EndPath(process.dqmSaver)
else:
    process.dqmOutput = cms.OutputModule("DQMRootOutputModule", fileName=cms.untracked.string("run_output_shift_357734.root"))
    process.end_path = cms.EndPath(process.dqmOutput)


process.schedule = cms.Schedule(process.ALL,   process.end_path)

#print(process.dumpPython())
