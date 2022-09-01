
import FWCore.ParameterSet.Config as cms
import copy

process = cms.Process('PPSTiming2')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100000)
)

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = ''
process.MessageLogger.cerr.FwkReport.reportEvery = 1000




process.source = cms.Source ("PoolSource",
                             fileNames = cms.untracked.vstring(
    "file:/eos/home-m/mobrzut/357/734/00000/73c5eeec-578d-438f-92f1-3e1fff2ed390.root",
    "file:/eos/home-m/mobrzut/357/734/00000/08e5ee70-ebff-4571-b14c-806961a3d9dd.root",
    "file:/eos/home-m/mobrzut/357/734/00000/d42eaf2c-57bb-48fe-9a22-317053f89115.root"
)
)

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
# 2022 prompt: to be updated
process.GlobalTag = GlobalTag(process.GlobalTag, "124X_dataRun3_Prompt_v4") 


#JH - use new tag for timing calibrations



                                        
process.load('CondCore.CondDB.CondDB_cfi')
process.CondDB.connect = 'sqlite_file:/afs/cern.ch/user/m/mobrzut/automation/timing_calibration/CMSSW_12_4_6/src/Analyzer/DiamondTimingAnalyzer/temp/ppsDiamondTiming_calibration357734.sqlite' # SQLite input
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
                                        

# JH - rerun reco sequence with new timing conditions
process.load("RecoPPS.Configuration.recoCTPPS_cff")
process.ctppsDiamondLocalTracks.recHitsTag = cms.InputTag("ctppsDiamondRecHitsAlCaRecoProducer")
process.ctppsLocalTrackLiteProducer.tagDiamondTrack = cms.InputTag("ctppsDiamondLocalTracks","","PPSTiming2")
process.ctppsProtons.tagLocalTrackLite = cms.InputTag("ctppsLocalTrackLiteProducer","","PPSTiming2")
process.ctppsLocalTrackLiteProducer.includeDiamonds = cms.bool(True)
process.ctppsLocalTrackLiteProducer.tagPixelTrack = cms.InputTag("ctppsPixelLocalTracksAlCaRecoProducer")
process.ctppsLocalTrackLiteProducer.includePixels = cms.bool(True)
process.ctppsDiamondRecHits.timingCalibrationTag=cms.string("PoolDBESSource:PPSTestCalibration")
process.ctppsDiamondRecHits.digiTag= cms.InputTag("ctppsDiamondRawToDigiAlCaRecoProducer:TimingDiamond")

from DQMServices.Core.DQMEDAnalyzer import DQMEDAnalyzer
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



process.ALL = cms.Path(
    # Re-run the PPS local+proton reconstruction starting from AOD with new timing calibrations
    process.ctppsDiamondRecHits *
    process.ctppsDiamondLocalTracks *
    process.ctppsLocalTrackLiteProducer *
    process.ctppsProtons *
    process.diamondTimingWorker
                       )

process.dqmOutput = cms.OutputModule("DQMRootOutputModule", fileName=cms.untracked.string("temp_run_output.root"))

process.end_path = cms.EndPath(process.dqmOutput)


process.schedule = cms.Schedule(process.ALL)

#print(process.dumpPython())
