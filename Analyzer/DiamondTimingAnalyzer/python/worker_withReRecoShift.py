import FWCore.ParameterSet.Config as cms
import copy
from FWCore.ParameterSet.VarParsing import VarParsing
from Geometry.VeryForwardGeometry.geometryRPFromDB_cfi import * # TODO it looks like it can be deleted
from Configuration.Eras.Era_Run3_cff import Run3


process = cms.Process('TIMINGSTUDY', Run3) 
options = VarParsing ('analysis')

process.verbosity = cms.untracked.PSet( input = cms.untracked.int32(-1) ) # TODO: confirm if this is needed.

process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.threshold = ''
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

options.register ('outputFileName',
		  'run_output.root',
		  VarParsing.multiplicity.singleton,
                  VarParsing.varType.string,
                  "output file name")

options.register ('calibInput',
				  '',
				  VarParsing.multiplicity.singleton,
				  VarParsing.varType.string,
				  "Calibration input file for this iteration")

options.register ('planesConfig',
				  'planes.json',
				  VarParsing.multiplicity.singleton,
				  VarParsing.varType.string,
				  "planes config file")

options.register ('validOOT',
				  -1,
				  VarParsing.multiplicity.singleton,
				  VarParsing.varType.int,
				  "valid OOT slice")	
options.register('useDB',  # default source of time shift tag
                    '',
				  VarParsing.multiplicity.singleton,
                  VarParsing.varType.string,
                  'useDB'
)

options.register('sqlFileName', 
                    '',
				  VarParsing.multiplicity.singleton,
                  VarParsing.varType.string,
                  'sqlFileNameForCalib'
)

options.register('maxEventsToProcess', 
                    -1,
				  VarParsing.multiplicity.singleton,
                  VarParsing.varType.int,
                  'maxEvents to process'
)		  
options.parseArguments()

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(options.maxEventsToProcess))

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Geometry.VeryForwardGeometry.geometryRPFromDB_cfi')

process.source = cms.Source ("PoolSource",
    fileNames = cms.untracked.vstring(
        #  RUN 440 - "357440": [[1,354]],
        "/store/data/Run2022C/AlCaPPSPrompt/ALCARECO/PPSCalMaxTracks-PromptReco-v1/000/357/440/00000/53746080-d21b-4fdf-9d19-67c6dae347ac.root",
        "/store/data/Run2022C/AlCaPPSPrompt/ALCARECO/PPSCalMaxTracks-PromptReco-v1/000/357/440/00000/b8d3f012-0b59-43b4-adbc-811bcb59c9c4.root",
        "/store/data/Run2022C/AlCaPPSPrompt/ALCARECO/PPSCalMaxTracks-PromptReco-v1/000/357/440/00000/cfa9f82a-5296-49cf-8353-11f8761b675c.root",
        "/store/data/Run2022C/AlCaPPSPrompt/ALCARECO/PPSCalMaxTracks-PromptReco-v1/000/357/440/00000/db3c6c33-9ad0-4d61-9a64-a1aa9578550e.root"    
    )
)

# Manipulate GT to pick up the right time shift
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
# 2022 prompt: to be updated
process.GlobalTag = GlobalTag(process.GlobalTag, '124X_dataRun3_Prompt_frozen_v4', '')

use_sqlite_file = options.sqlFileName != ''
if (use_sqlite_file):
    print('Using SQL file')                                        
    process.load('CondCore.CondDB.CondDB_cfi')
    process.CondDB.connect = options.sqlFileName 
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
elif options.calibInput != '':
    print('Using CalibInput file ')
    process.ppsTimingCalibrationESSource = cms.ESSource('PPSTimingCalibrationESSource',
        calibrationFile = cms.string(options.calibInput),
        subDetector = cms.uint32(2),
        appendToDataLabel = cms.string('')
    )

else: # default use db
    print('Using db') 
    process.GlobalTag.toGet = cms.VPSet()
    process.GlobalTag.toGet.append(
    cms.PSet(record = cms.string("PPSTimingCalibrationRcd"),
            tag =  cms.string("PPSDiamondTimingCalibration_Run3_recovered_v1"),
            label = cms.untracked.string('PPSTestCalibration'),
            connect = cms.string("frontier://FrontierPrep/CMS_CONDITIONS")
        )
    )

process.load("RecoPPS.Configuration.recoCTPPS_cff")
process.ctppsDiamondRecHits.digiTag=cms.InputTag('ctppsDiamondRawToDigiAlCaRecoProducer:TimingDiamond')
process.ctppsDiamondRecHits.timingCalibrationTag=cms.string("GlobalTag:PPSTestCalibration")

# Configure rolling calibration worker
tagLocalTrack_ = cms.InputTag("ctppsDiamondLocalTracks","" ,"TIMINGSTUDY")
tagRecHit_ =  cms.InputTag("ctppsDiamondRecHits","" ,"TIMINGSTUDY")

from DQMServices.Core.DQMEDAnalyzer import DQMEDAnalyzer
if(options.calibInput != ''):
    print('Using calib input')
    process.diamondTimingWorker = DQMEDAnalyzer("DiamondTimingWorker",
        tagDigi = cms.InputTag("ctppsDiamondRawToDigiAlCaRecoProducer", "TimingDiamond"),
        tagRecHit = tagRecHit_, # changed
        tagPixelLocalTrack = cms.InputTag("ctppsPixelLocalTracksAlCaRecoProducer", "", "RECO"),
        tagLocalTrack = tagLocalTrack_, #changed
        timingCalibrationTag=cms.string(":"),
        tagValidOOT = cms.int32(options.validOOT), 
        planesConfig = cms.string(options.planesConfig),
        Ntracks_Lcuts = cms.vint32([-1,1,-1,1]), # minimum number of tracks in pots [45-210, 45-220, 56-210, 56-220]
        Ntracks_Ucuts = cms.vint32([-1,6,-1,6]), # maximum number of tracks in pots [45-210, 45-220, 56-210, 56-220]
    ) 
elif (use_sqlite_file):
    process.diamondTimingWorker = DQMEDAnalyzer("DiamondTimingWorker",
            tagDigi = cms.InputTag("ctppsDiamondRawToDigiAlCaRecoProducer", "TimingDiamond"),
            tagRecHit =tagRecHit_,
            tagPixelLocalTrack = cms.InputTag("ctppsPixelLocalTracksAlCaRecoProducer"),
            timingCalibrationTag=cms.string("PoolDBESSource:PPSTestCalibration"),
            tagLocalTrack = tagLocalTrack_,
            tagValidOOT = cms.int32(options.validOOT, 
            planesConfig = cms.string(options.planesConfig), 
            Ntracks_Lcuts = cms.vint32([-1,1,-1,1]), # minimum number of tracks in pots [45-210, 45-220, 56-210, 56-220]
            Ntracks_Ucuts = cms.vint32([-1,6,-1,6]), # maximum number of tracks in pots [45-210, 45-220, 56-210, 56-220]
        )
else:
    process.diamondTimingWorker = DQMEDAnalyzer("DiamondTimingWorker",
                tagDigi = cms.InputTag("ctppsDiamondRawToDigiAlCaRecoProducer", "TimingDiamond",""),
                tagRecHit = tagRecHit_,
                tagPixelLocalTrack = cms.InputTag("ctppsPixelLocalTracksAlCaRecoProducer"),
                timingCalibrationTag=cms.string("GlobalTag:PPSTestCalibration"),
                tagLocalTrack =tagLocalTrack_,
                tagValidOOT = cms.int32(-1), 
                planesConfig = cms.string(options.planesConfig),  
                Ntracks_Lcuts = cms.vint32([-1,1,-1,1]), # minimum number of tracks in pots [45-210, 45-220, 56-210, 56-220]
                Ntracks_Ucuts = cms.vint32([-1,6,-1,6]), # maximum number of tracks in pots [45-210, 45-220, 56-210, 56-220]
            )

process.ALL = cms.Path(
    process.ctppsDiamondRecHits *
    process.ctppsDiamondLocalTracks *
    process.diamondTimingWorker *
)

process.dqmOutput = cms.OutputModule("DQMRootOutputModule", fileName=cms.untracked.string(options.outputFileName))
process.end_path = cms.EndPath(process.dqmOutput)

process.schedule = cms.Schedule(process.ALL,   process.end_path)
