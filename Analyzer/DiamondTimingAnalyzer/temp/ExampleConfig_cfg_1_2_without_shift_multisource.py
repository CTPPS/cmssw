import FWCore.ParameterSet.Config as cms
import copy
from FWCore.ParameterSet.VarParsing import VarParsing
from Geometry.VeryForwardGeometry.geometryRPFromDB_cfi import * # TODO it looks like it can be deleted


# TEMP config:
apply_shift = False                                        
saveToDQM = False


process = cms.Process('PPSTiming2') # or TIMINGSTUDY??
options = VarParsing ('analysis')

process.verbosity = cms.untracked.PSet( input = cms.untracked.int32(-1) ) # TODO: confirm if this is needed.

process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.threshold = ''
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

options.register ('rootInput',
		  '',
		  VarParsing.multiplicity.singleton,
                  VarParsing.varType.string,
                  "root input file name")


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
options.register('useDB', 
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
process.load('Geometry.VeryForwardGeometry.geometryRPFromDB_cfi') #TODO: use geometry form DB not from file 

process.source = cms.Source ("PoolSource",
                             fileNames = cms.untracked.vstring(options.rootInput)
                             )

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
# 2022 prompt: to be updated
process.GlobalTag = GlobalTag(process.GlobalTag, '124X_dataRun3_Prompt_frozen_v4', '')


#JH - use new tag for timing calibration

### SET TIME SHIFT SOURCE BEGIN
if options.useDB=='True':
    use_db=True
elif options.useDB=='False':
    use_db=False
elif options.useDB!='':
    assert 'UseDB paramter is not valid. It should be True or False (case sensitive)'
else: 
    use_db=False


if((options.sqlFileName != '') ^  (options.calibInput  != '') ^ (options.useDB !='')):
    assert 'Please specify exactly one source of time shift paramiters. One from set {calib.json, SQLLite File or globalDB}'

use_sqlite_file = options.sqlFileName != ''
if (use_sqlite_file):
        print('Using SQL file')                                        
        process.load('CondCore.CondDB.CondDB_cfi')
        process.CondDB.connect = options.sqlFileName # SQLite input TODO: migrate to using tag
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

if options.calibInput != '':
    print('Using CalibInput file ')
    process.ppsTimingCalibrationESSource = cms.ESSource('PPSTimingCalibrationESSource',
        calibrationFile = cms.string(options.calibInput),
        subDetector = cms.uint32(2),
        appendToDataLabel = cms.string('')
    )

if(use_db):
    print('Using db') 
    # TODO: uncomment below when delete sqlite file dependency 
    process.GlobalTag.toGet = cms.VPSet()
    process.GlobalTag.toGet.append(
    cms.PSet(record = cms.string("PPSTimingCalibrationRcd"),
            tag =  cms.string("PPSDiamondTimingCalibration_Run3_recovered_v1"),
            label = cms.untracked.string('PPSTestCalibration'),
            connect = cms.string("frontier://FrontierPrep/CMS_CONDITIONS")
        )
    )
    # process.GlobalTag.toGet = cms.VPSet(
    #     cms.PSet(record = cms.string('PPSTimingCalibrationRcd'),
    #                 tag = cms.string('PPSDiamondTimingCalibration_Run3_v1_hlt'), # working tag: PPSDiamondTimingCalibration_Run3_v1_hlt
    #                 #TODO: old tag PPSDiamondTimingCalibration_v1  - to delete
    #             connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS")
    #             )
    # )
### SET TIME SHIFT SOURCE END


# TODO: we miss tjhis from worker.py vvvvv
#  rechits production
process.load('RecoPPS.Local.ctppsDiamondLocalReconstruction_cff')
# ^^^^^^^^^^^^^^^^^^^^^^^^^

tagLocalTrack_ = cms.InputTag("ctppsDiamondLocalTracks","" ,"RECODQM")
tagRecHit_ =  cms.InputTag("ctppsDiamondRecHits","" ,"RECODQM")
from DQMServices.Core.DQMEDAnalyzer import DQMEDAnalyzer
if(options.calibInput != ''):
    print('Using calib input')
    process.diamondTimingWorker = DQMEDAnalyzer("DiamondTimingWorker",
        tagDigi = cms.InputTag("ctppsDiamondRawToDigiAlCaRecoProducer", "TimingDiamond"),
        tagRecHit = tagRecHit_, # changed
        tagPixelLocalTrack = cms.InputTag("ctppsPixelLocalTracksAlCaRecoProducer", "", "RECO"),
        tagLocalTrack = tagLocalTrack_, #changed
        timingCalibrationTag=cms.string(":"),
        tagValidOOT = cms.int32(-1), #TODO: remove parameter from options or don't hardcode it. 
        planesConfig = cms.string("planes.json"),
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
            tagValidOOT = cms.int32(-1), #TODO: remove parameter from options or don't hardcode it. 
            planesConfig = cms.string("planes.json"), #TODO: remove parameter from options or don't hardcode it. 
            Ntracks_Lcuts = cms.vint32([-1,1,-1,1]), # minimum number of tracks in pots [45-210, 45-220, 56-210, 56-220]
            Ntracks_Ucuts = cms.vint32([-1,6,-1,6]), # maximum number of tracks in pots [45-210, 45-220, 56-210, 56-220]
        )
elif(use_db):
    process.diamondTimingWorker = DQMEDAnalyzer("DiamondTimingWorker",
                tagDigi = cms.InputTag("ctppsDiamondRawToDigiAlCaRecoProducer", "TimingDiamond"),
                tagRecHit = tagRecHit_,
                tagPixelLocalTrack = cms.InputTag("ctppsPixelLocalTracksAlCaRecoProducer"),
                timingCalibrationTag=cms.string("GlobalTag:PPSTestCalibration"),
                tagLocalTrack =tagLocalTrack_,
                tagValidOOT = cms.int32(-1), #TODO: remove parameter from options or don't hardcode it. 
                planesConfig = cms.string("planes.json"), #TODO: remove parameter from options or don't hardcode it. 
                Ntracks_Lcuts = cms.vint32([-1,1,-1,1]), # minimum number of tracks in pots [45-210, 45-220, 56-210, 56-220]
                Ntracks_Ucuts = cms.vint32([-1,6,-1,6]), # maximum number of tracks in pots [45-210, 45-220, 56-210, 56-220]
            )
else: 
    assert "need to provide timing calibration tag from json, slq file or db"
process.ALL = cms.Path(process.diamondTimingWorker)


process.dqmOutput = cms.OutputModule("DQMRootOutputModule", fileName=cms.untracked.string(options.outputFileName))
process.end_path = cms.EndPath(process.dqmOutput)


process.schedule = cms.Schedule(process.ALL,   process.end_path)

#print(process.dumpPython())
