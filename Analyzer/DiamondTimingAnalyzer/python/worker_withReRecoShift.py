import FWCore.ParameterSet.Config as cms
import copy
from FWCore.ParameterSet.VarParsing import VarParsing
from Geometry.VeryForwardGeometry.geometryRPFromDB_cfi import * # TODO it looks like it can be deleted


process = cms.Process('TIMINGSTUDY') # or TIMINGSTUDY??
options = VarParsing ('analysis')

process.verbosity = cms.untracked.PSet( input = cms.untracked.int32(-1) ) # TODO: confirm if this is needed.

process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.threshold = ''
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

options.register('globalTag',
                '',
                VarParsing.multiplicity.singleton,
                VarParsing.varType.string,
                "GT to use")

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
				  '',
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

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100000))

if len(options.inputFiles) != 0:
    fileList = [f'file:{f}' if not (f.startswith('/store/') or f.startswith('file:')) else f for f in options.inputFiles]
    inputFiles = cms.untracked.vstring(fileList)
    print('Input files:')
    print(inputFiles)
else:
    fileList = ["/store/data/Run2023D/AlCaPPSPrompt/ALCARECO/PPSCalMaxTracks-PromptReco-v1/000/369/956/00000/36ea1bb9-f031-41da-a7a2-3bd258d0e255.root",
          "/store/data/Run2023D/AlCaPPSPrompt/ALCARECO/PPSCalMaxTracks-PromptReco-v1/000/369/956/00000/8d7b7d9a-a891-4027-aeeb-e198037e1a0b.root",
          "/store/data/Run2023D/AlCaPPSPrompt/ALCARECO/PPSCalMaxTracks-PromptReco-v1/000/369/956/00000/a0f8f0b6-ee2b-4b0d-aa59-72ab1d74235d.root",
          "/store/data/Run2023D/AlCaPPSPrompt/ALCARECO/PPSCalMaxTracks-PromptReco-v1/000/369/956/00000/79c30c13-96bf-4193-a432-48e83616b553.root",
          "/store/data/Run2023D/AlCaPPSPrompt/ALCARECO/PPSCalMaxTracks-PromptReco-v1/000/369/956/00000/3f60fd90-b91c-44d3-a4cb-e667cebf1137.root",
          "/store/data/Run2023D/AlCaPPSPrompt/ALCARECO/PPSCalMaxTracks-PromptReco-v1/000/369/956/00000/aefb71f3-3aa3-4fcd-bdd9-513e3620a3c8.root",
          "/store/data/Run2023D/AlCaPPSPrompt/ALCARECO/PPSCalMaxTracks-PromptReco-v1/000/369/956/00000/84e9108d-c99d-4bb6-ae75-59992ab55e02.root",
          "/store/data/Run2023D/AlCaPPSPrompt/ALCARECO/PPSCalMaxTracks-PromptReco-v1/000/369/956/00000/0f00e63b-b735-42d4-b3b1-13184cae5f1e.root"]  # test file run 369956
    inputFiles = cms.untracked.vstring( *fileList)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Geometry.VeryForwardGeometry.geometryRPFromDB_cfi') #TODO: use geometry form DB not from file 

print("Input files:", inputFiles)

process.source = cms.Source ("PoolSource",
                             fileNames = inputFiles
                             )

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
# 2022 prompt: to be updated

if options.globalTag != '':
    gt = options.globalTag
else:
    gt = '130X_dataRun3_Prompt_v4'

print('Using GT:',gt)

process.GlobalTag = GlobalTag(process.GlobalTag, gt, '')

# WHEN CHANGING THIS CHANGE ALSO THE LINES AT ELSE STATEMENT AROUND 153-156
process.GlobalTag.toGet.append(
  cms.PSet(record = cms.string("PPSTimingCalibrationRcd"),
           tag =  cms.string("CTPPPSTimingCalibration_HPTDC_byPCL_v1_prompt"),
           label = cms.untracked.string('PPSTimingCalibrationbyPCL'),
           connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS")
	)
)




#JH - use new tag for timing calibration

### SET TIME SHIFT SOURCE BEGIN
# if options.useDB=='True':
#     use_db=True
# elif options.useDB=='False':
#     use_db=False
# elif options.useDB!='':
#     assert 'UseDB paramter is not valid. It should be True or False (case sensitive)'
# else: 
#     use_db=False


# if((options.sqlFileName != '') ^  (options.calibInput  != '') ^ (options.useDB !='')):
    # assert 'Please specify exactly one source of time shift paramiters. One from set {calib.json, SQLLite File or globalDB}'

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

elif options.calibInput != '':
    print('Using CalibInput file ')
    process.ppsTimingCalibrationESSource = cms.ESSource('PPSTimingCalibrationESSource',
        calibrationFile = cms.string(options.calibInput),
        subDetector = cms.uint32(2),
        appendToDataLabel = cms.string('')
    )

else: # default use db
    print('Using db')
    # TODO: uncomment below when delete sqlite file dependency 
    process.GlobalTag.toGet = cms.VPSet()
    process.GlobalTag.toGet.append(
    cms.PSet(record = cms.string("PPSTimingCalibrationRcd"),
            tag =  cms.string("CTPPPSTimingCalibration_HPTDC_byPCL_v1_prompt"),
            label = cms.untracked.string('PPSTimingCalibrationbyPCL'),
            connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS")
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
process.ctppsDiamondRecHits.timingCalibrationTag="GlobalTag:PPSTimingCalibrationbyPCL"
process.ctppsDiamondRecHits.digiTag="ctppsDiamondRawToDigiAlCaRecoProducer:TimingDiamond"
# ^^^^^^^^^^^^^^^^^^^^^^^^^

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
        tagValidOOT = cms.int32(-1), #TODO: remove parameter from options or don't hardcode it. 
        planesConfig = cms.string("/afs/cern.ch/user/m/molkowsk/work/CMSSW_13_0_9/src/Analyzer/DiamondTimingAnalyzer/planes.json"),
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
            planesConfig = cms.string("/afs/cern.ch/user/m/molkowsk/work/CMSSW_13_0_9/src/Analyzer/DiamondTimingAnalyzer/planes.json"), #TODO: remove parameter from options or don't hardcode it. 
            Ntracks_Lcuts = cms.vint32([-1,1,-1,1]), # minimum number of tracks in pots [45-210, 45-220, 56-210, 56-220]
            Ntracks_Ucuts = cms.vint32([-1,6,-1,6]), # maximum number of tracks in pots [45-210, 45-220, 56-210, 56-220]
        )
else:
    process.diamondTimingWorker = DQMEDAnalyzer("DiamondTimingWorker",
                tagDigi = cms.InputTag("ctppsDiamondRawToDigiAlCaRecoProducer", "TimingDiamond"),
                tagRecHit = tagRecHit_,
                tagPixelLocalTrack = cms.InputTag("ctppsPixelLocalTracksAlCaRecoProducer"),
                timingCalibrationTag=cms.string("GlobalTag:PPSTimingCalibrationbyPCL"),
                tagLocalTrack =tagLocalTrack_,
                tagValidOOT = cms.int32(-1), #TODO: remove parameter from options or don't hardcode it. 
                planesConfig = cms.string("/afs/cern.ch/user/m/molkowsk/work/CMSSW_13_0_9/src/Analyzer/DiamondTimingAnalyzer/planes.json"), #TODO: remove parameter from options or don't hardcode it. 
                Ntracks_Lcuts = cms.vint32([-1,1,-1,1]), # minimum number of tracks in pots [45-210, 45-220, 56-210, 56-220]
                Ntracks_Ucuts = cms.vint32([-1,6,-1,6]), # maximum number of tracks in pots [45-210, 45-220, 56-210, 56-220]
            )
# else: 
    # assert "need to provide timing calibration tag from json, slq file or db"
process.content = cms.EDAnalyzer("EventContentAnalyzer") 
process.ALL = cms.Path(
    process.ctppsDiamondLocalReconstruction *
    #process.content*
    process.diamondTimingWorker    
    ) 


process.dqmOutput = cms.OutputModule("DQMRootOutputModule", fileName=cms.untracked.string(options.outputFileName))
process.end_path = cms.EndPath(process.dqmOutput)

process.schedule = cms.Schedule(process.ALL,   process.end_path)

#print(process.dumpPython())
