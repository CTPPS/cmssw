import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
from RecoPPS.Local.ctppsDiamondLocalReconstruction_cff import *
from Geometry.VeryForwardGeometry.geometryRPFromDB_cfi import *


process = cms.Process("TIMINGSTUDY")
options = VarParsing ('analysis')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.verbosity = cms.untracked.PSet( input = cms.untracked.int32(-1) )

# minimum of logs
process.load("FWCore.MessageService.MessageLogger_cfi")
# process.MessageLogger.cerr.threshold = "INFO"
# process.MessageLogger.cerr.INFO.limit = -1
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

options.register('globalTag',
                '',
                VarParsing.multiplicity.singleton,
                VarParsing.varType.string,
                "GT to use")

options.register ('rootFiles',
				  'file:run_output.root', 
                  VarParsing.multiplicity.list,
				  VarParsing.varType.string,
				  "root files produced by DQMWorker")

options.register ('outputDirectoryRoot',
				  'OutputHarvester', 
                  VarParsing.multiplicity.singleton,
				  VarParsing.varType.string,
				  "root output dircetory")

# options.register ('calibInput',
# 				  '',
# 				  VarParsing.multiplicity.singleton,
# 				  VarParsing.varType.string,
# 				  "Input file for calibration from this iteration")

options.register ('calibOutput',
				  'calib.json',
				  VarParsing.multiplicity.singleton,
				  VarParsing.varType.string,
				  "Output file for calibration from this iteration")

options.register ('calibInput',
				  '',
				  VarParsing.multiplicity.singleton,
				  VarParsing.varType.string,
				  "Input file")

options.register ('calibFiles',
				  [],
				  VarParsing.multiplicity.list,
				  VarParsing.varType.string,
				  "Calibration files from previous iterations in ascending order")

options.register ('loopIndex',
				  0,
				  VarParsing.multiplicity.singleton,
				  VarParsing.varType.int,
				  "loop iteration")

options.register ('meanMax',
				  1,
				  VarParsing.multiplicity.singleton,
				  VarParsing.varType.float,
				  "max mean of t distribution")

options.register ('rmsMax',
				  1,
				  VarParsing.multiplicity.singleton,
				  VarParsing.varType.float,
				  "max rms of t distribution")

options.register('sqlFileName', 
                    '',
				  VarParsing.multiplicity.singleton,
                  VarParsing.varType.string,
                  'sqlFileNameForCalib'
)

options.register('useDB', 
                    '',
				  VarParsing.multiplicity.singleton,
                  VarParsing.varType.string,
                  'use data base for calibraiton'
)
		
				  		  
options.parseArguments()
				  

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

#SOURCE
process.source = cms.Source("DQMRootSource",
    fileNames = cms.untracked.vstring(options.rootFiles)
)

from DQMServices.Core.DQMEDHarvester import DQMEDHarvester
from Configuration.AlCa.GlobalTag import GlobalTag

if options.globalTag != '':
    gt = options.globalTag
else:
    gt = '130X_dataRun3_Prompt_v4'

process.GlobalTag = GlobalTag(process.GlobalTag, gt, '')

### SET TIME SHIFT SOURCE BEGIN
if options.useDB=='True':
    use_db=True
elif options.useDB=='False':
    use_db=False
elif options.useDB!='':
    assert 'UseDB paramter is not valid. It should be True or False (case sensitive)'
else: 
    use_db=False


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
else: #default use db
    print('Using db') 
    process.GlobalTag.toGet = cms.VPSet()
    process.GlobalTag.toGet.append(
    cms.PSet(record = cms.string("PPSTimingCalibrationRcd"),
            tag =  cms.string("PPSDiamondTimingCalibration_Run3_recovered_v1"),
            label = cms.untracked.string('PPSTestCalibration'),
            connect = cms.string("frontier://FrontierPrep/CMS_CONDITIONS")
        )
    )
    #TODO: uncomment below when delete sqlite file dependency 
    # process.GlobalTag.toGet = cms.VPSet(
    #     cms.PSet(record = cms.string('PPSTimingCalibrationRcd'),
    #                 tag = cms.string('PPSDiamondTimingCalibration_Run3_v1_hlt'), # working tag: PPSDiamondTimingCalibration_Run3_v1_hlt
    #                 #TODO: old tag PPSDiamondTimingCalibration_v1  - to delete
    #             connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS")
    #             )
    # )
### SET TIME SHIFT SOURCE END

# rechits production
process.load("DQM.Integration.config.environment_cfi")
process.load("DQMServices.Components.DQMEnvironment_cfi")
process.load('Geometry.VeryForwardGeometry.geometryRPFromDB_cfi') #TODO: use geometry form DB not from file 

# PROCESS HARVESTER:
if(options.calibInput != ''):
	process.diamondTimingHarvester = DQMEDHarvester("DiamondTimingHarvester",
	#    timingCalibrationTag=cms.string("PoolDBESSource:PPSTestCalibration"),
	#    timingCalibrationTag=cms.string("PPSTimingCalibrationESSource:PPSTiming"),
	# timingCalibrationTag=cms.string("PoolDBESSource:PPSTestCalibration"),
	timingCalibrationTag=cms.string(":"),
	calib_json_output = cms.string(options.calibOutput),
	calibFiles = cms.vstring(options.calibFiles),
	loopIndex = cms.int32(options.loopIndex),
	meanMax = cms.double(options.meanMax),
	rmsMax = cms.double(options.rmsMax)
	)
elif (use_sqlite_file):
	process.diamondTimingHarvester = DQMEDHarvester("DiamondTimingHarvester",
	timingCalibrationTag=cms.string("PoolDBESSource:PPSTestCalibration"),
   calib_json_output = cms.string(options.calibOutput),
   calibFiles = cms.vstring(options.calibFiles),
   loopIndex = cms.int32(options.loopIndex),
   meanMax = cms.double(options.meanMax),
   rmsMax = cms.double(options.rmsMax)
)
else: # defualt use db
    process.diamondTimingHarvester = DQMEDHarvester("DiamondTimingHarvester",
    timingCalibrationTag=cms.string("GlobalTag:PPSTestCalibration"),
	calib_json_output = cms.string(options.calibOutput),
	calibFiles = cms.vstring(options.calibFiles),
	loopIndex = cms.int32(options.loopIndex),
	meanMax = cms.double(options.meanMax),
	rmsMax = cms.double(options.rmsMax)
	)
# else: 
#     assert "need to provide timing calibration tag from json, slq file or db"

# process.diamondTimingHarvester.timingCalibrationTag=cms.string("PoolDBESSource:PPSTestCalibration")


#CONFIGURE DQM Saver
process.dqmEnv.subSystemFolder = "CalibPPS"
process.dqmSaver.convention = 'Offline'
process.dqmSaver.workflow = "/AnalyzerWithShift9/DiamondTimingAnalyzer/CMSSW_12_4_6"
# process.dqmSaver.dirName = './OutputFiles/' # TODO: it was commented out - revert the change before merge

process.path = cms.Path(
   process.diamondTimingHarvester
)

process.end_path = cms.EndPath(process.dqmSaver)

process.schedule = cms.Schedule(
    process.path,
    process.end_path
)
