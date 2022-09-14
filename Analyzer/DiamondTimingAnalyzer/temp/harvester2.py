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
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

options.register ('rootFiles',
				  '', 
                  VarParsing.multiplicity.singleton, #TODO: it was designed to be a list revert change before merging
				  VarParsing.varType.string,
				  "root files produced by DQMWorker")

options.register ('outputDirectoryRoot',
				  'OutputFilesRolling', 
                  VarParsing.multiplicity.singleton, #TODO: it was designed to be a list revert change before merging
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
process.GlobalTag = GlobalTag(process.GlobalTag, '124X_dataRun3_v9', '')

useSQLFile = True
if useSQLFile:
	process.load('CondCore.CondDB.CondDB_cfi')
	process.CondDB.connect =  options.sqlFileName # SQLite input
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


# if options.calibInput == '':
#     process.PoolDBESSource = cms.ESSource('PoolDBESSource',
#         process.CondDB,
#         DumpStats = cms.untracked.bool(True),
#         toGet = cms.VPSet(
#             cms.PSet(
#                 record = cms.string('PPSTimingCalibrationRcd'),
#                 tag = cms.string('DiamondTimingCalibration'),
# 				label = cms.untracked.string('PPSTestCalibration')
#             )
#         )
#     )
# else:
if options.calibInput != '':
	process.ppsTimingCalibrationESSource = cms.ESSource('PPSTimingCalibrationESSource',
		calibrationFile = cms.string(options.calibInput),
		subDetector = cms.uint32(2),
		appendToDataLabel = cms.string('')
	)	

# rechits production
process.load("DQM.Integration.config.environment_cfi")
process.load("DQMServices.Components.DQMEnvironment_cfi")
process.load('Geometry.VeryForwardGeometry.geometryRPFromDB_cfi') #TODO: use geometry form DB not from file 

process.diamondTimingHarvester = DQMEDHarvester("DiamondTimingHarvester",
#    timingCalibrationTag=cms.string("PoolDBESSource:PPSTestCalibration"),
#    timingCalibrationTag=cms.InputTag("PoolDBESSource:PPSTestCalibration"),

   calib_json_output = cms.string(options.calibOutput),
   calibFiles = cms.vstring(options.calibFiles),
   loopIndex = cms.int32(options.loopIndex),
   meanMax = cms.double(options.meanMax),
   rmsMax = cms.double(options.rmsMax)
)

# process.diamondTimingHarvester.timingCalibrationTag=cms.string("PoolDBESSource:PPSTestCalibration")


#CONFIGURE DQM Saver
process.dqmEnv.subSystemFolder = "CalibPPS"
process.dqmSaver.convention = 'Offline'
process.dqmSaver.workflow = f'/AnalyzerWithShiftRolling/DiamondTimingAnalyzer/CMSSW_12_4_6_loop{options.loopIndex}'
process.dqmSaver.dirName = f'./{options.outputDirectoryRoot}/' # TODO: it was commented out - revert the change before merge

process.path = cms.Path(
   process.diamondTimingHarvester
)

process.end_path = cms.EndPath(process.dqmSaver)

process.schedule = cms.Schedule(
    process.path,
    process.end_path
)
