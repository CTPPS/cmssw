import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
from RecoPPS.Local.ctppsDiamondLocalReconstruction_cff import *

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

options.register ('geometryFile',
				  'Geometry.VeryForwardGeometry.geometryRPFromDD_2022_cfi',
				  VarParsing.multiplicity.singleton,
				  VarParsing.varType.string,
				  "Geometry input file")

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

if options.calibInput == '':
    process.GlobalTag.toGet = cms.VPSet(
        cms.PSet(record = cms.string('PPSTimingCalibrationRcd'),
                    tag = cms.string('PPSDiamondTimingCalibration_Run3_v1_hlt'),
                connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS")
                )
    )
else:
	process.ppsTimingCalibrationESSource = cms.ESSource('PPSTimingCalibrationESSource',
		calibrationFile = cms.string(options.calibInput),
		subDetector = cms.uint32(2),
		appendToDataLabel = cms.string('')
	)

# rechits production
process.load("DQM.Integration.config.environment_cfi")
process.load("DQMServices.Components.DQMEnvironment_cfi")
process.load(options.geometryFile)

process.diamondTimingHarvester = DQMEDHarvester("DiamondTimingHarvester",
   calib_json_output = cms.string(options.calibOutput),
   calibInput = cms.string(options.calibInput),
   calibFiles = cms.vstring(options.calibFiles),
   loopIndex = cms.int32(options.loopIndex),
   meanMax = cms.double(options.meanMax),
   rmsMax = cms.double(options.rmsMax)
)

#CONFIGURE DQM Saver
process.dqmEnv.subSystemFolder = "CalibPPS"
process.dqmSaver.convention = 'Offline'
process.dqmSaver.workflow = "/AnalyzerWithShift9/DiamondTimingAnalyzer/CMSSW_12_4_6"
process.dqmSaver.dirName = './OutputFiles/' # TODO: it was commented out - revert the change before merge

process.path = cms.Path(
   process.diamondTimingHarvester
)

process.end_path = cms.EndPath(process.dqmSaver)

process.schedule = cms.Schedule(
    process.path,
    process.end_path
)
