import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

process = cms.Process('test')

options = VarParsing ('analysis')
options.register ('calibInput',
				  '',
				  VarParsing.multiplicity.singleton,
				  VarParsing.varType.string,
				  "Calibration input file for this iteration")
options.register ('outputName',
				  '',
				  VarParsing.multiplicity.singleton,
				  VarParsing.varType.string,
				  "calib.sqlite")
options.register ('IOVFirst',
				  1,
				  VarParsing.multiplicity.singleton,
				  VarParsing.varType.int,
				  "")
options.register ('IOVLast',
				  1,
				  VarParsing.multiplicity.singleton,
				  VarParsing.varType.int,
				  "")
options.register ('IOVInterval',
				  1,
				  VarParsing.multiplicity.singleton,
				  VarParsing.varType.int,
				  "")
options.parseArguments()

from CondFormats.PPSObjects.PPSTimingDetEnum_cff import PPSTimingDetEnum

process.source = cms.Source('EmptyIOVSource',
    timetype = cms.string('runnumber'),
    firstValue = cms.uint64(options.IOVFirst),
    lastValue = cms.uint64(options.IOVLast),
    interval = cms.uint64(options.IOVInterval)
)

# load calibrations from JSON file
process.load('CalibPPS.ESProducers.ppsTimingCalibrationESSource_cfi')
process.ppsTimingCalibrationESSource.calibrationFile = cms.string(options.calibInput)
process.ppsTimingCalibrationESSource.subDetector = PPSTimingDetEnum.PPS_DIAMOND

# output service for database
process.load('CondCore.CondDB.CondDB_cfi')
process.CondDB.connect = 'sqlite_file:'+options.outputName # SQLite output

process.PoolDBOutputService = cms.Service('PoolDBOutputService',
    process.CondDB,
    timetype = cms.untracked.string('runnumber'),
    toPut = cms.VPSet(
        cms.PSet(
            record = cms.string('PPSTimingCalibrationRcd'),
            tag = cms.string('PPSDiamondTimingCalibration'),
        )
    )
)

process.ppsTimingCalibrationWriter = cms.EDAnalyzer('PPSTimingCalibrationWriter')

process.path = cms.Path(
    process.ppsTimingCalibrationWriter
)