import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing


process = cms.Process('test')
options = VarParsing()

options.register('runNumber',
                     -1,
                      VarParsing.multiplicity.singleton,
                      VarParsing.varType.int,
                      "Run number")

options.register('sqlitePath',
                     'ppsDiamondTiming_calibration.sqlite',
                      VarParsing.multiplicity.singleton,
                      VarParsing.varType.string,
                      "Path to the SQLite file")

options.register('record',
                     'PPSTimingCalibrationRcd',
                      VarParsing.multiplicity.singleton,
                      VarParsing.varType.string,
                      "Record")

options.register('tag',
                     'PPSDiamondTimingCalibration',
                      VarParsing.multiplicity.singleton,
                      VarParsing.varType.string,
                      "Tag")

options.parseArguments()

# minimum logging
process.MessageLogger = cms.Service("MessageLogger",
    cerr = cms.untracked.PSet(
        enable = cms.untracked.bool(False)
    ),
    cout = cms.untracked.PSet(
        enable = cms.untracked.bool(True),
        threshold = cms.untracked.string('INFO')
    )
)

process.source = cms.Source('EmptyIOVSource',
    timetype = cms.string('runnumber'),
    firstValue = cms.uint64(options.runNumber),
    lastValue = cms.uint64(options.runNumber),
    interval = cms.uint64(1)
)

# load calibrations from database
process.load('CondCore.CondDB.CondDB_cfi')
process.CondDB.connect = f'sqlite_file:{options.sqlitePath}' # SQLite input

process.PoolDBESSource = cms.ESSource('PoolDBESSource',
    process.CondDB,
    DumpStat = cms.untracked.bool(True),
    toGet = cms.VPSet(
        cms.PSet(
            record = cms.string(options.record),
            tag = cms.string(options.tag)
        )
    )
)

process.ppsTimingCalibrationAnalyzer = cms.EDAnalyzer('PPSTimingCalibrationAnalyzer')

process.path = cms.Path(
    process.ppsTimingCalibrationAnalyzer
)

