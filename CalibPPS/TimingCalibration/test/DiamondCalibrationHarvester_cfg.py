import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

from Configuration.AlCa.GlobalTag import GlobalTag
from Configuration.StandardSequences.Eras import eras


process = cms.Process("harvester", eras.Run3)

options = VarParsing.VarParsing()

options.register("globalTag",
    "",
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.string,
    "Global Tag"
)
options.register("inputFiles",
    "",
    VarParsing.VarParsing.multiplicity.list,
    VarParsing.VarParsing.varType.string
)
options.register("tVsLsFilename",
    "",
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.string,
    "ROOT filename with t vs LS histogram for double peak correction"
)

options.parseArguments()

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1_000_000

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.EventContent.EventContent_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag = GlobalTag(process.GlobalTag, options.globalTag, "")

# Source (histograms)
fileList = [
    f"file:{f}"
        if not (f.startswith("/store/") or f.startswith("file:") or f.startswith("root:"))
        else f for f in options.inputFiles
]

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(fileList)
)

# output service for database
process.load("CondCore.CondDB.CondDB_cfi")
process.CondDB.connect = "sqlite_file:ppsDiamondTiming_calibration.sqlite" # SQLite output

process.PoolDBOutputService = cms.Service("PoolDBOutputService",
    process.CondDB,
    timetype = cms.untracked.string("runnumber"),
    toPut = cms.VPSet(
        cms.PSet(
            record = cms.string("PPSTimingCalibrationRcd_HPTDC"),
            tag = cms.string("DiamondTimingCalibration")
        )
    )
)

process.load("CalibPPS.TimingCalibration.ppsTimingCalibrationPCLHarvester_cfi")
process.ppsTimingCalibrationPCLHarvester.tVsLsFilename = options.tVsLsFilename

# load DQM framework
process.load("DQMServices.Core.DQMStore_cfi")
process.load("DQMServices.Components.DQMEnvironment_cfi")
process.dqmEnv.subSystemFolder = "CalibPPS"
process.dqmSaver.convention = "Offline"
process.dqmSaver.workflow = "/CalibPPS/TimingCalibration/CMSSW_14_1_4"
process.dqmSaver.saveByRun = -1
process.dqmSaver.saveAtJobEnd = True

process.DQMStore = cms.Service("DQMStore")

process.dqmOutput = cms.OutputModule("DQMRootOutputModule",
    fileName = cms.untracked.string("harvester_output.root")
)

process.load("Geometry.VeryForwardGeometry.geometryRPFromDB_cfi")
process.load("DQMServices.Components.EDMtoMEConverter_cff")
process.EDMtoMEConverter.lumiInputTag = "MEtoEDMConvertPPSTimingCalib:MEtoEDMConverterLumi"
process.EDMtoMEConverter.runInputTag = "MEtoEDMConvertPPSTimingCalib:MEtoEDMConverterRun"

process.p = cms.Path(
    process.EDMtoMEConverter *
    process.ppsTimingCalibrationPCLHarvester
)

process.end_path = cms.EndPath(
    process.dqmEnv +
    process.dqmSaver
)

process.schedule = cms.Schedule(
    process.p,
    process.end_path
)
