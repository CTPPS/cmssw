import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
import FWCore.PythonUtilities.LumiList as LumiList

from Configuration.AlCa.GlobalTag import GlobalTag


process = cms.Process("worker")

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
options.register("jsonFileName",
    "",
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.string,
    "Certification JSON filename"
)

options.parseArguments()

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1_000_000

process.load("RecoPPS.Local.totemTimingLocalReconstruction_cff")
process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.a1 = cms.EDAnalyzer("StreamThingAnalyzer",
    product_to_get = cms.string("m1")
)

process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.EventContent.EventContent_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag = GlobalTag(process.GlobalTag, options.globalTag, "")
process.load("RecoPPS.Configuration.recoCTPPS_cff")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(options.inputFiles)
)

if options.jsonFileName != "":
    process.source.lumisToProcess = LumiList.LumiList(filename=options.jsonFileName).getVLuminosityBlockRange()
    print(f"Using JSON file: {options.jsonFileName}")

process.load("CalibPPS.TimingCalibration.ppsTimingCalibrationPCLWorker_cfi")
process.ppsTimingCalibrationPCLWorker.tVsLsFilename = options.tVsLsFilename

process.DQMStore = cms.Service("DQMStore")

process.dqmOutput = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string("file:worker_output.root"),
    outputCommands = cms.untracked.vstring(
        "drop *",
        "keep *_MEtoEDMConvertPPSTimingCalib_*_*"
    )
)

process.load("CalibPPS.TimingCalibration.ALCARECOPromptCalibProdPPSTimingCalib_cff")

process.path = cms.Path(
    process.ppsTimingCalibrationPCLWorker *
    process.MEtoEDMConvertPPSTimingCalib
)

process.end_path = cms.EndPath(
    process.dqmOutput
)

process.schedule = cms.Schedule(
    process.path,
    process.end_path
)
