import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

process = cms.Process("worker")

options = VarParsing.VarParsing ()
options.register('globalTag',
                     '', 
                      VarParsing.VarParsing.multiplicity.singleton,
                      VarParsing.VarParsing.varType.string,
                      "Global Tag")
                      
options.register('inputFiles',
                 '',
                 VarParsing.VarParsing.multiplicity.list,
                 VarParsing.VarParsing.varType.string)
options.parseArguments()


process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('RecoPPS.Local.totemTimingLocalReconstruction_cff')
process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.a1 = cms.EDAnalyzer("StreamThingAnalyzer",
    product_to_get = cms.string('m1')
)

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
from Configuration.AlCa.autoCond import autoCond
#process.GlobalTag = GlobalTag(process.GlobalTag, autoCond['run3_data_prompt'], '')
process.GlobalTag = GlobalTag(process.GlobalTag, options.globalTag,'')
#process.load("EventFilter.CTPPSRawToDigi.ctppsRawToDigi_cff")
process.load("RecoPPS.Configuration.recoCTPPS_cff")

#process.load('CondCore.CondDB.CondDB_cfi')
#process.CondDB.connect = 'sqlite_file:ppsDiamondTiming_calibration.sqlite' # SQLite input
#process.PoolDBESSource = cms.ESSource('PoolDBESSource',
#        process.CondDB,
#        DumpStat = cms.untracked.bool(True),
#        toGet = cms.VPSet(
#            cms.PSet(
#                record = cms.string('PPSTimingCalibrationRcd'),
#                tag = cms.string('DiamondTimingCalibration')
#        )
#    )
#)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring (options.inputFiles))

process.load("CalibPPS.TimingCalibration.ppsTimingCalibrationPCLWorker_cfi")
process.DQMStore = cms.Service("DQMStore")

process.dqmOutput = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string("file:worker_output.root"),
    outputCommands = cms.untracked.vstring(
    'drop *',
    'keep *_MEtoEDMConvertPPSTimingCalib_*_*',
    ))

process.load("CalibPPS.TimingCalibration.ALCARECOPromptCalibProdPPSTimingCalib_cff")

#process.ctppsPixelDigis.inputLabel = "hltPPSCalibrationRaw"
#process.ctppsDiamondRawToDigi.rawDataTag = "hltPPSCalibrationRaw"
#process.totemRPRawToDigi.rawDataTag = "hltPPSCalibrationRaw"
#process.totemTimingRawToDigi.rawDataTag = "hltPPSCalibrationRaw"

process.path = cms.Path(
    #process.a1* 
    #process.ctppsRawToDigi *
    #process.recoCTPPS *
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
