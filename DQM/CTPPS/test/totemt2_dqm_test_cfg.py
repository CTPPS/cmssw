import FWCore.ParameterSet.Config as cms
import string

from Configuration.Eras.Era_Run3_cff import Run3

process = cms.Process('RECODQM', Run3)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(6000) )
process.verbosity = cms.untracked.PSet( input = cms.untracked.int32(-1) )

# minimum of logs
process.MessageLogger = cms.Service("MessageLogger",
    statistics = cms.untracked.vstring(),
    destinations = cms.untracked.vstring('cerr'),
    cerr = cms.untracked.PSet(
        threshold = cms.untracked.string('WARNING')
    )
)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

# load DQM framework
process.load("DQM.Integration.config.environment_cfi")
process.dqmEnv.subSystemFolder = "CTPPS"
process.dqmEnv.eventInfoFolder = "EventInfo"
process.dqmSaver.path = ""
process.dqmSaver.tag = "CTPPS"
# raw data source
process.source = cms.Source("NewEventStreamFileReader",
    fileNames = cms.untracked.vstring(
        'file:../../../../../../CMSSW_13_1_0_pre1/t2v2p1/run364983_ls0001_streamA_StorageManager.dat',
        'file:../../../../../../CMSSW_13_1_0_pre1/t2v2p1/run364983_ls0002_streamA_StorageManager.dat',
#        '/store/data/Run2018D/ZeroBias/RAW/v1/000/324/747/00000/97A72F4B-786F-5A48-B97E-C596DD73BD77.root',
    )
)

#process.source = cms.Source('PoolSource',
#    fileNames = cms.untracked.vstring(
#        '/store/data/Run2018D/ZeroBias/RAW/v1/000/324/747/00000/97A72F4B-786F-5A48-B97E-C596DD73BD77.root',
#    ),
#)

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '130X_dataRun3_HLT_v2', '')

#Raw-to-digi
process.load('EventFilter.CTPPSRawToDigi.ctppsRawToDigi_cff')

#process.load('CalibPPS.ESProducers.totemT2DAQMapping_cff')

process.load('Geometry.ForwardCommonData.totemT22021V2XML_cfi')
process.load('Geometry.ForwardGeometry.totemGeometryESModule_cfi')
process.load('RecoPPS.Local.totemT2RecHits_cfi')
process.load('DQM.CTPPS.totemT2DQMSource_cfi')
process.totemDAQMappingESSourceXML_TotemT2.verbosity = cms.untracked.uint32(0)
process.totemT2Digis.RawUnpacking.verbosity = cms.untracked.uint32(0)
process.totemT2Digis.RawToDigi.verbosity = cms.untracked.uint32(0)
process.totemT2Digis.RawToDigi.printUnknownFrameSummary = cms.untracked.uint32(0)
process.totemT2Digis.RawToDigi.printErrorSummary = cms.untracked.uint32(0)
process.totemDAQMappingESSourceXML_TotemT2.multipleChannelsPerPayload = cms.untracked.bool(True)

process.path = cms.Path(
    process.ctppsRawToDigi *
    process.totemT2Digis *
    process.totemT2RecHits *
    process.totemT2DQMSource
)

process.end_path = cms.EndPath(
    process.dqmEnv +
    process.dqmSaver
)

process.schedule = cms.Schedule(
    process.path,
    process.end_path
)
