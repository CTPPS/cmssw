import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
from RecoPPS.Local.ctppsDiamondLocalReconstruction_cff import *
from Geometry.VeryForwardGeometry.geometryRPFromDB_cfi import *

process = cms.Process("TIMINGSTUDY")
options = VarParsing ('analysis')
saveToDQM = False

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )
process.verbosity = cms.untracked.PSet( input = cms.untracked.int32(-1) )

# minimum of logs
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

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
				  'planes.json',
				  VarParsing.multiplicity.singleton,
				  VarParsing.varType.string,
				  "planes config file")

options.register ('validOOT',
				  -1,
				  VarParsing.multiplicity.singleton,
				  VarParsing.varType.int,
				  "valid OOT slice")				  
				  
options.parseArguments()
				  

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
# process.load('Geometry.VeryForwardGeometry.geometryRPFromDB_cfi') #TODO: use geometry form DB not from file 
process.load('Geometry.VeryForwardGeometry.geometryRPFromDD_2022_cfi')


process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    # "file:/eos/cms/store/data/Run2022B/AlCaPPS/RAW/v1/000/355/207/00000/c23440f4-49c0-44aa-b8f6-f40598fb4705.root"
    "file:/eos/home-m/mobrzut/357/734/00000/73c5eeec-578d-438f-92f1-3e1fff2ed390.root"
))

from DQMServices.Core.DQMEDAnalyzer import DQMEDAnalyzer
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '124X_dataRun3_v9', '')

if options.calibInput == '':
    process.GlobalTag.toGet = cms.VPSet(
        cms.PSet(record = cms.string('PPSTimingCalibrationRcd'),
                    tag = cms.string('PPSDiamondTimingCalibration_Run3_v1_hlt'), # working tag: PPSDiamondTimingCalibration_Run3_v1_hlt
                    #TODO: old tag PPSDiamondTimingCalibration_v1  - to delete
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
process.load('RecoPPS.Local.ctppsDiamondLocalReconstruction_cff')

process.diamondTimingWorker = DQMEDAnalyzer("DiamondTimingWorker",
    tagDigi = cms.InputTag("ctppsDiamondRawToDigiAlCaRecoProducer", "TimingDiamond"),
    tagRecHit = cms.InputTag("ctppsDiamondRecHitsAlCaRecoProducer"),
    tagPixelLocalTrack = cms.InputTag("ctppsPixelLocalTracksAlCaRecoProducer"),
    tagLocalTrack = cms.InputTag("ctppsDiamondLocalTracksAlCaRecoProducer"),
    tagValidOOT = cms.int32(-1),
    planesConfig = cms.string(options.planesConfig),
    Ntracks_Lcuts = cms.vint32([-1,1,-1,1]),
    Ntracks_Ucuts = cms.vint32([-1,6,-1,6]),
)

#process.TFileService = cms.Service("TFileService",
#    fileName = cms.string(options.outputFile)
#)
if(saveToDQM):
    process.load("DQMServices.Core.DQMStore_cfi")
    process.load("DQMServices.Components.DQMEnvironment_cfi")
    process.dqmEnv.subSystemFolder = "CalibPPS"
    process.dqmSaver.convention = 'Offline'
    process.dqmSaver.workflow = "/CalibPPS/TimingCalibration/CMSSW_12_0_0_pre2"
    process.dqmSaver.saveByRun = -1
    process.dqmSaver.saveAtJobEnd = True
    process.dqmSaver.forceRunNumber = 999999


    process.dqmOutput = cms.OutputModule("DQMRootOutputModule",
        fileName = cms.untracked.string("harvester_output.root")
    )

    process.load("DQMServices.Components.EDMtoMEConverter_cff")

else: 
    process.dqmOutput = cms.OutputModule("DQMRootOutputModule", fileName=cms.untracked.string(options.outputFileName))

process.DQMStore = cms.Service("DQMStore")

#print(vars(process))

process.path = cms.Path(
   process.diamondTimingWorker
)

if(saveToDQM):
    process.end_path = cms.EndPath(
    process.dqmEnv +
    process.dqmSaver)
else:
    process.end_path = cms.EndPath(process.dqmOutput)

process.schedule = cms.Schedule(
    process.path,
    process.end_path
)
