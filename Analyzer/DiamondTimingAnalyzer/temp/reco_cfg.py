#start with the https://github.com/cms-sw/cmssw/blob/master/RecoPPS/Local/test/diamonds_reco_cfg.py
# 1. start with AOD file 
# 2. apply SQL lite 
# 3. new AOD with shift?? -> check produce DQM files


import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

process = cms.Process('Shift')

options = VarParsing ('analysis')
options.register('sqlFileName', 
                    'sqlite_file:/afs/cern.ch/user/m/mobrzut/public/timing_sqlcalib/ppsDiamondTiming_calibration357440.sqlite',
				  VarParsing.multiplicity.singleton,
                  VarParsing.varType.string,
                  'sqlFileNameForCalib'
)

options.register('maxEventsToProcess', 
                    10000,
				  VarParsing.multiplicity.singleton,
                  VarParsing.varType.int,
                  'maxEvents to process'
)


			  				  
options.parseArguments()



# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')


#load GT:
from Configuration.AlCa.GlobalTag import GlobalTag
# 2022 prompt: to be updated
process.GlobalTag = GlobalTag(process.GlobalTag, "124X_dataRun3_v9") 

#load condDB
process.load('CondCore.CondDB.CondDB_cfi')
process.CondDB.connect =  options.sqlFileName
process.PoolDBESSource = cms.ESSource('PoolDBESSource',
    process.CondDB,
    DumpStats = cms.untracked.bool(True),
    toGet = cms.VPSet(
        cms.PSet(
            record = cms.string('PPSTimingCalibrationRcd'),
            tag = cms.string('DiamondTimingCalibration'), 
            label = cms.untracked.string('PPSTestCalibration'), # todo added by Chrys 

        )
    )
)

process.source = cms.Source ("PoolSource",
                             fileNames = cms.untracked.vstring(
# RUN 440 - "357440": [[1,354]],
"/store/data/Run2022C/AlCaPPSPrompt/ALCARECO/PPSCalMaxTracks-PromptReco-v1/000/357/440/00000/53746080-d21b-4fdf-9d19-67c6dae347ac.root",
"/store/data/Run2022C/AlCaPPSPrompt/ALCARECO/PPSCalMaxTracks-PromptReco-v1/000/357/440/00000/b8d3f012-0b59-43b4-adbc-811bcb59c9c4.root",
"/store/data/Run2022C/AlCaPPSPrompt/ALCARECO/PPSCalMaxTracks-PromptReco-v1/000/357/440/00000/cfa9f82a-5296-49cf-8353-11f8761b675c.root",
"/store/data/Run2022C/AlCaPPSPrompt/ALCARECO/PPSCalMaxTracks-PromptReco-v1/000/357/440/00000/db3c6c33-9ad0-4d61-9a64-a1aa9578550e.root"
)
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEventsToProcess)
)

process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 1000 )



# local RP reconstruction chain with standard settings
process.load("RecoPPS.Configuration.recoCTPPS_cff")

# rechits production
from Geometry.VeryForwardGeometry.geometryRPFromDB_cfi import * # TODO it looks like it can be deleted
process.load('Geometry.VeryForwardGeometry.geometryRPFromDB_cfi') #TODO: use geometry form DB not from file 

process.load('RecoPPS.Local.ctppsDiamondRecHits_cfi')

# local tracks fitter
process.load('RecoPPS.Local.ctppsDiamondLocalTracks_cfi')

process.output = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string("file:AOD.root"),
    outputCommands = cms.untracked.vstring(
        'drop *',
        'keep *_ctpps*_*_*',
    ),
)

# execution configuration
process.p = cms.Path(
    process.ctppsDiamondLocalReconstruction
)

process.outpath = cms.EndPath(process.output)