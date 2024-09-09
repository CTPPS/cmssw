import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras
from Configuration.Eras.Modifier_ctpps_directSim_cff import ctpps_directSim
process = cms.Process('PPS', ctpps_directSim, eras.Run2_2018)

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('SimPPS.Configuration.directSimPPS_cff')
process.load('RecoPPS.Configuration.recoCTPPS_cff')

# minimum of logs
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.threshold = cms.untracked.string('')
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(100)

# global tag
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '131X_mcRun2_asymptotic_v3', '')

# raw data source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring("FIXME"),
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)

process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
    beamDivergenceVtxGenerator = cms.PSet(initialSeed = cms.untracked.uint32(3849)),
    ppsDirectProtonSimulation = cms.PSet(initialSeed = cms.untracked.uint32(4981))
)

from SimPPS.DirectSimProducer.matching_cff import matchDirectSimOutputsAOD
matchDirectSimOutputsAOD(process)

print('Verbosity of DirectSimProducer = '+str(process.ppsDirectProtonSimulation.verbosity))
print('Loading verbosity if still not loaded')
process.ppsDirectProtonSimulation.verbosity = cms.untracked.uint32(1)

process.p = cms.Path(
    process.directSimPPS
    * process.recoDirectSimPPS
)

# output configuration
from RecoPPS.Configuration.RecoCTPPS_EventContent_cff import RecoCTPPSAOD
RecoCTPPSAOD.outputCommands.extend(cms.untracked.vstring(
    'keep *_genPUProtons_*_*',
    'keep *_genParticles_*_*'
    )
)

process.output = cms.OutputModule('PoolOutputModule',
    fileName = cms.untracked.string('file:output.root'),
    outputCommands = RecoCTPPSAOD.outputCommands
)

process.outpath = cms.EndPath(process.output)
