import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras
process = cms.Process('PPS', eras.Run3)

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
process.GlobalTag = GlobalTag(process.GlobalTag, '140X_mcRun3_2024_realistic_v20', '')

# raw data source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #'/store/mc/RunIISummer20UL16MiniAODAPVv2/GGToMuMu_Pt-25_Elastic_13TeV-lpair/MINIAODSIM/106X_mcRun2_asymptotic_preVFP_v11-v1/280000/3870E880-4A47-7440-B122-C76062D2290F.root',
#        '/store/mc/RunIISummer20UL18MiniAODv2/GGToMuMu_Pt-25_Elastic_13TeV-lpair/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/260000/6EC2EE65-77C0-3C43-A51B-7B94FE441894.root',
        '/store/mc/Run3Summer22DRPremix/GGToMuMu_PT-25_El-El_13p6TeV_lpair/AODSIM/124X_mcRun3_2022_realistic_v12-v2/80000/4b641be6-f3ca-4a8c-9425-c73c7f97675f.root'
    ),
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
    beamDivergenceVtxGenerator = cms.PSet(initialSeed = cms.untracked.uint32(3849)),
    ppsDirectProtonSimulation = cms.PSet(initialSeed = cms.untracked.uint32(4981))
)

from SimPPS.DirectSimProducer.matching_cff import matchDirectSimOutputsAOD
matchDirectSimOutputsAOD(process)

process.ctppsBeamParametersFromLHCInfoESSource.vtxOffsetX45 = -0.1078 * 1E-1
process.ctppsBeamParametersFromLHCInfoESSource.vtxOffsetY45 = -0.4189 * 1E-1
process.ctppsBeamParametersFromLHCInfoESSource.vtxOffsetZ45 = -0.2488 * 1E-1
process.ctppsBeamParametersFromLHCInfoESSource.vtxOffsetX56 = -0.1078 * 1E-1
process.ctppsBeamParametersFromLHCInfoESSource.vtxOffsetY56 = -0.4189 * 1E-1
process.ctppsBeamParametersFromLHCInfoESSource.vtxOffsetZ56 = -0.2488 * 1E-1

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
