"""
PPS G4 Ntuple to Reconstruction configuration.
Reads G4 simulation ROOT ntuple and runs PPS local track reconstruction.

Usage:
    cmsRun run_reco_cfg.py inputFile=g4.root outputFile=reco.root maxEvents=-1
"""

import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

# Command-line options
options = VarParsing('analysis')
options.register('inputFile', 'g4_output.root', VarParsing.multiplicity.singleton,
                 VarParsing.varType.string, 'Input G4 ntuple')
options.register('recoFile', 'reco_output.root', VarParsing.multiplicity.singleton,
                 VarParsing.varType.string, 'Output EDM file')
options.register('trackDumpFile', 'track_dump.root', VarParsing.multiplicity.singleton,
                 VarParsing.varType.string, 'Output flat track dump')
options.register('arm', 0, VarParsing.multiplicity.singleton, VarParsing.varType.int,
                 'PPS arm: 0=sector45, 1=sector56')
options.register('station', 2, VarParsing.multiplicity.singleton, VarParsing.varType.int,
                 'Station: 0=210m, 2=220m')
options.register('rp', 3, VarParsing.multiplicity.singleton, VarParsing.varType.int,
                 'Roman Pot index')
options.register('verbosity', 0, VarParsing.multiplicity.singleton, VarParsing.varType.int, '')
options.parseArguments()

# Process
from Configuration.Eras.Era_Run3_cff import Run3
process = cms.Process('PPSReco', Run3)

# Message logger
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100 if options.maxEvents > 100 else 1

# Source
process.source = cms.Source('EmptySource')
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(options.maxEvents))

# Geometry and conditions
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '141X_dataRun3_Prompt_v3', '')
process.GlobalTag.connect = 'frontier://FrontierProd/CMS_CONDITIONS'

# PPS geometry
process.load('Geometry.VeryForwardGeometry.geometryRPFromDB_cfi')

# G4 ntuple producer
process.ppsG4NtupleToRecHitProducer = cms.EDProducer('PPSG4NtupleToRecHitProducer',
    inputFileName = cms.string(options.inputFile),
    arm = cms.uint32(options.arm),
    station = cms.uint32(options.station),
    rp = cms.uint32(options.rp),
    sigmaHitX = cms.double(0.029),
    sigmaHitY = cms.double(0.043),
    verbosity = cms.untracked.uint32(options.verbosity)
)

# PPS reconstruction
process.load('RecoPPS.Configuration.recoCTPPS_cff')
process.ctppsPixelLocalTracks.tag = cms.InputTag('ppsG4NtupleToRecHitProducer')

# Pixel-only config
process.ctppsLocalTrackLiteProducer.includePixels = cms.bool(True)
process.ctppsLocalTrackLiteProducer.includeStrips = cms.bool(False)
process.ctppsLocalTrackLiteProducer.includeDiamonds = cms.bool(False)

# Track dumper for flat ROOT output
process.ppsTrackDumper = cms.EDAnalyzer('PPSTrackDumper',
    tracksTag = cms.InputTag('ctppsPixelLocalTracks')
)

# TFileService for track dumper
process.TFileService = cms.Service('TFileService',
    fileName = cms.string(options.trackDumpFile)
)

# Path
process.recoPath = cms.Path(
    process.ppsG4NtupleToRecHitProducer *
    process.ctppsPixelLocalTracks *
    process.ppsTrackDumper
)

# Output
process.out = cms.OutputModule('PoolOutputModule',
    fileName = cms.untracked.string(options.recoFile),
    outputCommands = cms.untracked.vstring(
        'drop *',
        'keep *_ppsG4NtupleToRecHitProducer_*_*',
        'keep *_ctppsPixelLocalTracks_*_*'
    )
)
# todo kinematics of the proton

process.outpath = cms.EndPath(process.out)

# Print config
print(f"Input:  {options.inputFile}")
print(f"Output: {options.recoFile}")
print(f"Tracks: {options.trackDumpFile}")
print(f"Events: {options.maxEvents}")
print(f"Detector: arm={options.arm} station={options.station} rp={options.rp}")
