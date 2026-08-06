"""
Simplified test configuration for PPSG4NtupleToRecHitProducer.
Uses inline EDProducer definition to avoid import issues.
"""

import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

# Command-line options
options = VarParsing('analysis')
options.register('inputFile',
                 'test_rechits.root',
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 'Path to G4 simulation ROOT ntuple')
options.register('verbosity', 0, VarParsing.multiplicity.singleton, VarParsing.varType.int, '')
options.parseArguments()

# Process definition
from Configuration.Eras.Era_Run3_cff import Run3
process = cms.Process('G4NtupleReco', Run3)

# Message logger
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1

# Empty source
process.source = cms.Source('EmptySource')
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(options.maxEvents))

# Geometry and conditions
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
# Use explicit GT to avoid site-local-config dependency
process.GlobalTag = GlobalTag(process.GlobalTag, '141X_dataRun3_Prompt_v3', '')
# Set frontier connection directly
process.GlobalTag.connect = 'frontier://FrontierProd/CMS_CONDITIONS'

# PPS geometry
process.load('Geometry.VeryForwardGeometry.geometryRPFromDB_cfi')

# G4 ntuple producer - defined inline to avoid import issues
process.ppsG4NtupleToRecHitProducer = cms.EDProducer('PPSG4NtupleToRecHitProducer',
    inputFileName = cms.string(options.inputFile),
    # Run3 PPS pixel detector: arm=0/1, station=2 (220m), rp=3
    arm = cms.uint32(0),
    station = cms.uint32(2),
    rp = cms.uint32(3),
    sigmaHitX = cms.double(0.029),
    sigmaHitY = cms.double(0.043),
    verbosity = cms.untracked.uint32(options.verbosity)
)

# PPS reconstruction
process.load('RecoPPS.Configuration.recoCTPPS_cff')

# Rewire pixel tracks to use our RecHits
process.ctppsPixelLocalTracks.tag = cms.InputTag('ppsG4NtupleToRecHitProducer')

# Pixel-only reconstruction
process.ctppsLocalTrackLiteProducer.includePixels = cms.bool(True)
process.ctppsLocalTrackLiteProducer.includeStrips = cms.bool(False)
process.ctppsLocalTrackLiteProducer.includeDiamonds = cms.bool(False)

# Path
process.recoPath = cms.Path(
    process.ppsG4NtupleToRecHitProducer *
    process.ctppsPixelLocalTracks
)

# Output
process.out = cms.OutputModule('PoolOutputModule',
    fileName = cms.untracked.string('g4_reco_output.root'),
    outputCommands = cms.untracked.vstring(
        'drop *',
        'keep *_ppsG4NtupleToRecHitProducer_*_*',
        'keep *_ctppsPixelLocalTracks_*_*'
    )
)
process.outpath = cms.EndPath(process.out)

print(f"Config: inputFile={options.inputFile}, maxEvents={options.maxEvents}")
