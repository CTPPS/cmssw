"""
Test configuration for PPSG4NtupleToRecHitProducer.

This config reads a G4 simulation ROOT ntuple and runs the PPS reconstruction chain
starting from RecHits (bypassing digitization and clustering).

Usage:
    cmsRun template_cfg.py inputFile=path/to/g4_output.root maxEvents=100

The output ROOT file contains:
    - CTPPSPixelLocalTrack: fitted local tracks
    - CTPPSLocalTrackLite: lite version for proton reconstruction
    - reco::ForwardProton: reconstructed protons (if optics is configured)
    - edm::HepMCProduct: MC truth from the G4 simulation
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
options.register('arm',
                 0,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.int,
                 'PPS arm: 0 = sector 45, 1 = sector 56')
options.register('station',
                 0,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.int,
                 'Station: 0 = 210m, 2 = 220m')
options.register('rp',
                 3,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.int,
                 'Roman Pot index within station')
options.register('verbosity',
                 0,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.int,
                 'Verbosity level (0=quiet, 1=info, 2=debug)')
options.parseArguments()

# Process definition
from Configuration.Eras.Era_Run3_cff import Run3
process = cms.Process('G4NtupleReco', Run3)

# Message logger
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100
if options.verbosity > 0:
    process.MessageLogger.cerr.threshold = 'INFO'

# Empty source - events are driven by the G4 ntuple
process.source = cms.Source('EmptySource')
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(options.maxEvents))

# Geometry and conditions
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
# Use a Run3 global tag - adjust as needed for your specific conditions
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run3_data', '')

# PPS geometry from database
process.load('Geometry.VeryForwardGeometry.geometryRPFromDB_cfi')

# Load the G4 ntuple producer
# Use process.load() which uses cfipython correctly
process.load('SimPPS.G4NtupleProducer.ppsG4NtupleToRecHitProducer_cfi')
process.ppsG4NtupleToRecHitProducer = process.ppsG4NtupleToRecHitProducer.clone(
    inputFileName = cms.string(options.inputFile),
    arm = cms.uint32(options.arm),
    station = cms.uint32(options.station),
    rp = cms.uint32(options.rp),
    verbosity = cms.untracked.uint32(options.verbosity)
)

# Load PPS reconstruction (starts from RecHits)
process.load('RecoPPS.Configuration.recoCTPPS_cff')

# Rewire reconstruction to use G4 ntuple RecHits
# Inline the matching logic to avoid import issues
process.ctppsPixelLocalTracks.tag = cms.InputTag('ppsG4NtupleToRecHitProducer')

# For pixel-only reconstruction, we need to exclude strip and diamond processing
# Configure the lite producer to only use pixel tracks (no diamond/strip)
process.ctppsLocalTrackLiteProducer.includePixels = cms.bool(True)
process.ctppsLocalTrackLiteProducer.includeStrips = cms.bool(False)
process.ctppsLocalTrackLiteProducer.includeDiamonds = cms.bool(False)

# Only run pixel track fitting (not strip/diamond)
process.pixelOnlyRecoTask = cms.Task(
    process.ctppsPixelLocalTracks
)

# Add track dumper to write simple flat ROOT file for visualization
process.ppsTrackDumper = cms.EDAnalyzer('PPSTrackDumper',
    tracksTag = cms.InputTag('ctppsPixelLocalTracks')
)

# TFileService for track dumper output
process.TFileService = cms.Service('TFileService',
    fileName = cms.string('track_dump.root')
)

# If proton reconstruction is desired, add it (requires optics setup)
# For now, just do local track fitting
# process.pixelOnlyRecoTask.add(process.ctppsProtons)

# Define the path
process.g4NtuplePath = cms.Path(
    process.ppsG4NtupleToRecHitProducer,
    process.pixelOnlyRecoTask
)

# Separate path for track dumper (runs after reconstruction)
process.dumpPath = cms.Path(process.ppsTrackDumper)

# Output module
process.out = cms.OutputModule('PoolOutputModule',
    fileName = cms.untracked.string('g4_reco_output.root'),
    outputCommands = cms.untracked.vstring(
        'drop *',
        'keep *_ppsG4NtupleToRecHitProducer_*_*',
        'keep *_ctppsPixelLocalTracks_*_*'
    )
)

process.outpath = cms.EndPath(process.out)

# Schedule
process.schedule = cms.Schedule(process.g4NtuplePath, process.dumpPath, process.outpath)

# Print configuration summary
print(f"""
============================================
G4 Ntuple to PPS Reconstruction Configuration
============================================
Input file:  {options.inputFile}
Arm:         {options.arm} ({'sector 45' if options.arm == 0 else 'sector 56'})
Station:     {options.station} ({'210m' if options.station == 0 else '220m'})
RP:          {options.rp}
Max events:  {options.maxEvents}
Verbosity:   {options.verbosity}
============================================
""")
