"""
PPS G4 Ntuple to Multi-RP Reconstruction configuration.
Reads G4 dual-station simulation ROOT ntuple (12 planes) and runs PPS
multi-RP proton reconstruction using both 210m and 220m stations.

Usage:
    cmsRun run_reco_multirp_cfg.py inputFile=g4_dual.root outputFile=reco.root maxEvents=-1
"""

import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

# Command-line options
options = VarParsing('analysis')
options.register('inputFile', 'g4_output.root', VarParsing.multiplicity.singleton,
                 VarParsing.varType.string, 'Input G4 ntuple (dual-station)')
options.register('recoFile', 'reco_output.root', VarParsing.multiplicity.singleton,
                 VarParsing.varType.string, 'Output EDM file')
options.register('trackDumpFile', 'track_dump.root', VarParsing.multiplicity.singleton,
                 VarParsing.varType.string, 'Output flat track dump')
options.register('arm', 0, VarParsing.multiplicity.singleton, VarParsing.varType.int,
                 'PPS arm: 0=sector45, 1=sector56')
options.register('rp', 3, VarParsing.multiplicity.singleton, VarParsing.varType.int,
                 'Roman Pot index')
options.register('verbosity', 0, VarParsing.multiplicity.singleton, VarParsing.varType.int, '')
options.register('nThreads', 1, VarParsing.multiplicity.singleton, VarParsing.varType.int,
                 'Number of threads for multi-threaded processing')
options.parseArguments()

# Process
from Configuration.Eras.Era_Run3_cff import Run3
process = cms.Process('PPSRecoMultiRP', Run3)

# Multi-threading
if options.nThreads > 1:
    process.options.numberOfThreads = cms.untracked.uint32(options.nThreads)
    process.options.numberOfStreams = cms.untracked.uint32(options.nThreads)

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
process.GlobalTag = GlobalTag(process.GlobalTag, '140X_mcRun3_2024_realistic_v20', '')
# Direct frontier connection (no proxy) - for running outside CERN
process.GlobalTag.connect = 'frontier://(serverurl=http://cmsfrontier.cern.ch:8000/FrontierProd)/CMS_CONDITIONS'

# PPS geometry
process.load('Geometry.VeryForwardGeometry.geometryRPFromDB_cfi')

# G4 ntuple producer - DUAL STATION MODE
# Reads planes 1-6 (station 0, 210m) and planes 7-12 (station 2, 220m)
process.ppsG4NtupleToRecHitProducer = cms.EDProducer('PPSG4NtupleToRecHitProducer',
    inputFileName = cms.string(options.inputFile),
    arm = cms.uint32(options.arm),
    station = cms.uint32(0),  # Ignored in dual station mode
    rp = cms.uint32(options.rp),
    dualStation = cms.bool(True),  # Enable dual-station (multi-RP) mode
    sigmaHitX = cms.double(0.029),
    sigmaHitY = cms.double(0.043),
    verbosity = cms.untracked.uint32(options.verbosity)
)

# PPS reconstruction
process.load('RecoPPS.Configuration.recoCTPPS_cff')
process.ctppsPixelLocalTracks.tag = cms.InputTag('ppsG4NtupleToRecHitProducer')

# Pixel-only config for lite tracks
process.ctppsLocalTrackLiteProducer.includePixels = cms.bool(True)
process.ctppsLocalTrackLiteProducer.includeStrips = cms.bool(False)
process.ctppsLocalTrackLiteProducer.includeDiamonds = cms.bool(False)

# Proton reconstruction - MULTI-RP MODE
# Both single-RP and multi-RP are enabled for comparison
process.ctppsProtons.doSingleRPReconstruction = cms.bool(True)
process.ctppsProtons.doMultiRPReconstruction = cms.bool(True)

# Custom association cuts for G4 simulation
# We use very relaxed cuts since our G4 geometry doesn't perfectly match CMSSW
from CalibPPS.ESProducers.ppsAssociationCutsESSource_cfi import ppsAssociationCutsESSource as _ppsAssociationCutsESSource

process.ppsAssociationCutsESSource = _ppsAssociationCutsESSource.clone(
    configuration = cms.VPSet(
        cms.PSet(
            validityRange = cms.EventRange("0:min - 999999:max"),
            # Sector 45 (arm 0) - very relaxed cuts for G4 simulation
            association_cuts_45 = cms.PSet(
                # Position cuts: allow up to 10 mm difference
                x_cut_mean = cms.string("0.0"),
                x_cut_threshold = cms.string("10.0"),
                y_cut_mean = cms.string("0.0"),
                y_cut_threshold = cms.string("10.0"),
                # Xi cut: allow up to 0.1 difference
                xi_cut_mean = cms.string("0.0"),
                xi_cut_threshold = cms.string("0.1"),
                # Theta_y cut: allow up to 100 mrad difference (very relaxed for G4 sim)
                th_y_cut_mean = cms.string("0.0"),
                th_y_cut_threshold = cms.string("0.1"),
                # Timing-tracking association (not used for pixel-only)
                ti_tr_min = cms.double(-100.0),
                ti_tr_max = cms.double(100.0),
            ),
            # Sector 56 (arm 1) - same relaxed cuts
            association_cuts_56 = cms.PSet(
                x_cut_mean = cms.string("0.0"),
                x_cut_threshold = cms.string("10.0"),
                y_cut_mean = cms.string("0.0"),
                y_cut_threshold = cms.string("10.0"),
                xi_cut_mean = cms.string("0.0"),
                xi_cut_threshold = cms.string("0.1"),
                th_y_cut_mean = cms.string("0.0"),
                th_y_cut_threshold = cms.string("0.1"),
                ti_tr_min = cms.double(-100.0),
                ti_tr_max = cms.double(100.0),
            ),
        )
    )
)

# Tell ctppsProtons to use our custom association cuts
process.ctppsProtons.ppsAssociationCutsLabel = cms.string("")

# Prefer our custom association cuts over the GlobalTag ones
process.es_prefer_ppsAssociationCuts = cms.ESPrefer("PPSAssociationCutsESSource", "ppsAssociationCutsESSource")

# Track dumper for flat ROOT output (local track level)
process.ppsTrackDumper = cms.EDAnalyzer('PPSTrackDumper',
    tracksTag = cms.InputTag('ctppsPixelLocalTracks')
)

# Proton dumper for full kinematics (xi, theta_x, theta_y, t, etc.)
# Both single-RP and multi-RP protons are dumped
process.ppsProtonDumper = cms.EDAnalyzer('PPSProtonDumper',
    useSingleRP = cms.bool(True),
    useMultiRP = cms.bool(True),
    singleRPTag = cms.InputTag('ctppsProtons', 'singleRP'),
    multiRPTag = cms.InputTag('ctppsProtons', 'multiRP')
)

# TFileService for track dumper
process.TFileService = cms.Service('TFileService',
    fileName = cms.string(options.trackDumpFile)
)

# Path: RecHits -> Local Tracks -> Track Lite -> Protons
process.recoPath = cms.Path(
    process.ppsG4NtupleToRecHitProducer *
    process.ctppsPixelLocalTracks *
    process.ctppsLocalTrackLiteProducer *
    process.ctppsProtons *
    process.ppsTrackDumper *
    process.ppsProtonDumper
)

# Output
process.out = cms.OutputModule('PoolOutputModule',
    fileName = cms.untracked.string(options.recoFile),
    outputCommands = cms.untracked.vstring(
        'drop *',
        'keep *_ppsG4NtupleToRecHitProducer_*_*',
        'keep *_ctppsPixelLocalTracks_*_*',
        'keep *_ctppsLocalTrackLiteProducer_*_*',
        'keep *_ctppsProtons_*_*'  # Both singleRP and multiRP collections
    )
)

process.outpath = cms.EndPath(process.out)

# Print config
print(f"=== Multi-RP Reconstruction Configuration ===")
print(f"Input:  {options.inputFile}")
print(f"Output: {options.recoFile}")
print(f"Tracks: {options.trackDumpFile}")
print(f"Events: {options.maxEvents}")
print(f"Threads: {options.nThreads}")
print(f"Detector: arm={options.arm}, rp={options.rp}")
print(f"Stations: 0 (210m near) + 2 (220m far) - dual station mode")
print(f"Reconstruction: single-RP + multi-RP")
print(f"=============================================")
