import FWCore.ParameterSet.Config as cms

ppsG4NtupleToRecHitProducer = cms.EDProducer('PPSG4NtupleToRecHitProducer',
    # Path to the G4 simulation ROOT ntuple
    inputFileName = cms.string('g4_output.root'),

    # Detector identification
    # arm: 0 = sector 45 (negative z), 1 = sector 56 (positive z)
    # station: 0 = 210m station, 2 = 220m station
    # rp: Roman Pot index within station (typically 3 for pixel RP)
    arm = cms.uint32(0),
    station = cms.uint32(0),
    rp = cms.uint32(3),

    # Hit resolution (sigma in mm)
    # Default values: pixel pitch / sqrt(12)
    # Pixel pitch: 100 um (X) x 150 um (Y)
    sigmaHitX = cms.double(0.029),  # 100 um / sqrt(12) in mm
    sigmaHitY = cms.double(0.043),  # 150 um / sqrt(12) in mm

    verbosity = cms.untracked.uint32(0)
)
