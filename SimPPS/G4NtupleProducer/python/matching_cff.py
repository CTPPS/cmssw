import FWCore.ParameterSet.Config as cms

def matchG4NtupleOutputs(process):
    """
    Rewire the PPS reconstruction chain to use RecHits from the G4 ntuple producer.
    This is analogous to matchDirectSimOutputs in SimPPS/DirectSimProducer/python/matching_cff.py
    """
    # Match sources of rechits with G4 ntuple producer outputs
    # Note: Only pixel tracks are produced from the G4 ntuple (no strips, no diamonds)
    process.ctppsPixelLocalTracks.tag = cms.InputTag('ppsG4NtupleToRecHitProducer')

    # Handle clashes between simulation and GT conditions
    # These ES preferences ensure the right conditions are used
    if hasattr(process, 'ctppsCompositeESSource'):
        process.es_prefer_composrc = cms.ESPrefer('CTPPSCompositeESSource', 'ctppsCompositeESSource')
    if hasattr(process, 'ppsPixelTopologyESSource'):
        process.es_prefer_pixtopo = cms.ESPrefer('PPSPixelTopologyESSource', 'ppsPixelTopologyESSource')
    if hasattr(process, 'ctppsBeamParametersFromLHCInfoESSource'):
        process.es_prefer_lhcinfo = cms.ESPrefer('CTPPSBeamParametersFromLHCInfoESSource', 'ctppsBeamParametersFromLHCInfoESSource')
    if hasattr(process, 'ppsAssociationCutsESSource'):
        process.es_prefer_assocuts = cms.ESPrefer('PPSAssociationCutsESSource', 'ppsAssociationCutsESSource')

    return process


def matchG4NtupleOutputsMinimal(process):
    """
    Minimal rewiring for pixel-only reconstruction from G4 ntuple.
    Only sets the pixel local tracks input tag.
    """
    process.ctppsPixelLocalTracks.tag = cms.InputTag('ppsG4NtupleToRecHitProducer')
    return process
