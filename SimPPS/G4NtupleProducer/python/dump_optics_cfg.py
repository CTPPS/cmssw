"""
CMSSW config to dump optical functions from the Global Tag.

Run with:
    cmsRun SimPPS/G4NtupleProducer/python/dump_optics_cfg.py

This extracts the optics that CMSSW uses for proton reconstruction.
"""

import FWCore.ParameterSet.Config as cms

process = cms.Process("OPTICS")

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1))
process.source = cms.Source("EmptySource")

# Global Tag - this is where optics come from
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, "140X_mcRun3_2024_realistic_v20", "")

# Direct frontier connection (no proxy)
process.GlobalTag.connect = 'frontier://(serverurl=http://cmsfrontier.cern.ch:8000/FrontierProd)/CMS_CONDITIONS'

# The optics dumper
process.opticsDumper = cms.EDAnalyzer("PPSOpticsDumper")

# Path
process.p = cms.Path(process.opticsDumper)
