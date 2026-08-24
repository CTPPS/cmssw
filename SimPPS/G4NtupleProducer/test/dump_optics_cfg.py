"""
Dump optical functions from CMSSW conditions database.
Compares with ROOT file values to check consistency.
"""

import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Run3_cff import Run3
process = cms.Process('OpticsDump', Run3)

process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.source = cms.Source('EmptySource')
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1))

# Geometry and conditions
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '140X_mcRun3_2024_realistic_v20', '')
process.GlobalTag.connect = 'frontier://(serverurl=http://cmsfrontier.cern.ch:8000/FrontierProd)/CMS_CONDITIONS'

# PPS geometry and optics
process.load('Geometry.VeryForwardGeometry.geometryRPFromDB_cfi')

# Optics dumper
process.ppsOpticsDumper = cms.EDAnalyzer('PPSOpticsDumper')

process.p = cms.Path(process.ppsOpticsDumper)
