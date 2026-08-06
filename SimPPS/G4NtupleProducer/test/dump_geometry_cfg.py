"""
Dump PPS geometry for arm=0, station=2, rp=3
"""

import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Run3_cff import Run3
process = cms.Process('GeomDump', Run3)

process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.source = cms.Source('EmptySource')
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1))

# Geometry and conditions
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '141X_dataRun3_Prompt_v3', '')
process.GlobalTag.connect = 'frontier://FrontierProd/CMS_CONDITIONS'

# PPS geometry
process.load('Geometry.VeryForwardGeometry.geometryRPFromDB_cfi')

# Geometry dumper
process.ppsGeometryDumper = cms.EDAnalyzer('PPSGeometryDumper',
    arm = cms.uint32(0),
    station = cms.uint32(2),
    rp = cms.uint32(3)
)

process.p = cms.Path(process.ppsGeometryDumper)
