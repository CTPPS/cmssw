import FWCore.ParameterSet.Config as cms

from DQMServices.Core.DQMEDAnalyzer import DQMEDAnalyzer
totemT2DQMSource = DQMEDAnalyzer('TotemT2DQMSource',
    digisTag = cms.InputTag('totemT2Digis', 'TotemT2'),
    nbinsx = cms.uint32(25),
    nbinsy = cms.uint32(25),
    specialRunT2 = cms.bool(False),
    windowsNum = cms.uint32(8),

    perLSsaving = cms.untracked.bool(False), #driven by DQMServices/Core/python/DQMStore_cfi.py
)
