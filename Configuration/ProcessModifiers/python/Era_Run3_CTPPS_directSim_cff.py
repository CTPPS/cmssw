import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Run3_cff import Run3
from Configuration.Eras.Era_Run3_2023_cff import Run3_2023
from Configuration.Eras.Modifier_ctpps_directSim_cff import ctpps_directSim

Run3_CTPPS_directSim = cms.ModifierChain(Run3,ctpps_directSim)
Run3_2023_CTPPS_directSim = cms.ModifierChain(Run3_2023,ctpps_directSim)