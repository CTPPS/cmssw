import FWCore.ParameterSet.Config as cms

# import default alignment settings
from CalibPPS.ESProducers.ctppsAlignment_cff import *

# import default optics settings
from CalibPPS.ESProducers.ctppsOpticalFunctions_cff import *

# import and adjust proton-reconstructions settings
from RecoPPS.ProtonReconstruction.ctppsProtons_cfi import *

print("in file ctppsProtons_cfi")
# print(ctppsProtons)

ctppsProtons.lhcInfoLabel = ctppsLHCInfoLabel
ctppsProtons.lhcInfoPerLSLabel =   ctppsLHCInfoPerLSLabel
ctppsProtons.lhcInfoPerFillLabel = ctppsLHCInfoPerFillLabel

ctppsProtons.pixelDiscardBXShiftedTracks = True
ctppsProtons.default_time = -999.
