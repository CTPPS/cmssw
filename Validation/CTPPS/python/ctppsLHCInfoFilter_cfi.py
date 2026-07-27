from Validation.CTPPS.ctppsLHCInfoFilterDefault_cfi import ctppsLHCInfoFilterDefault as _ctppsLHCInfoFilterDefault
ctppsLHCInfoFilter = _ctppsLHCInfoFilterDefault.clone()
    
from Configuration.Eras.Modifier_run3_common_cff import run3_common
run3_common.toModify(ctppsLHCInfoFilter, useNewLHCInfo = True)

from Configuration.Eras.Modifier_ctpps_directSim_cff import ctpps_directSim
ctpps_directSim.toModify(ctppsLHCInfoFilter, useNewLHCInfo = False)
