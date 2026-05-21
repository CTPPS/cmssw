## Temporary setup instructions for 2023 simulation

These are temporary instructions to setup the direct simulation and run it

### Installation

Follow the commands below on lxplus9 to setup a working environment for the direct simulation

```bash
export SCRAM_ARCH=el9_amd64_gcc14
cmsrel CMSSW_16_1_0_pre4
cd CMSSW_16_1_0_pre4/src
cmsenv
git cms-init
git cms-merge-topic CTPPS:ab-dev-dirsim
scram b -j10
```

### Producing a config and running

With the following commands, one can easily produce and test-run configuration files for each year of data-taking. 

```bash
cd Validation/CTPPS/test/simu
./run_multiple
```

The config files created by the script are called `simu_[year]_cfg.py` and consist in a simple process that generates protons from a uniform $\xi$ distribution, runs the proton direct simulation and recontruction, and finally runs some proton-related plotters (which produce the other `simu_[year]_*.root` files). 

In particular, the code block at the bottom (example below) and the specified era modifier are the crucial items for the direct simulation to run correctly:

```python
from SimPPS.Configuration.Utils import setupPPSDirectSim
setupPPSDirectSim(process)

process.ctppsCompositeESSource.seed = 12345 # This shall be set randomly for each production
process.ctppsBeamParametersFromLHCInfoESSource.vtxOffsetX45 = 0.
process.ctppsBeamParametersFromLHCInfoESSource.vtxOffsetY45 = 0.
process.ctppsBeamParametersFromLHCInfoESSource.vtxOffsetZ45 = 0.
process.source.numberEventsInLuminosityBlock = process.ctppsCompositeESSource.generateEveryNEvents

```

The rest of the code can be adapted to the user's needs, in order to plug in the needed generator and analysers, together with any other simulation modules. Conversely, this snippet can be integrated in a user-supplied config to implement the direct simulation. 

Better setup tools, relying on `cmsDriver.py` are also present, but currently not tested with 2023 settings. 

**Please mind that this shouldn't be taken as a stable and final setup. It's an active development branch, that may change and/or break.** 
**User discretion is advised.**