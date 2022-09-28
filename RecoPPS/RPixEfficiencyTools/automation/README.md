# Automation module
Automation module contains code which is executed by Jenkins.
This module uses `automation_module` ([gitlab here]()) for submitting jobs to CRAB, HTCondor, checking the database status etc.    

## Setup 
1. Export the CRAB authentication key to the file with `--out` option. 
```
ex. voms-proxy-init -voms cms -rfc  -valid 192:00 --out ~/public/jenkins_proxy/jenkins.pem
```

## Structure 
This folder contains:
- **CrabConfigs** for submitting the EA (Efficiency Analysis) and RA (Reference Analysis) Workers to the CRAB. 
- **Engine** The engine for . You can modify this engine to create any workflow with CRAB and HTCondor tasks. More about the engine can be found in [documentation](https://indico.cern.ch/event/1075717/contributions/4523828/attachments/2312956/3951051/documentation.pdf)  

## Description of temporary setup:
- **TempSteps** contains python scripts for manual testing of single steps from EfficiencyAnalysisEngine

**Running test setup** 
```
cmsenv
source TempSetup.sh
source /cvmfs/cms.cern.ch/common/crab-setup.sh 
```
##TODO/Noticed problems
- [ ] Template parameter is named badly. It should be renamed to crab_config_file_path.   (Background: running submit_task_to_crab method shows error in CrabWrapper no template - which is confusing. )





