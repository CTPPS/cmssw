#jc_lhcInfo_split

##Set CMSSW workspace and compile

```bash
source /cvmfs/cms.cern.ch/cmsset_default.sh
scram project CMSSW CMSSW_12_5_0_pre3
cd CMSSW_12_5_0_pre3/src/
cmsenv
git cms-init
git cms-checkout-topic CTPPS:jc_lhcInfo_split
scram b -j 8
```

##Running tests
assuming being in the CMSSW_12_5_0_pre3/src directory and having run cmsenv 

```bash
cd CondTools/RunInfo/test
cmsRun LHCInfoPerLSWriter_cfg.py 
cmsRun LHCInfoPerFillWriter_cfg.py 
cmsRun LHCInfoPerLSAnalyzer_cfg.py 
cmsRun LHCInfoPerFillAnalyzer_cfg.py 
```