## Project setup 
Instruction on how to setup and run efficiency suite in CMSSW_12_4_0 environment:

1. Preprare CMSSW_12_4_0 environment with: `cmsrel CMSSW_12_4_0`
2. Go to newly created directory: `cd CMSSW_12_4_0`
3. Setup cms environmental variables(repeat every time you log in): `cmsenv`
4. Merge with the main branch containing efficiency suite: `git cms-merge-checkout CTPPS:MO/timing/migration_to_CMSSW_12_4_6_clean_branch` #(TODO change it to the main )
5. Compile the solution `scram b -j10`
6. Switch directory to the one containing the source code: `cd src/Analyzer/DiamondTimingAnalyzer`

## Run rolling calibrations: 


## More Details:
### Applying time shift 
### Finding Test File


## Project Structure
### Config.json
Config.json describes the structure of detectors. 
- Sector - left of right arm
- station - the pot in which the detector are located
- plane - the detector contains planes - 
- channel - each plane is divided into channels. They are counted in U-shape, due to the phisical connection structure. 
## Manual tests
1. Running worker  

```bash
cmsRun python/worker.py
```

2. Running harvester
```bash 
cmsRun python/harvester.py rootFiles=file:run_output.root 
```

Tip - logs to file - append: 
```
 2>&1 | tee Logs/worker_log_5.txt
```

3. Submit to CRAB:
```cmsemv```
```voms-proxy-init -voms cms -rfc  -valid 192:00```
```crab submit submitTimingToCRAB.py ```
