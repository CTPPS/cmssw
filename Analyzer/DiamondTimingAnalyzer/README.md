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
