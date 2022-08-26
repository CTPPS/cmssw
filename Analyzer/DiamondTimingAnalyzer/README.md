
## Manual tests
1. Running worker  

```bash
cmsRun python/worker.py
```

2. Running harvester
```bash 
cmsRun python/harvester.py rootFiles=file:\run_output.root 
```

Tip - logs to file - append: `2>&1 | Logs/tee harvester_log_4.txt`