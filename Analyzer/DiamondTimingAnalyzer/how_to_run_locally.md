HOW TO RUN THE CODE LOCALLY\
-cmsrel CMSSW_13_0_9\
-cd CMSSW_13_0_9/src\
-cmsenv\
-git cms-merge-topic Milosz-O:rolling_calib_hptdc_timing_13_0_9\
-cd Analyzer/DiamondTimingAnalyzer\
-scram b -j8

#1st iteration:\
-cmsRun python/worker_withReRecoShift.py\
#By default number of events for the worker is set to all the events, we can specify the number of events we want with, for example maxEventsToProcess=10000\
-cmsRun python/harvester.py

#Next iterations:\
#Rename the calibration file and the output file produced by the harvester after each iteration so they don't get overwritten in the next one\
-cmsRun python/worker_withReRecoShift.py calibInput={previous calibration file name}.json\
-cmsRun python/harvester.py loopIndex={1 if it's the second iteration, then 2,3 etc.} calibFiles={name of the calibration file nr1}.json,{name of the calibration file nr2}.json,...

