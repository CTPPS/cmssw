#config
runNumber=357440
maxEvents=450000
experimentName=WORKER_ON_SHIFTED_FILE_by_chrysFile_7
resultsFolderName=test_run_${runNumber}_maxEv${maxEvents}_${experimentName}
# sqlFileName=sqlite_file:/afs/cern.ch/user/m/mobrzut/public/timing_sqlcalib/ppsDiamondTiming_calibration${runNumber}.sqlite

# worker setup
workerFile=temp/ExampleConfig_cfg_1_2_without_shift_multisource.py
shiftExperimentPrefix=JUST_SHIFT_357440_maxEv90000
rootInput=file:${PWD}/${shiftExperimentPrefix}__diamondDQM_chrys/try3/AOD.root
workerPrefix=${resultsFolderName}/worker_timeshift_${runNumber}_maxEv_${maxEvents}

# harvester setup

# scritpt
mkdir $resultsFolderName
echo "######################################## ##########################################"
echo "WORKER 0:"
worker0="cmsRun ${workerFile} " 
worker0+="outputFileName=${workerPrefix}_loop_0.root "
worker0+="rootInput=${rootInput} "
worker0+="useDB=True "
# worker0+="calibInput=${PWD}/calib_${runNumber}.json "
# worker0+="sqlFileName=${sqlFileName} "
worker0+="maxEventsToProcess=${maxEvents} "
worker0+="2>&1 | tee ${resultsFolderName}/worker_log_0.txt "
echo $worker0
echo "HARVESTER 0:"
harvester0="cmsRun temp/harvester2_multisource.py "
harvester0+="rootFiles=file:${workerPrefix}_loop_0.root "
harvester0+="calibOutput=${resultsFolderName}/calib_loop_0.json "
harvester0+="useDB=True "
# harvester0+="calibInput=${resultsFolderName}/calib_${runNumber}.json "
# harvester0+="sqlFileName=${sqlFileName} "
harvester0+="outputDirectoryRoot=${resultsFolderName} "
harvester0+="2>&1 | tee ${resultsFolderName}/harvester_log_0.txt "

echo $harvester0


calibFilesList="${resultsFolderName}/calib_loop_0.json"
for i in 1 2
do
    prev_loop=$(($i-1))
    if [ $prev_loop -gt 0 ]  
    then        
        calibFilesList="${calibFilesList},${resultsFolderName}/calib_loop_${prev_loop}.json"
    fi
    echo "**************************************************************************************"
    echo "WORKER $i"
    workeri="cmsRun ${workerFile} "
    workeri+="rootInput=${rootInput} "
    # workeri+="rootInput=file:${workerPrefix}_loop_0.root" # it looks like it should be always the same 
    workeri+="outputFileName=${workerPrefix}_loop_$i.root "
    # workeri+="useDB=True "
    workeri+="calibInput=${resultsFolderName}/calib_loop_$prev_loop.json "
    # workeri+="sqlFileName=${sqlFileName} "
    workeri+="maxEventsToProcess=${maxEvents} " 
    workeri+="2>&1 | tee ${resultsFolderName}/worker_log_$i.txt "
    echo $workeri

    echo "HARVESTER $i"
    harvesteri="cmsRun temp/harvester2_multisource.py "
    harvesteri+="rootFiles=file:${workerPrefix}_loop_$i.root "
    harvesteri+="calibOutput=${resultsFolderName}/calib_loop_$i.json " 
    harvesteri+="loopIndex=$i "
    harvesteri+="calibFiles=$calibFilesList " 
    # harvesteri+="useDB=True "
    harvesteri+="calibInput=${resultsFolderName}/calib_loop_$prev_loop.json "
    # harvesteri+="sqlFileName=${sqlFileName} "
    harvesteri+="outputDirectoryRoot=${resultsFolderName} "
    harvesteri+="2>&1 | tee ${resultsFolderName}/harvester_log_$i.txt "
    echo $harvesteri
done
echo "--------------------------------------------------------------------------------------------------"
