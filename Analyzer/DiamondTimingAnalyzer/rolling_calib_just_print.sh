
#config
runNumber=357440
maxEvents=90000
experimentName=loop0shift_workerUseCalib
resultsFolderName=test_run_${runNumber}_maxEv${maxEvents}_${experimentName}
#local variables
workerPrefix=${resultsFolderName}/worker_timeshift_${runNumber}_maxEv_${maxEvents}
sqlFileName=sqlite_file:/afs/cern.ch/user/m/mobrzut/public/timing_sqlcalib/ppsDiamondTiming_calibration${runNumber}.sqlite
# scritpt
mkdir $resultsFolderName
echo "######################################## ##########################################"
echo "WORKER 0:"
worker0="cmsRun temp/ExampleConfig_cfg_2_2_with_shift_multisource.py " 
worker0+="outputFileName=${workerPrefix}_loop_0.root "
worker0+="calibInput=${PWD}/calib_${runNumber}.json "
# worker0+="sqlFileName=${sqlFileName} "
worker0+="maxEventsToProcess=${maxEvents}"
echo $worker0
echo "HARVESTER 0:"
harvester0="cmsRun temp/harvester2_multisource.py "
harvester0+="rootFiles=file:${workerPrefix}_loop_0.root "
harvester0+="calibOutput=${resultsFolderName}/calib_loop_0.json "
harvester0+="calibInput=${PWD}/calib_${runNumber}.json "
# harvester0+="sqlFileName=${sqlFileName} "
harvester0+="outputDirectoryRoot=${resultsFolderName}"
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
    workeri="cmsRun temp/ExampleConfig_cfg.py "
    workeri+="outputFileName=${workerPrefix}_loop_$i.root "
    # workeri+="calibInput=calib_loop_$prev_loop.json "
    workeri+="useSQLFile=True "
    workeri+="sqlFileName=${sqlFileName} "
    workeri+="maxEventsToProcess=${maxEvents}" 
    echo $workeri

    echo "HARVESTET $i"
    harvesteri="cmsRun temp/harvester2.py "
    harvesteri+="rootFiles=file:${workerPrefix}_loop_$i.root "
    # harvesteri+="calibInput=calib_loop_$prev_loop.json "
    harvesteri+="calibOutput=${resultsFolderName}/calib_loop_$i.json " 
    harvesteri+="loopIndex=$i "
    harvesteri+="calibFiles=$calibFilesList " 
    harvesteri+="sqlFileName=${sqlFileName} "
    harvesteri+="outputDirectoryRoot=${resultsFolderName} "
    echo $harvesteri
done
echo "--------------------------------------------------------------------------------------------------"


