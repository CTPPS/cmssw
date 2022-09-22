#config
runNumber=357440
maxEvents=90000
constPREFIX=JUST_SHIFT
experimentName=_diamondDQM_chrys
try=try4
resultsFolderName=${constPREFIX}_${runNumber}_maxEv${maxEvents}_${experimentName}
#local variables
# workerPrefix=${resultsFolderName}/worker_${constPREFIX}_${runNumber}_maxEv_${maxEvents}
sqlFileName=sqlite_file:/afs/cern.ch/user/m/mobrzut/public/timing_sqlcalib/ppsDiamondTiming_calibration${runNumber}.sqlite
# scritpt
mkdir $resultsFolderName
mkdir $resultsFolderName/$try
echo "######################################## ##########################################"
echo "shift 0:"
echo "ls *.root" # make sure that the directory doesn't contain any root files
shift="cmsRun temp/produceAOD_shift_fromRoot_justShift/reco_cfg_4_diamondDQM_test.py " 
# shift="cmsRun python/reco_cfg_AOD_DQM.py " 

# shift+="outputFileName=${workerPrefix}_loop_0.root "
# shift+="calibInput=calib_357440_90k.json "
# shift+="sqlFileName=${sqlFileName} "
# shift+="maxEventsToProcess=${maxEvents}"
echo $shift

echo "mv *.root ${resultsFolderName}/${try}"


echo "--------------------------------------------------------------------------------------------------"
