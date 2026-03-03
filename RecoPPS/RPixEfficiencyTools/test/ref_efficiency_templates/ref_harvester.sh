#!/bin/bash

set -x

CMSSW_BASE=${1}
INFILE=${2}
OFILE=${3}

# Setup CMSSW and move back to the running dir
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd $CMSSW_BASE/src
eval $(scram runtime -sh)
cd -

# Print the working directory
echo "Working directory: $(pwd)"

# Run the job
cmsRun $CMSSW_BASE/src/RecoPPS/RPixEfficiencyTools/python/ReferenceAnalysisDQMHarvester_cfg.py inputFiles=$INFILE
RETCMSSW=$?

if [ "$RETCMSSW" == "0" ]
then
    mv DQM*.root $OFILE
    RETCOPY=$?
else
    RETCOPY=1
fi

RET=$(echo "$RETCMSSW+$RETCOPY" | bc)

# if [ "$RET" == "0" ]
# then
#     ecalautomation.py $TASK jobctrl --id $JOBID --done --fields "output:${EOSDIR}/${OFILE}"
# else
#     ecalautomation.py $TASK jobctrl --id $JOBID --failed
# fi

exit $RET