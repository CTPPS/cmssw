#!/bin/bash

set -x

CMSSW_BASE=${1}
JOBID=${2}
INFILE=${3}
EFFICIENCYFILE=${4}
ODIR=${5}

# Setup CMSSW and move back to the running dir
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd $CMSSW_BASE
eval $(scram runtime -sh)
cd -

# Print the working directory
echo "Working directory: $(pwd)"
OFILE=outputReferenceAnalysisDQMWorker_$JOBID.root

# Run the job
cmsRun $CMSSW_BASE/src/RecoPPS/RPixEfficiencyTools/python/ReferenceAnalysisDQMWorker_cfg.py outputFileName=$OFILE inputFiles=$INFILE efficiencyFileName=$EFFICIENCYFILE # globalTag=$GT
RETCMSSW=$?

if [ "$RETCMSSW" == "0" ]
then
    mkdir -p $ODIR
    mv $OFILE $ODIR/$OFILE
    RETCOPY=$?
else
    RETCOPY=1
fi

RET=$(echo "$RETCMSSW+$RETCOPY" | bc)

exit $RET
