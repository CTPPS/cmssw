#!/bin/bash
function die { echo $1: status $2; exit $2; }

customise_commands="process.GlobalTag.toGet = cms.VPSet()\
\nprocess.GlobalTag.toGet.append(cms.PSet(record = cms.string(\"CTPPSOpticsRcd\"),tag =  cms.string(\"PPSOpticalFunctions_test\"), connect = cms.string(\"sqlite_file:/afs/cern.ch/user/w/wcarvalh/public/CTPPS/optical_functions/validation/Aug22/PPSOpticalFunctions_singleIOV_2022_with_verticals.db\")))\
\nprocess.GlobalTag.toGet.append(cms.PSet(record = cms.string(\"LHCInfoRcd\"),tag =  cms.string(\"LHCInfo_end_2\"), connect = cms.string(\"sqlite_file:/afs/cern.ch/user/j/jchyczyn/public/reco_test_dbs/lhcinfo_8083.db\")))\
\nprocess.GlobalTag.toGet.append(cms.PSet(record = cms.string(\"LHCInfoPerLSRcd\"),tag =  cms.string(\"LHCInfoPerLS_end_2\"), connect = cms.string(\"sqlite_file:/afs/cern.ch/user/j/jchyczyn/public/reco_test_dbs/lhcinfo_8083.db\")))\
\nprocess.GlobalTag.toGet.append(cms.PSet(record = cms.string(\"LHCInfoPerFillRcd\"),tag =  cms.string(\"LHCInfoPerFill_end_2\"), connect = cms.string(\"sqlite_file:/afs/cern.ch/user/j/jchyczyn/public/reco_test_dbs/lhcinfo_8083.db\")))\
\nprocess.GlobalTag.toGet.append(cms.PSet(record = cms.string(\"AlCaRecoTriggerBitsRcd\"),tag =  cms.string(\"AlCaRecoHLTpaths_PPS2022_express_v1\"), connect = cms.string(\"frontier://FrontierProd/CMS_CONDITIONS\")))\
\nprocess.GlobalTag.toGet.append(cms.PSet(record = cms.string(\"PPSTimingCalibrationLUTRcd\"),tag =  cms.string(\"PPSDiamondTimingCalibrationLUT_test\"), connect = cms.string(\"frontier://FrontierProd/CMS_CONDITIONS\")))"
# \nprocess.ALCARECOPPSCalMaxTracksFilter.TriggerResultsTag = cms.InputTag(\"TriggerResults\",\"\",\"HLT\")"

INPUTFILE="/store/data/Run2022C/ZeroBias/RAW/v1/000/356/615/00000/a4076bff-3bec-459e-bd44-b66538bf5c2b.root"
# INPUTFILE="/store/group/alca_global/pps_alcareco_producer_tests/outputALCAPPS_single.root"
COMMMAND=`xrdfs cms-xrd-global.cern.ch locate $INPUTFILE`
STATUS=$?
echo "xrdfs command status = "$STATUS
if [ $STATUS -eq 0 ]; then
    echo "Using file ${INPUTFILE}. Running in ${LOCAL_TEST_DIR}."
    (cmsDriver.py testExpressPPSAlCaRecoProducer -s ALCAPRODUCER:PPSCalMaxTracks,ENDJOB \
    --process ALCARECO \
    --scenario pp \
    --era Run3 \
    --conditions auto:run3_data_express \
    --data  \
    --datatier ALCARECO \
    --eventcontent ALCARECO \
    --nThreads 8 \
    --number 10000 --filein ${INPUTFILE} \
    --fileout file:outputALCAPPS_RECO_express.root \
    --customise_commands="$customise_commands") || die 'failed running test_express_AlCaRecoProducer' $?
else
    die "SKIPPING test, file ${INPUTFILE} not found" 0
fi
