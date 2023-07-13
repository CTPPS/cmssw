#!/bin/bash
function die { echo $1: status $2; exit $2; }

customise_commands="process.GlobalTag.toGet = cms.VPSet()\
\nprocess.GlobalTag.toGet.append(cms.PSet(record = cms.string(\"CTPPSOpticsRcd\"),tag =  cms.string(\"PPSOpticalFunctions_2023_v1_validation\"), connect = cms.string(\"frontier://FrontierProd/CMS_CONDITIONS\")))\
\nprocess.GlobalTag.toGet.append(cms.PSet(record = cms.string(\"LHCInfoRcd\"),tag =  cms.string(\"LHCInfo_PopCon_test\"), connect = cms.string(\"sqlite_file:/afs/cern.ch/user/a/akulczyc/public/DBsRecoTestFill9019/lhcinfo_pop_unit_test_old.db\")))\
\nprocess.GlobalTag.toGet.append(cms.PSet(record = cms.string(\"LHCInfoPerLSRcd\"),tag =  cms.string(\"ls_end_test\"), connect = cms.string(\"sqlite_file:/afs/cern.ch/user/a/akulczyc/public/DBsRecoTestFill9019/lhcinfo_pop_unit_test_new.db\")))\
\nprocess.GlobalTag.toGet.append(cms.PSet(record = cms.string(\"LHCInfoPerFillRcd\"),tag =  cms.string(\"fill_end_test\"), connect = cms.string(\"sqlite_file:/afs/cern.ch/user/a/akulczyc/public/DBsRecoTestFill9019/lhcinfo_pop_unit_test_new.db\")))"
# \nprocess.ALCARECOPPSCalMaxTracksFilter.TriggerResultsTag = cms.InputTag(\"TriggerResults\",\"\",\"HLT\")"

INPUTFILE="/store/data/Run2023D/ZeroBias/RAW/v1/000/369/956/00000/33d5acec-484f-4ac0-9b83-c6a3104ddd2b.root"
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
