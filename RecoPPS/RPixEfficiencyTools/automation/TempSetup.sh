# This scirpt is for a setup of environment variables. 
# The full setup is done by Jenkins script
# export WORKFLOW=pps_workflow_for_cmssw12
# export CAMPAIGN=mobrzut_test_1
# export DATASET=/EGamma/Run2018B-12Nov2019_UL2018-v2/AOD
# export PROXY=/afs/cern.ch/user/m/mobrzut/public/jenkins_proxy/jenkins.pem # should it be .x509 or .pem file  
export PYTHONPATH=\$PYTHONPATH:`pwd`/automation_control				  
source /cvmfs/cms.cern.ch/common/crab-setup.sh #todo why I still need to run this manually? 
