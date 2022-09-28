import CRABClient
from CRABClient.UserUtilities import config
config = config()

InputDataset ="/EGamma/Run2018B-12Nov2019_UL2018-v2/AOD"
number = 26

config.General.transferOutputs = True
config.General.transferLogs = True

config.General.requestName = f'mobrzut_test_EA_DQM_Worker_{number}'
config.General.workArea = '/afs/cern.ch/user/m/mobrzut/automation/environment_based_on_jenkins_script/pps_workflow_for_cmssw12/mobrzut_test_1/CMSSW_12_4_0/src/RecoPPS/RPixEfficiencyTools' # TODO: your path here



config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '/afs/cern.ch/user/m/mobrzut/automation/environment_based_on_jenkins_script/pps_workflow_for_cmssw12/mobrzut_test_1/CMSSW_12_4_0/src/RecoPPS/RPixEfficiencyTools/python/EfficiencyAnalysisDQMWorker_cfg.py' # TODO: your path here
config.JobType.pyCfgParams = ["sourceFileList=/afs/cern.ch/user/m/mobrzut/public/Era.dat", "outputFileName=tmp.root"]
config.Data.inputDataset = InputDataset

config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
# config.Data.splitting = 'Automatic'                                                                                         

config.Data.unitsPerJob = 20
config.Data.publication = False
config.Data.outLFNDirBase = '/store/group/dpg_ctpps/comm_ctpps/2018_PixelEfficiency'
config.Data.outputDatasetTag = f'CRAB3_tmobrzut_test_EA_DQM_Worker_{number}'
config.Data.runRange = '317080'


config.Site.storageSite = 'T2_CH_CERN'
