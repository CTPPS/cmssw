import CRABClient
from CRABClient.UserUtilities import config
config = config()

# InputDataset ="/AlCaPPSPrompt/Run2022C-PPSCalMaxTracks-PromptReco-v1/ALCARECO"
InputDataset = "/AlCaPPSPrompt/Run2022C-PPSCalMaxTracks-PromptReco-v1/ALCARECO" # on german DC
requestNumber = 205

config.General.transferOutputs = True
config.General.transferLogs = True

config.General.requestName = f'DiamondTimingRes_mobrzut_{requestNumber}'



config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '/afs/cern.ch/user/m/mobrzut/public/timing_calibration_env/CMSSW_12_4_6/src/Analyzer/DiamondTimingAnalyzer/temp/ExampleConfig_cfg.py'
# config.JobType.pyCfgParams = ["sourceFileList=/afs/cern.ch/user/m/mobrzut/public/Era.dat", "outputFileName=tmp.root"]
config.Data.inputDataset = InputDataset
config.JobType.outputFiles = ["run_output.root"]
config.JobType.priority = 40

config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 1000
config.Data.publication = False
config.Data.outLFNDirBase = '/store/group/dpg_ctpps/comm_ctpps/2018_PixelEfficiency'
config.Data.outputDatasetTag = f'DiamondTimingRes_mobrzut_{requestNumber}'
config.Data.runRange = '357442'
config.Data.lumiMask ='/afs/cern.ch/user/m/mobrzut/public/allrunsSB-PPS-forCalib.json'


config.Site.storageSite = 'T2_CH_CERN'
