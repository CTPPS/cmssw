from CRABAPI.RawCommand import crabCommand
from CRABClient.UserUtilities import config

datasets = [
    # ZeroBias test
    # '/ZeroBias/Run2023B-22Sep2023-v1/MINIAOD',
    # '/ZeroBias/Run2023C-22Sep2023_v1-v1/MINIAOD',
    # '/ZeroBias/Run2023C-22Sep2023_v2-v1/MINIAOD',
    # '/ZeroBias/Run2023C-22Sep2023_v3-v1/MINIAOD',
    '/ZeroBias/Run2023C-22Sep2023_v4-v1/MINIAOD',
    # '/ZeroBias/Run2023D-22Sep2023_v1-v1/MINIAOD',
    # '/ZeroBias/Run2023D-22Sep2023_v2-v1/MINIAOD'

    # EGamma dataset
    # '/EGamma0/Run2023B-22Sep2023-v1/MINIAOD',
    # '/EGamma0/Run2023C-22Sep2023_v1-v1/MINIAOD',
    # '/EGamma0/Run2023C-22Sep2023_v2-v1/MINIAOD',
    # '/EGamma0/Run2023C-22Sep2023_v3-v1/MINIAOD',
    # '/EGamma0/Run2023C-22Sep2023_v4-v1/MINIAOD',
    # '/EGamma0/Run2023D-22Sep2023_v1-v1/MINIAOD',
    # '/EGamma0/Run2023D-22Sep2023_v2-v1/MINIAOD',
    # '/EGamma1/Run2023B-22Sep2023-v1/MINIAOD',
    # '/EGamma1/Run2023C-22Sep2023_v1-v1/MINIAOD',
    # '/EGamma1/Run2023C-22Sep2023_v2-v1/MINIAOD',
    # '/EGamma1/Run2023C-22Sep2023_v3-v1/MINIAOD',
    # '/EGamma1/Run2023C-22Sep2023_v4-v1/MINIAOD',
    # '/EGamma1/Run2023D-22Sep2023_v1-v1/MINIAOD',
    # '/EGamma1/Run2023D-22Sep2023_v2-v1/MINIAOD'
]

for ds in datasets:
    cfg = config()

    # General configuration
    cfg.General.requestName = 'LHCInfoPlotter' +ds.replace('/','_')
    cfg.General.workArea = 'crab_projects'
    cfg.General.transferOutputs = True
    cfg.General.transferLogs = True

    # JobType configuration
    cfg.JobType.pluginName = 'Analysis'
    cfg.JobType.psetName = 'testLHCInfoPlotter_cfg.py'

    # Pass arguments to VarParsing
    cfg.JobType.pyCfgParams = [
        'outputFileName=LHCInfoPlots.root',
        'globalTag=150X_dataRun3_v3',
        'hltPath=HLT_ZeroBias_v*',
        # 'hltPath=HLT_Ele32_WPTight_Gsf_v*',
        'crab=True'
    ]

    # Outputs to transfer
    cfg.JobType.outputFiles = ['LHCInfoPlots.root']

    # Data configuration
    cfg.Data.inputDataset = ds
    cfg.Data.inputDBS = 'global'
    cfg.Data.splitting = 'FileBased'
    cfg.Data.unitsPerJob = 15   # adjust depending on dataset size

    # This ensures lumi filtering consistency with your cfg
    cfg.Data.lumiMask = '/eos/project-c/ctpps/Operations/DataExternalConditions/2023/PPS-July2023_track_only_AND_CMS_Golden.json'
    cfg.Data.outLFNDirBase = '/store/group/phys_pps/abellora/dirsim/2023_lhcinfo'
    cfg.Data.publication = False
    
    # Site configuration
    cfg.Site.storageSite = 'T2_CH_CERN' 
    
    # Submit the job
    crabCommand('submit', config=cfg)