from CRABAPI.RawCommand import crabCommand
from CRABClient.UserUtilities import config

datasets = [
    # ZeroBias test
    # '/ZeroBias/Run2023B-22Sep2023-v1/MINIAOD',
    # '/ZeroBias/Run2023C-22Sep2023_v1-v1/MINIAOD',
    # '/ZeroBias/Run2023C-22Sep2023_v2-v1/MINIAOD',
    # '/ZeroBias/Run2023C-22Sep2023_v3-v1/MINIAOD',
    # '/ZeroBias/Run2023C-22Sep2023_v4-v1/MINIAOD',
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
    
    # AlCaPPS dataset
    "/AlCaPPSPrompt/Run2023B-PPSCalMaxTracks-PromptReco-v1/ALCARECO",
    # "/AlCaPPSPrompt/Run2023C-PPSCalMaxTracks-PromptReco-v1/ALCARECO",
    # "/AlCaPPSPrompt/Run2023C-PPSCalMaxTracks-PromptReco-v2/ALCARECO",
    # "/AlCaPPSPrompt/Run2023C-PPSCalMaxTracks-PromptReco-v3/ALCARECO",
    # "/AlCaPPSPrompt/Run2023C-PPSCalMaxTracks-PromptReco-v4/ALCARECO",
    # "/AlCaPPSPrompt/Run2023D-PPSCalMaxTracks-PromptReco-v1/ALCARECO",
    # "/AlCaPPSPrompt/Run2023D-PPSCalMaxTracks-PromptReco-v2/ALCARECO",
]

for ds in datasets:
    cfg = config()

    # General configuration
    cfg.General.requestName = 'ProtonRecoPlotters' +ds.replace('/','_')+'_xangle_159_161'
    cfg.General.workArea = 'crab_projects'
    cfg.General.transferOutputs = True
    cfg.General.transferLogs = True

    # JobType configuration
    cfg.JobType.pluginName = 'Analysis'
    cfg.JobType.psetName = 'testProtonRecoPlotters_cfg.py'

    # Pass arguments to VarParsing
    cfg.JobType.pyCfgParams = [
        'globalTag=150X_dataRun3_v5',
        # 'hltPath=HLT_ZeroBias_v*',
        # 'hltPath=HLT_Ele32_WPTight_Gsf_v*',
        'hltPath=HLT_PPSMaxTracksPerRP4_v*',
        'minCrossingAngle=159',
        'maxCrossingAngle=161',
        'crab=True'
    ]

    # Outputs to transfer
    cfg.JobType.outputFiles = ['LHCInfoPlots.root', 'ProtonRecoPlots.root']

    # Data configuration
    cfg.Data.inputDataset = ds
    cfg.Data.inputDBS = 'global'
    cfg.Data.splitting = 'FileBased'
    cfg.Data.unitsPerJob = 15   # adjust depending on dataset size

    # This ensures lumi filtering consistency with your cfg
    cfg.Data.lumiMask = '/eos/project-c/ctpps/Operations/DataExternalConditions/2023/PPS-July2023_track_only_AND_CMS_Golden.json'
    cfg.Data.outLFNDirBase = '/store/group/phys_pps/abellora/dirsim/2023_protonreco'
    cfg.Data.publication = False
    
    # Site configuration
    cfg.Site.storageSite = 'T2_CH_CERN' 
    
    # Submit the job
    crabCommand('submit', config=cfg)