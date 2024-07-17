process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '124X_dataRun3_v9', '')

process.load("CalibPPS.TimingCalibration.ppsTimingCalibrationPCLWorker_cfi")
process.ppsTimingCalibrationPCLWorker.diamondRecHitTags=cms.InputTag("ctppsDiamondRecHits","","PPSTiming2")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        "/store/data/Run2022C/AlCaPPSPrompt/ALCARECO/PPSCalMaxTracks-PromptReco-v1/000/357/440/00000/53746080-d21b-4fdf-9d19-67c6dae347ac.root",
"/store/data/Run2022C/AlCaPPSPrompt/ALCARECO/PPSCalMaxTracks-PromptReco-v1/000/357/440/00000/b8d3f012-0b59-43b4-adbc-811bcb59c9c4.root",
"/store/data/Run2022C/AlCaPPSPrompt/ALCARECO/PPSCalMaxTracks-PromptReco-v1/000/357/440/00000/cfa9f82a-5296-49cf-8353-11f8761b675c.root",
"/store/data/Run2022C/AlCaPPSPrompt/ALCARECO/PPSCalMaxTracks-PromptReco-v1/000/357/440/00000/db3c6c33-9ad0-4d61-9a64-a1aa9578550e.root"
    # 'file:LOCAL_FILE_NAME', #local file have to be declared with file prefix
    # 'REMOTE_FILE_PATH' #but remote files - for example T2_UK have to be declared without the prefix
))

# process.load("EventFilter.CTPPSRawToDigi.ctppsRawToDigi_cff")
# process.load("RecoPPS.Configuration.recoCTPPS_cff")
# 

# process.ctppsDiamondRawToDigi.rawDataTag = cms.InputTag("INSERT_PRODUCT_TAG")

process.path = cms.Path(
    # process.a1* #when running on the streamer files
    # process.ctppsRawToDigi *
    # process.recoCTPPS *
    process.ppsTimingCalibrationPCLWorker
)