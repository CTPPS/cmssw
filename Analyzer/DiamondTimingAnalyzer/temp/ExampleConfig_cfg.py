
import FWCore.ParameterSet.Config as cms
import copy

process = cms.Process('PPSTiming2')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(30000)
)

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = ''
process.MessageLogger.cerr.FwkReport.reportEvery = 1000




process.source = cms.Source ("PoolSource",
                             fileNames = cms.untracked.vstring(
    # "file:/eos/home-m/mobrzut/357/734/00000/73c5eeec-578d-438f-92f1-3e1fff2ed390.root",   
    # "file:/eos/home-m/mobrzut/357/734/00000/08e5ee70-ebff-4571-b14c-806961a3d9dd.root",
    # "file:/eos/home-m/mobrzut/357/734/00000/d42eaf2c-57bb-48fe-9a22-317053f89115.root"
"file:/eos/home-m/mobrzut/357/479/00000/08527d3f-280d-4e30-8a6c-59e0e49e229f.root",
"file:/eos/home-m/mobrzut/357/479/00000/08e370cc-8447-4dad-8b03-37c69a9c4a9a.root",
"file:/eos/home-m/mobrzut/357/479/00000/0b666b18-e81a-4761-85d3-e7933739a707.root",
"file:/eos/home-m/mobrzut/357/479/00000/118f67aa-1bf5-44bd-bb8e-72f6b1e5199f.root",
"file:/eos/home-m/mobrzut/357/479/00000/12ff4fd8-3deb-43c6-b4d0-f543a832e0e4.root",
"file:/eos/home-m/mobrzut/357/479/00000/144ac447-939f-49ad-aa12-0611b5ddc59f.root",
"file:/eos/home-m/mobrzut/357/479/00000/20bc35c7-8f38-45b2-b78c-1096ac3ec406.root",
"file:/eos/home-m/mobrzut/357/479/00000/21533e2b-671d-497a-ad6b-7c3475c76373.root",
"file:/eos/home-m/mobrzut/357/479/00000/23bfd3d4-c39e-4764-ac77-068522debdb6.root",
"file:/eos/home-m/mobrzut/357/479/00000/288b5e2f-175f-45cf-b3d0-663a48784a79.root",
"file:/eos/home-m/mobrzut/357/479/00000/28a54f2b-0c7f-4bf2-9c8e-9a1a54b1068f.root",
"file:/eos/home-m/mobrzut/357/479/00000/28c05462-ba33-4d96-9003-7f8e7bbda5be.root",
"file:/eos/home-m/mobrzut/357/479/00000/2cbeab2b-242a-4166-801f-47aa219ac337.root",
"file:/eos/home-m/mobrzut/357/479/00000/2dc99004-1314-4d6e-8772-2c8162ba0b39.root",
"file:/eos/home-m/mobrzut/357/479/00000/36d49a7d-933f-4b1f-b5cd-8b17d1602fc3.root",
"file:/eos/home-m/mobrzut/357/479/00000/3763c457-56a6-4c5e-b8db-da118a9b64a6.root",
"file:/eos/home-m/mobrzut/357/479/00000/3a78ccd5-d768-4fef-8c19-a883f7235f49.root",
"file:/eos/home-m/mobrzut/357/479/00000/4077ebaa-ed45-466e-8475-58c1e7c7788e.root",
"file:/eos/home-m/mobrzut/357/479/00000/412550de-2e35-4937-a255-f62bb411adb9.root",
"file:/eos/home-m/mobrzut/357/479/00000/4636f253-374e-408c-8626-4a40ec1c4cdd.root",
"file:/eos/home-m/mobrzut/357/479/00000/46e829bb-f69c-42f9-9a2d-1f9b85d68f54.root",
"file:/eos/home-m/mobrzut/357/479/00000/487cd5dc-364d-4c10-a133-80044105470a.root",
"file:/eos/home-m/mobrzut/357/479/00000/49f659be-9825-4df9-be8b-03fe2ef26bf4.root",
"file:/eos/home-m/mobrzut/357/479/00000/4b63445d-a6a1-43f1-b691-16f77464a30b.root",
"file:/eos/home-m/mobrzut/357/479/00000/53a43613-e1fb-49b0-b1a0-877a701dc9e1.root",
"file:/eos/home-m/mobrzut/357/479/00000/564cd45c-8486-4d6c-88dc-aeee534ca177.root",
"file:/eos/home-m/mobrzut/357/479/00000/57168ee6-c48d-4046-8814-e2d9db8b35b3.root",
"file:/eos/home-m/mobrzut/357/479/00000/57dcd973-5eb6-4aa9-a04a-b8ba8ef26099.root",
"file:/eos/home-m/mobrzut/357/479/00000/626a3f06-7db8-4d77-9434-f6d71906ab8f.root",
"file:/eos/home-m/mobrzut/357/479/00000/63e4e673-83e2-4482-a425-f3fefe9a7f14.root",
"file:/eos/home-m/mobrzut/357/479/00000/68a5a64a-b756-4cf5-9a2a-00181e34f501.root",
"file:/eos/home-m/mobrzut/357/479/00000/69634dd6-28bd-4c7a-bb26-16b60cd3c963.root",
"file:/eos/home-m/mobrzut/357/479/00000/6bf1c57b-8a37-4d25-9f7e-100383e31e89.root",
"file:/eos/home-m/mobrzut/357/479/00000/71c9c5eb-0510-41b2-bb86-70d406e311c9.root",
"file:/eos/home-m/mobrzut/357/479/00000/74586414-5bdf-4b8e-9150-38992711303f.root",
"file:/eos/home-m/mobrzut/357/479/00000/7542d52b-b72f-4d14-b715-68a634d91cf0.root",
"file:/eos/home-m/mobrzut/357/479/00000/80e68cf8-f666-4193-9216-25e7ceecaa1b.root",
"file:/eos/home-m/mobrzut/357/479/00000/8115b18d-1667-4b68-a757-4c6ffcb4d2da.root",
"file:/eos/home-m/mobrzut/357/479/00000/82329b1d-38a7-41fd-8dd2-182527ebd92a.root",
"file:/eos/home-m/mobrzut/357/479/00000/8344c1dd-2d5d-4478-92f8-718073deffd2.root",
"file:/eos/home-m/mobrzut/357/479/00000/83546199-aabb-4cd6-8a65-49b9cda74287.root",
"file:/eos/home-m/mobrzut/357/479/00000/84344b82-5bd9-4a75-a1c9-cc6115502d8e.root",
"file:/eos/home-m/mobrzut/357/479/00000/84ee7daa-082b-41cd-bb38-5319b8d4477a.root",
"file:/eos/home-m/mobrzut/357/479/00000/85cdc11d-7278-49d0-be1f-d79669a3e52c.root",
"file:/eos/home-m/mobrzut/357/479/00000/85f91b64-0ccf-419e-a811-34e2f82daf85.root",
"file:/eos/home-m/mobrzut/357/479/00000/86548615-30f7-44fa-a26d-e2d6e4c2d491.root",
"file:/eos/home-m/mobrzut/357/479/00000/876635f4-a030-4fd0-ab01-c56ee52e08de.root",
"file:/eos/home-m/mobrzut/357/479/00000/8c26f951-e390-49f4-a00b-4e679444efca.root",
"file:/eos/home-m/mobrzut/357/479/00000/8d050fc6-ee48-403f-ab93-9db0e4e90dab.root",
"file:/eos/home-m/mobrzut/357/479/00000/8d99bb61-3482-405d-b097-06e33d7489bd.root",
"file:/eos/home-m/mobrzut/357/479/00000/8e22fff6-0083-4978-9e86-d3764e1a794c.root",
"file:/eos/home-m/mobrzut/357/479/00000/97d96fb2-cdcd-42f9-ad98-30b696baa99b.root",
"file:/eos/home-m/mobrzut/357/479/00000/99215ea4-0f71-4e1e-b109-764db7d48fb9.root",
"file:/eos/home-m/mobrzut/357/479/00000/994e839c-72f6-4664-b5fd-e952d90f1038.root",
"file:/eos/home-m/mobrzut/357/479/00000/9a0f6223-a50c-4637-8e5e-38f9588528db.root",
"file:/eos/home-m/mobrzut/357/479/00000/9a203008-1da8-4c97-a1f8-2eea98848997.root",
"file:/eos/home-m/mobrzut/357/479/00000/9aefd88f-7eb1-431e-ae70-3c0a7a6d1106.root",
"file:/eos/home-m/mobrzut/357/479/00000/9c9f172e-9cb5-4b33-88cb-78660ee97c4b.root",
"file:/eos/home-m/mobrzut/357/479/00000/9e5d933e-290d-4ffd-adf1-d35a5fb05dae.root",
"file:/eos/home-m/mobrzut/357/479/00000/a3dfad41-8bf8-41f2-a3fd-1b82c4f4c413.root",
"file:/eos/home-m/mobrzut/357/479/00000/a450432a-1f65-424b-aa19-aaad1fff258c.root",
"file:/eos/home-m/mobrzut/357/479/00000/aa4bb8fd-e03f-4fb9-ba78-18aea70fe277.root",
"file:/eos/home-m/mobrzut/357/479/00000/b2c0fcb2-8ed7-459d-af57-09b74ca250f3.root",
"file:/eos/home-m/mobrzut/357/479/00000/b4654b43-e800-412b-afed-7f64657e3af1.root",
"file:/eos/home-m/mobrzut/357/479/00000/b475fca2-bced-4ee1-9aa9-cb15c4ae1a62.root",
"file:/eos/home-m/mobrzut/357/479/00000/b4ef7261-4282-498a-9e74-0588d29fba97.root",
"file:/eos/home-m/mobrzut/357/479/00000/b7448794-d7ad-43df-bdd4-1f0edbf00824.root",
"file:/eos/home-m/mobrzut/357/479/00000/bcce134f-be55-4fc7-860d-307d8fe99b50.root",
"file:/eos/home-m/mobrzut/357/479/00000/bfb9965d-1c77-4ed7-8be5-0fbd091d5b16.root",
"file:/eos/home-m/mobrzut/357/479/00000/c1745025-633d-426a-a3a4-abb6138c81be.root",
"file:/eos/home-m/mobrzut/357/479/00000/c1db7e32-1e39-43b3-8c96-b42b8bae2303.root",
"file:/eos/home-m/mobrzut/357/479/00000/c7a3b355-af9f-4312-81b9-013d80bcfd8f.root",
"file:/eos/home-m/mobrzut/357/479/00000/c804e4b2-fc97-44fa-8bac-466c514ca149.root",
"file:/eos/home-m/mobrzut/357/479/00000/c8c8e55a-8126-4007-ae53-26fee9b611c8.root",
"file:/eos/home-m/mobrzut/357/479/00000/c8e6fab9-13d5-4df0-83e5-ded0443cb7f6.root",
"file:/eos/home-m/mobrzut/357/479/00000/ca9a67b7-b672-4d3a-9b32-b5dca0266f3f.root",
"file:/eos/home-m/mobrzut/357/479/00000/d19c8512-9fe0-4007-97a0-c84a29dbb1ce.root",
"file:/eos/home-m/mobrzut/357/479/00000/d4118130-d057-4fbe-9ea1-63e411df25a1.root",
"file:/eos/home-m/mobrzut/357/479/00000/d63b0bdf-8443-49e1-a3d9-99fce9b86833.root",
"file:/eos/home-m/mobrzut/357/479/00000/d7040537-edaa-4fe9-bd87-29b44bef4973.root",
"file:/eos/home-m/mobrzut/357/479/00000/d8c59cf0-67a3-4310-94f7-79e49b464d83.root",
"file:/eos/home-m/mobrzut/357/479/00000/dbda4868-8d7f-4111-b7d0-e002543143c7.root",
"file:/eos/home-m/mobrzut/357/479/00000/dca4f0b6-5a72-4b22-a77e-0d4aac589b17.root",
"file:/eos/home-m/mobrzut/357/479/00000/de0844d5-20af-4eae-a01a-c2a9170c4f98.root",
"file:/eos/home-m/mobrzut/357/479/00000/e0232315-a0e9-4c2d-bd92-91bcaa5f3e54.root",
"file:/eos/home-m/mobrzut/357/479/00000/e31ebdd9-daed-4666-a198-45dd7a5717a4.root",
"file:/eos/home-m/mobrzut/357/479/00000/e427550e-c886-4ef5-a357-19153639c140.root",
"file:/eos/home-m/mobrzut/357/479/00000/e490c59f-0cc6-436b-af9b-687ff2e6554f.root",
"file:/eos/home-m/mobrzut/357/479/00000/e52992c5-b43b-4c66-9624-4c2b2c6850d6.root",
"file:/eos/home-m/mobrzut/357/479/00000/e7d08e18-8ead-4e68-b339-9e21988d2624.root",
"file:/eos/home-m/mobrzut/357/479/00000/f01a091a-6171-4f5f-958e-a08dc64887f0.root",
"file:/eos/home-m/mobrzut/357/479/00000/f06f5294-d0b9-434c-9aac-2d6ddf6b74e8.root",
"file:/eos/home-m/mobrzut/357/479/00000/f2bab587-9806-45d3-9f84-636d8a2d52b3.root",
"file:/eos/home-m/mobrzut/357/479/00000/f7906aa9-ee7d-44b4-a529-f6b8b5b14543.root",
"file:/eos/home-m/mobrzut/357/479/00000/fa423635-30d1-44e3-b5c6-8708df198434.root",
"file:/eos/home-m/mobrzut/357/479/00000/fab66f14-5e52-4e11-b732-d7353714ede7.root",
"file:/eos/home-m/mobrzut/357/479/00000/fc155bd5-6346-4b10-86e3-5a1801da2fb0.root",
"file:/eos/home-m/mobrzut/357/479/00000/fc247753-aeae-4a03-b930-9af9d59d0808.root",
"file:/eos/home-m/mobrzut/357/479/00000/fc368df8-6037-4966-aa6e-5d78f4d2c262.root" 

)
)

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
# 2022 prompt: to be updated
process.GlobalTag = GlobalTag(process.GlobalTag, "124X_dataRun3_v9") 


#JH - use new tag for timing calibrations



                                        
process.load('CondCore.CondDB.CondDB_cfi')
process.CondDB.connect = 'sqlite_file:/afs/cern.ch/user/m/mobrzut/automation/timing_calibration/CMSSW_12_4_6/src/Analyzer/DiamondTimingAnalyzer/temp/ppsDiamondTiming_calibration357734.sqlite' # SQLite input
process.PoolDBESSource = cms.ESSource('PoolDBESSource',
        process.CondDB,
        DumpStats = cms.untracked.bool(True),
        toGet = cms.VPSet(
            cms.PSet(
                record = cms.string('PPSTimingCalibrationRcd'),
                tag = cms.string('DiamondTimingCalibration'),
                           label = cms.untracked.string('PPSTestCalibration'), 
            )
        )
    )

apply_shift = True                                        

# JH - rerun reco sequence with new timing conditions
process.load("RecoPPS.Configuration.recoCTPPS_cff")
process.ctppsDiamondLocalTracks.recHitsTag = cms.InputTag("ctppsDiamondRecHits","","PPSTiming2")
process.ctppsLocalTrackLiteProducer.tagDiamondTrack = cms.InputTag("ctppsDiamondLocalTracks","","PPSTiming2")
process.ctppsProtons.tagLocalTrackLite = cms.InputTag("ctppsLocalTrackLiteProducer","","PPSTiming2")
process.ctppsLocalTrackLiteProducer.includeDiamonds = cms.bool(True)
process.ctppsLocalTrackLiteProducer.tagPixelTrack = cms.InputTag("ctppsPixelLocalTracksAlCaRecoProducer")
process.ctppsLocalTrackLiteProducer.includePixels = cms.bool(True)
process.ctppsDiamondRecHits.timingCalibrationTag=cms.string("PoolDBESSource:PPSTestCalibration")
process.ctppsDiamondRecHits.digiTag= cms.InputTag("ctppsDiamondRawToDigiAlCaRecoProducer:TimingDiamond")

from DQMServices.Core.DQMEDAnalyzer import DQMEDAnalyzer
if(apply_shift):
    process.diamondTimingWorker = DQMEDAnalyzer("DiamondTimingWorker",
        tagDigi = cms.InputTag("ctppsDiamondRawToDigiAlCaRecoProducer", "TimingDiamond"),
        tagRecHit = cms.InputTag("ctppsDiamondRecHits","","PPSTiming2"),
        tagPixelLocalTrack = cms.InputTag("ctppsPixelLocalTracksAlCaRecoProducer"),
        timingCalibrationTag=cms.string("PoolDBESSource:PPSTestCalibration"),
        tagLocalTrack = cms.InputTag("ctppsDiamondLocalTracks","","PPSTiming2"),
        tagValidOOT = cms.int32(-1),
        planesConfig = cms.string("planes.json"),
        Ntracks_Lcuts = cms.vint32([-1,1,-1,1]), # minimum number of tracks in pots [45-210, 45-220, 56-210, 56-220]
        Ntracks_Ucuts = cms.vint32([-1,6,-1,6]), # maximum number of tracks in pots [45-210, 45-220, 56-210, 56-220]
    )
else:
    process.diamondTimingWorker = DQMEDAnalyzer("DiamondTimingWorker",
    tagDigi = cms.InputTag("ctppsDiamondRawToDigiAlCaRecoProducer", "TimingDiamond"),
    tagRecHit = cms.InputTag("ctppsDiamondRecHitsAlCaRecoProducer"),
    tagPixelLocalTrack = cms.InputTag("ctppsPixelLocalTracksAlCaRecoProducer"),
    tagLocalTrack = cms.InputTag("ctppsDiamondLocalTracksAlCaRecoProducer"),
    timingCalibrationTag=cms.string("PoolDBESSource:PPSTestCalibration"),
    tagValidOOT = cms.int32(-1),
    planesConfig = cms.string("planes.json"),
    Ntracks_Lcuts = cms.vint32([-1,1,-1,1]), # minimum number of tracks in pots [45-210, 45-220, 56-210, 56-220]
    Ntracks_Ucuts = cms.vint32([-1,6,-1,6]), # maximum number of tracks in pots [45-210, 45-220, 56-210, 56-220]
) 

if(apply_shift):
    process.ALL = cms.Path(
        # Re-run the PPS local+proton reconstruction starting from AOD with new timing calibrations
        process.ctppsDiamondRecHits *
        process.ctppsDiamondLocalTracks *
        process.diamondTimingWorker
                        )
else: 
    process.ALL = cms.Path(process.diamondTimingWorker)

saveToDQM = False
if(saveToDQM):
    process.load("DQMServices.Core.DQMStore_cfi")
    process.load("DQMServices.Components.DQMEnvironment_cfi")
    process.dqmEnv.subSystemFolder = "CalibPPS"
    process.dqmSaver.convention = 'Offline'
    process.dqmSaver.workflow = "/CalibPPS/TimingCalibration/CMSSW_12_0_0_pre2"
    process.dqmSaver.saveByRun = -1
    process.dqmSaver.saveAtJobEnd = True
    process.dqmSaver.forceRunNumber = 999999


    process.dqmOutput = cms.OutputModule("DQMRootOutputModule",
        fileName = cms.untracked.string("worker_output.root")
    )

    process.load("DQMServices.Components.EDMtoMEConverter_cff")
    process.end_path = cms.EndPath(process.dqmSaver)
else:
    process.dqmOutput = cms.OutputModule("DQMRootOutputModule", fileName=cms.untracked.string("temp_run_output_just_diamond.root"))
    process.end_path = cms.EndPath(process.dqmOutput)


process.schedule = cms.Schedule(process.ALL,   process.end_path)

#print(process.dumpPython())
