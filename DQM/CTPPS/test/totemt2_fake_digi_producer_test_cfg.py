# test to check if producer works correctly

import FWCore.ParameterSet.Config as cms

MAX_NUMBER_OF_EVENTS = 6000
PATH_TO_FAKE_T2_DATA = "/eos/home-a/acwikla/data/nT2_Nino.root"

process = cms.Process('TESTDIGIPRODUCER')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(MAX_NUMBER_OF_EVENTS) )

process.source = cms.Source("EmptySource")

process.verbosity = cms.untracked.PSet( input = cms.untracked.int32(-1) )

# minimum of logs
process.MessageLogger = cms.Service("MessageLogger",
    statistics = cms.untracked.vstring(),
    destinations = cms.untracked.vstring('cerr'),
    cerr = cms.untracked.PSet(
        threshold = cms.untracked.string('WARNING')
    )
)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')

ids = [
    # arm 0
    # plane 0
    # tile 0    tile 1      tile 2      tile 3
    2113929216, 2113945600, 2113961984, 2113978368,
    # plane 1
    2113994752, 2114011136, 2114027520, 2114043904,
    # plane 2
    2114060288, 2114076672, 2114093056, 2114109440,
    # plane 3
    2114125824, 2114142208, 2114158592, 2114174976,
    # plane 4
    2114191360, 2114207744, 2114224128, 2114240512,
    # plane 5
    2114256896, 2114273280, 2114289664, 2114306048,
    # plane 6
    2114322432, 2114338816, 2114355200, 2114371584,
    # plane 7
    2114387968, 2114404352, 2114420736, 2114437120,

    # arm 1
    # plane 0
    # tile 0    tile 1      tile 2      tile 3
    2130706432, 2130722816, 2130739200, 2130755584,
    # plane 1
    2130771968, 2130788352, 2130804736, 2130821120,
    # plane 2
    2130837504, 2130853888, 2130870272, 2130886656,
    # plane 3
    2130903040, 2130919424, 2130935808, 2130952192,
    # plane 4
    2130968576, 2130984960, 2131001344, 2131017728,
    # plane 5
    2131034112, 2131050496, 2131066880, 2131083264,
    # plane 6
    2131099648, 2131116032, 2131132416, 2131148800,
    # plane 7
    2131165184, 2131181568, 2131197952, 2131214336
]

ids_arm_0 = ids[:len(ids)//2]
ids_arm_1 = ids[len(ids)//2:]

ids_arm_0_plane_0 = ids_arm_0[:4]
ids_arm_1_plane_1 = ids_arm_1[4:8]

ids_arm_0_plane_all_tile_0 = ids_arm_0[0::4]
ids_arm_0_plane_all_tile_1 = ids_arm_0[1::4]
ids_arm_0_plane_all_tile_2 = ids_arm_0[2::4]
ids_arm_0_plane_all_tile_3 = ids_arm_0[3::4]

ids_arm_1_plane_all_tile_0 = ids_arm_1[0::4]
ids_arm_1_plane_all_tile_1 = ids_arm_1[1::4]
ids_arm_1_plane_all_tile_2 = ids_arm_1[2::4]
ids_arm_1_plane_all_tile_3 = ids_arm_1[3::4]


test_cases = cms.VPSet(
		cms.PSet(detId = cms.vuint32(*ids_arm_0_plane_all_tile_0), eventLimit = cms.uint32(100)),
		cms.PSet(detId = cms.vuint32(*ids_arm_0_plane_all_tile_1), eventLimit = cms.uint32(200)),
		cms.PSet(detId = cms.vuint32(*ids_arm_0_plane_all_tile_2), eventLimit = cms.uint32(500)),
		cms.PSet(detId = cms.vuint32(*ids_arm_0_plane_all_tile_3), eventLimit = cms.uint32(1000)),
        
		cms.PSet(detId = cms.vuint32(*ids_arm_1_plane_all_tile_0), eventLimit = cms.uint32(10)),
		cms.PSet(detId = cms.vuint32(*ids_arm_1_plane_all_tile_1), eventLimit = cms.uint32(50)),
		cms.PSet(detId = cms.vuint32(*ids_arm_1_plane_all_tile_2), eventLimit = cms.uint32(100)),
		cms.PSet(detId = cms.vuint32(*ids_arm_1_plane_all_tile_3), eventLimit = cms.uint32(500)),
)

# fake digi producer
process.totemT2Digis = cms.EDProducer('TotemT2DigiProducer',
	t2DataFile=cms.string(PATH_TO_FAKE_T2_DATA),
	testCasesSet = test_cases
)

# fake digi analyzer
process.totemT2DigiAnalyzer = cms.EDAnalyzer('TotemT2DigiAnalyzer',
    digisTag = cms.InputTag('totemT2Digis', 'TotemT2'),
    numberOfAllEvents = cms.uint32(MAX_NUMBER_OF_EVENTS),
    testCasesSet = test_cases
)

process.path = cms.Path(
    process.totemT2Digis *
    process.totemT2DigiAnalyzer
)

process.schedule = cms.Schedule(
    process.path
)
