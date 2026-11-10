import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
from Configuration.StandardSequences.Eras import eras

process = cms.Process('ProtonRecoPlotter', eras.Run3_2023)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
    )

options = VarParsing.VarParsing ('analysis')

options.register('lhcInfoOutputFilename',
                'LHCInfoPlots.root',
                VarParsing.VarParsing.multiplicity.singleton,
                VarParsing.VarParsing.varType.string,
                "output ROOT file name")
options.register('protonOutputFilename',
                'ProtonRecoPlots.root',
                VarParsing.VarParsing.multiplicity.singleton,
                VarParsing.VarParsing.varType.string,
                "output ROOT file name")
options.register('jsonFileName',
                '',
                VarParsing.VarParsing.multiplicity.singleton,
                VarParsing.VarParsing.varType.string,
                "JSON file list name")
options.register('globalTag',
                '',
                VarParsing.VarParsing.multiplicity.singleton,
                VarParsing.VarParsing.varType.string,
                "GT to use")
options.register('hltPath',
                '',
                VarParsing.VarParsing.multiplicity.singleton,
                VarParsing.VarParsing.varType.string,
                "HLT path to use for event selection")
options.register('minCrossingAngle',
                 -1,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.float,
                 "Minimum crossing angle")
options.register('maxCrossingAngle',
                 -1,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.float,
                 "Maximum crossing angle")
options.register('crab',
                False,
                VarParsing.VarParsing.multiplicity.singleton,
                VarParsing.VarParsing.varType.bool,
                "Flag to indicate if running in crab")
options.register('nanoOutputFilename', 'PPSNanoAOD.root',
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.string, "nanoAOD output file name")
options.register('miniOutputFilename',
                 'PPSMiniAOD.root',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "MiniAOD output file name")
options.parseArguments()

# Limit the number of events
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

# Limit the logging
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

# Set INFO threshold for logs
process.MessageLogger.cerr.threshold = 'WARNING'
# process.MessageLogger.cerr.INFO.limit = -1

from Configuration.AlCa.GlobalTag import GlobalTag
if options.globalTag != '':
    gt = options.globalTag
else:
    gt = 'auto:run3_data_prompt'

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
print('Using GT:',gt)
process.GlobalTag = GlobalTag(process.GlobalTag, gt)

process.load('Geometry.VeryForwardGeometry.geometryRPFromDB_cfi')

# Check if input files are provided via command line arguments
if not options.crab:
    if not options.inputFiles:
        print('No input files provided. Please specify input files using inputFiles.')
        exit(1)
    fileList = [f'file:{f}' if not (f.startswith('/store/') or f.startswith('file:') or f.startswith('root:')) else f for f in options.inputFiles]
    inputFiles = cms.untracked.vstring(fileList)
    if len(inputFiles) == 0:
        print('No input files provided. Please specify input files using inputFiles.')
        exit(1)
    print('Input files:')
    print(inputFiles)

    process.source = cms.Source("PoolSource",
        fileNames = inputFiles,
        dropDescendantsOfDroppedBranches=cms.untracked.bool(False),
        inputCommands = cms.untracked.vstring(
            'drop *',
            'keep *_TriggerResults_*_*',
            'keep *_*Digi*_*_RECO',
            )
    )

    if options.jsonFileName:
        import FWCore.PythonUtilities.LumiList as LumiList
        jsonFileName = options.jsonFileName
        print('Using JSON file:', jsonFileName)
        process.source.lumisToProcess = LumiList.LumiList(filename = jsonFileName).getVLuminosityBlockRange()
else:
    process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring(),  # will be overridden by crab with the actual input files
        dropDescendantsOfDroppedBranches=cms.untracked.bool(False),
        inputCommands = cms.untracked.vstring(
            'drop *',
            'keep *_TriggerResults_*_*',
            'keep *_*Digi*_*_RECO',)
    )

# Temporarily manually implement the association cuts
pAssociationCuts = cms.PSet(
    validityRange=cms.EventRange("343890:min - 999999:max"),
    association_cuts_45=cms.PSet(
    # Filled from the latest analyzer output by update_assocuts_from_analyze.py.
        x_cut_mean=cms.string("- ( (0.1886599004902902+-0.019563398876999507 * [x_near]+0.0007653309974208171 * [x_near] * [x_near]) - (-0.028450862819029463+0.6299345376251246/[x_near]) * abs([y_near]) )"),
        x_cut_threshold=cms.string("1.5 * ( (0.3336556715260207+-0.00019050370618842877 * [x_near]) + (0.008737264460198742) * [y_near] * [y_near] )"),
        y_cut_mean=cms.string("- ( (3.826046933499323+-23.167508663013848/[x_near]) * TMath::Erf([y_near]/(-54.51529256817561+8.090154229044156 * [x_near])) )"),
        y_cut_threshold=cms.string("1.5 * ( (0.18530153234402777+2.2258475043878354/[x_near]/[x_near]) + (0.042701028234240955+91.33855712506406/pow([x_near],5)) * abs([y_near]) )"),
        ti_tr_min=cms.double(-1.5), 
        ti_tr_max=cms.double(2.0),
    ),
    association_cuts_56=cms.PSet(
        x_cut_mean=cms.string("- ( (0.06371888416162158+-0.019570009647643494 * [x_near]+0.0006094070505282202 * [x_near] * [x_near]) - (-0.028531198571901506+0.7118042948348349/[x_near]) * abs([y_near]) )"),
        x_cut_threshold=cms.string("1.5 * ( (0.2018301706199127+0.007529857241629447 * [x_near]) + (0.008666130062040476) * [y_near] * [y_near] )"),
        y_cut_mean=cms.string("- ( (1.68211583240235+-7.853516816096297/[x_near]) * TMath::Erf([y_near]/(-9.83022170905736+2.1434461875569877 * [x_near])) )"),
        y_cut_threshold=cms.string("1.5 * ( (0.15517271960650353+0.5811217032099263/[x_near]/[x_near]) + (0.017920830203140732+128.60886741524922/pow([x_near],5)) * abs([y_near]) )"),
        ti_tr_min=cms.double(-1.5), 
        ti_tr_max=cms.double(2.0),
    ),
)
process.load("CalibPPS.ESProducers.ppsAssociationCutsESSource_cfi")
process.ppsAssociationCutsESSource.configuration.append(pAssociationCuts)
process.es_prefer_ppsAssociationCutsESSource = cms.ESPrefer("PPSAssociationCutsESSource", "ppsAssociationCutsESSource")

# Load PPS AlcaReco to run re-reco
process.load("Calibration.PPSAlCaRecoProducer.ALCARECOPPSCalMaxTracks_cff")
process.recoPPSSequenceAlCaRecoProducer.remove(process.diamondSampicLocalReconstructionTaskAlCaRecoProducer)

process.load('Validation.CTPPS.ctppsLHCInfoPlotter_cfi')
process.ctppsLHCInfoPlotter.useNewLHCInfo = True
process.ctppsLHCInfoPlotter.outputFile = options.lhcInfoOutputFilename

process.load('Validation.CTPPS.ctppsLHCInfoFilter_cfi')
process.ctppsLHCInfoFilter.useNewLHCInfo = True
process.ctppsLHCInfoFilter.minCrossingAngle = options.minCrossingAngle # urad
process.ctppsLHCInfoFilter.maxCrossingAngle = options.maxCrossingAngle # urad
    
process.ctppsProtonReconstructionPlotter = cms.EDAnalyzer("CTPPSProtonReconstructionPlotter",
    tagTracks = cms.InputTag("ctppsLocalTrackLiteProducerAlCaRecoProducer","","ProtonRecoPlotter"),
    tagRecoProtonsSingleRP = cms.InputTag("ctppsProtonsAlCaRecoProducer", "singleRP","ProtonRecoPlotter"),
    tagRecoProtonsMultiRP = cms.InputTag("ctppsProtonsAlCaRecoProducer", "multiRP","ProtonRecoPlotter"),

    rpId_45_F = cms.uint32(23),
    rpId_45_N = cms.uint32(3),
    rpId_56_N = cms.uint32(103),
    rpId_56_F = cms.uint32(123),

    outputFile = cms.string(options.protonOutputFilename),
)

# PPS proton nanoAOD content
process.load('PhysicsTools.NanoAOD.protons_cff')

process.protonTable = process.protonTable.clone(
    tagRecoProtonsSingle = cms.InputTag("ctppsProtonsAlCaRecoProducer", "singleRP","ProtonRecoPlotter"),
    tagRecoProtonsMulti  = cms.InputTag("ctppsProtonsAlCaRecoProducer", "multiRP","ProtonRecoPlotter"),
    tagTrackLite         = cms.InputTag("ctppsLocalTrackLiteProducerAlCaRecoProducer","","ProtonRecoPlotter"),
)
process.multiRPTable = process.multiRPTable.clone(
    src = cms.InputTag("ctppsProtonsAlCaRecoProducer", "multiRP","ProtonRecoPlotter"),
)
process.singleRPTable = process.singleRPTable.clone(
    src = cms.InputTag("ctppsProtonsAlCaRecoProducer", "singleRP","ProtonRecoPlotter"),
)

# LHCInfo nanoAOD content
process.load('PhysicsTools.NanoAOD.lhcInfoProducer_cfi')
process.lhcInfoTable = process.lhcInfoProducer.clone()
(~eras.run3_common).toModify(
    process.lhcInfoTable, useNewLHCInfo=False
)

if options.hltPath:
    print('Using HLT path:', options.hltPath)
    process.load('HLTrigger.HLTfilters.hltHighLevel_cfi')
    process.hltHighLevel.HLTPaths = [options.hltPath]
    process.p = cms.Path(process.hltHighLevel * 
                         process.ctppsLHCInfoFilter *
                         process.recoPPSSequenceAlCaRecoProducer *
                         process.ctppsLHCInfoPlotter *
                         process.ctppsProtonReconstructionPlotter *
                         process.lhcInfoTable,
                         process.protonTablesTask 
    )
else:
    process.p = cms.Path(process.ctppsLHCInfoFilter *
                         process.recoPPSSequenceAlCaRecoProducer *
                         process.ctppsLHCInfoPlotter *
                         process.ctppsProtonReconstructionPlotter *
                         process.lhcInfoTable,
                         process.protonTablesTask 
    )

# MiniAOD output
process.load("Calibration.PPSAlCaRecoProducer.ALCARECOPPSCalMaxTracks_Output_cff")

process.MINIAODoutput = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string(options.miniOutputFilename),
    outputCommands = process.OutALCARECOPPSCalMaxTracks.outputCommands,
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring("p")
    ),
    compressionLevel = cms.untracked.int32(9),
    compressionAlgorithm = cms.untracked.string("LZMA"),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(""),
        dataTier = cms.untracked.string("MINIAOD")
    )
)

# nanoAOD output
process.load('PhysicsTools.NanoAOD.NanoAODEDMEventContent_cff')
process.NANOAODoutput = cms.OutputModule("NanoAODOutputModule",
    fileName = cms.untracked.string(options.nanoOutputFilename),
    outputCommands = process.NANOAODEventContent.outputCommands,
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('p')),
    compressionLevel = cms.untracked.int32(9),
    compressionAlgorithm = cms.untracked.string("LZMA"),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(""),
        dataTier = cms.untracked.string("NANOAOD")
    )
)

process.o = cms.EndPath(
    process.MINIAODoutput +
    process.NANOAODoutput
)