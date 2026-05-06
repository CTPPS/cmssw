import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

process = cms.Process('LHCInfoPlotter')

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
    )

options = VarParsing.VarParsing ('analysis')

options.register('outputFileName',
                'LHCInfoPlots.root',
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
options.register('crab',
                False,
                VarParsing.VarParsing.multiplicity.singleton,
                VarParsing.VarParsing.varType.bool,
                "Flag to indicate if running in crab")
options.parseArguments()

#CONFIGURE PROCESS
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

from Configuration.AlCa.GlobalTag import GlobalTag
if options.globalTag != '':
    gt = options.globalTag
else:
    gt = 'auto:run3_data_prompt'

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
print('Using GT:',gt)
process.GlobalTag = GlobalTag(process.GlobalTag, gt)

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
        inputCommands = cms.untracked.vstring(
            'drop *',
            'keep *_TriggerResults_*_*',)
    )

    if options.jsonFileName:
        import FWCore.PythonUtilities.LumiList as LumiList
        jsonFileName = options.jsonFileName
        print('Using JSON file:', jsonFileName)
        process.source.lumisToProcess = LumiList.LumiList(filename = jsonFileName).getVLuminosityBlockRange()
else:
    process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring(),  # will be overridden by crab with the actual input files
        inputCommands = cms.untracked.vstring(
            'drop *',
            'keep *_TriggerResults_*_*',)
    )

process.load('Validation.CTPPS.ctppsLHCInfoPlotter_cfi')
process.ctppsLHCInfoPlotter.outputFile = options.outputFileName
process.ctppsLHCInfoPlotter.useNewLHCInfo = True
    
if options.hltPath:
    print('Using HLT path:', options.hltPath)
    process.load('HLTrigger.HLTfilters.hltHighLevel_cfi')
    process.hltHighLevel.HLTPaths = [options.hltPath]
    process.p = cms.Path(process.hltHighLevel * process.ctppsLHCInfoPlotter)
else:
    process.p = cms.Path(process.ctppsLHCInfoPlotter)