input_file=['file:LOCAL_FILE_PATH'] #local or remote file path containing worker outputs, for the remote path delete the file prefix, ignored when running with automatic script
run = 357440 #insert run number

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '124X_dataRun3_v9', '')

# import FWCore.PythonUtilities.LumiList as LumiList
# process.source.lumisToProcess = LumiList.LumiList(filename = 'INSERT_JSON_MASK').getVLuminosityBlockRange() 