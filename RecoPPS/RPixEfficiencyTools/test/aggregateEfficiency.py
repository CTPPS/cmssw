# This should be executed from the test directory

import os
import subprocess
import sys
import shutil
import logging
import argparse

# Needed to suppress argparse warnings
os.environ["NUMEXPR_MAX_THREADS"] = "1"


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)

common_valid_options = {
    "yes": "yes",
    "no": "no",
    "skip": "skip",
}
def ask_user(prompt, valid_options=common_valid_options, default=None, case_sensitive=False):
    """
    prompt: str
    valid_options: dict mapping accepted inputs -> returned value
                   e.g. {"y": True, "yes": True, "n": False}
    default: value returned if user just presses Enter
    case_sensitive: whether matching is case-sensitive
    """

    # Prepare display string
    display_options = "/".join(sorted(set(valid_options.keys())))

    while True:
        raw = input(f"{prompt} [{display_options}]: ").strip()

        if not raw and default is not None:
            logging.info(f"User selected default: {default}")
            return default

        key = raw if case_sensitive else raw.lower()

        lookup = valid_options if case_sensitive else {
            k.lower(): v for k, v in valid_options.items()
        }

        if key in lookup:
            return lookup[key]

        logging.warning(f"Invalid option. Valid options: {display_options}")


# Parse from command line the basic parameters
parser = argparse.ArgumentParser(description='Aggregate efficiency results.', epilog='Made by A. Bellora - 03/03/2026')
parser.add_argument('--loglevel', default='INFO', help='Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL)')
parser.add_argument('--db-name', default='pps_test_v1', help='InfluxDB database name to use')
parser.add_argument('--campaign', required=True, help='Campaign name of the efficiency production')
parser.add_argument('--process-name', default='test', help='Process name to include in the output file name')
parser.add_argument('--begin-run', type=int, required=True, help='First run to include in the analysis')
parser.add_argument('--end-run', type=int, required=True, help='Last run to include in the analysis')

args = parser.parse_args()
process_name = args.process_name
begin_run = args.begin_run
end_run = args.end_run
dbname = args.db_name
campaign = args.campaign

# Set the logging level based on the command line argument
numeric_level = getattr(logging, args.loglevel.upper(), None)
if not isinstance(numeric_level, int):
    raise ValueError('Invalid log level: {}'.format(args.loglevel))
logging.getLogger().setLevel(numeric_level)

# Check that the credentials are set in the environment variables, otherwise stop
dbusr = os.environ.get('ECAL_INFLUXDB_USER')
dbpwd = os.environ.get('ECAL_INFLUXDB_PWD')
if not dbusr or not dbpwd:
    logging.critical('InfluxDB credentials not set in environment variables. Please set ECAL_INFLUXDB_USER and ECAL_INFLUXDB_PWD.')
    sys.exit(1)

# Need to load after setting user and pwd
from ecalautoctrl import RunCtrl

# Fixed parameters
log_dir = 'logs'
submit_files_dir = 'submit_files'
output_files_dir = 'output_files'
submit_files_dir = f'{submit_files_dir}_{process_name}'
output_files_dir = f'{output_files_dir}_{process_name}'

rctrl = RunCtrl(dbname=dbname, campaign=campaign)

# Status to include a run in the analysis
required_status = {
    "global": "done"
}

# Pull the runs from the DB
runs = rctrl.getRuns(status=required_status, active=False)
run_numbers = [run['run_number'] for run in runs]
run_numbers.sort()
# Include only runs between begin_run and end_run included
run_numbers = [run for run in run_numbers if begin_run <= int(run) <= end_run]
logging.info('Run numbers to analyse:'+','.join(map(str, run_numbers)))

# Get the outputs for those runs
first_worker_outputs = rctrl.getOutput(
    process='re-tracking-efficiency-analysis-worker-with-reco',
    runs=run_numbers[0:10]
)
logging.debug('First worker outputs: {}'.format(first_worker_outputs))

# Create input file list in text file
os.makedirs(submit_files_dir, exist_ok=True)
with open(os.path.join(submit_files_dir, 'input_files_analysis_harvester.txt'), 'w') as f:
    for file in first_worker_outputs:
        f.write('file:{}\n'.format(file))

# Run the analysis harvester
analysis_worker_cmd = [
    'cmsRun',
    '$CMSSW_BASE/src/RecoPPS/RPixEfficiencyTools/python/EfficiencyAnalysisDQMHarvester_cfg.py',
    f'inputFiles={os.path.join(submit_files_dir, "input_files_analysis_harvester.txt")}',
]
analysis_worker_cmd = [os.path.expandvars(c) for c in analysis_worker_cmd]
logging.info('Running command: {}'.format(' '.join(analysis_worker_cmd)))

efficiency_file = f'{output_files_dir}/EfficiencyAnalysis_{process_name}.root'

ans = ask_user('This will run the analysis harvester. Do you want to proceed?')
if ans == 'yes':
    logging.info('Proceeding with analysis harvester')
    result = subprocess.run(
        analysis_worker_cmd,
        env=os.environ.copy(),
        # capture_output=True,
        # text=True    
    )
    # Check the result of the command
    if result.returncode != 0:
        logging.critical('Command failed with return code {}'.format(result.returncode))
        sys.exit(1)
    # Make the output file directory if it does not exist and 
    # move the output file there with a name that includes the process name
    os.makedirs(output_files_dir, exist_ok=True)

    analysis_worker_mv_cmd = [
        'mv',
        'DQM_V0001_R000999999__efficiencyAnalysis__999999__9999999999.root',
        efficiency_file
    ]
    analysis_worker_mv_cmd = [os.path.expandvars(c) for c in analysis_worker_mv_cmd]
    logging.info('Running command: {}'.format(' '.join(analysis_worker_mv_cmd)))
    result = subprocess.run(analysis_worker_mv_cmd)
elif ans == 'no':
    logging.warning('Stopping aggregation here')
    sys.exit(0)
elif ans == 'skip':
    logging.warning('Skipping analysis harvester and moving to the next steps')
    # Check that the efficiency file exists, otherwise stop
    if not os.path.isfile(efficiency_file):
        logging.critical('Efficiency file {} does not exist, cannot proceed'.format(efficiency_file))
        sys.exit(1)
        

from reference_files import REFERENCE_FILE_SETS

ref_files = ask_user('Which reference file set do you want to use?', valid_options=REFERENCE_FILE_SETS)
logging.debug('Reference files: {}'.format(ref_files))
    
with open(os.path.join(submit_files_dir, 'input_files_reference_worker.txt'), 'w') as f:
    for file in ref_files:
        f.write('file:{}\n'.format(file))
        
logging.info('Preparing HTCondor submission files')

# Create a copy of the submission files 
shutil.copytree('ref_efficiency_templates', submit_files_dir, dirs_exist_ok=True)

# Make log dir 
os.makedirs(log_dir, exist_ok=True)

# Substitute the needed variables in the worker submission file
with open(os.path.join(submit_files_dir, 'ref_worker.sub'), 'r') as f:
    content = f.read()
content = content.replace('CMSSW_BASE', os.environ['CMSSW_BASE'])
content = content.replace('EFFICIENCY_FILE', os.path.abspath(efficiency_file))
content = content.replace('OUTPUT_DIR', os.path.abspath(output_files_dir))
content = content.replace('LOGDIR', log_dir)
content = content.replace('LOGFILE', f'log_reference_efficiency_worker_{process_name}')
content = content.replace('INPUT_FILE_LIST', 'input_files_reference_worker.txt')
with open(os.path.join(submit_files_dir, 'ref_worker.sub'), 'w') as f:
    f.write(content)

# Get the list of output files from the efficiency worker (outputReferenceAnalysisDQMWorker_$JOBID.root)
ref_efficiency_worker_outputs = []
for i in range(len(ref_files)):
    ref_worker_output_file = f'outputReferenceAnalysisDQMWorker_{i}.root'
    ref_worker_output_file_path = os.path.abspath(os.path.join(output_files_dir, ref_worker_output_file))
    ref_efficiency_worker_outputs.append(ref_worker_output_file_path)

logging.debug('Reference efficiency worker output files: {}'.format(ref_efficiency_worker_outputs))

ref_harvester_output_file = f'{output_files_dir}/ReferenceEfficiency_{process_name}.root'

# Substitute the needed variables in the harvester submission file
with open(os.path.join(submit_files_dir, 'ref_harvester.sub'), 'r') as f:
    content = f.read()
content = content.replace('CMSSW_BASE', os.environ['CMSSW_BASE'])
content = content.replace('OUTPUT_FILES', os.path.abspath(ref_harvester_output_file))
content = content.replace('LOGDIR', log_dir)
content = content.replace('LOGFILE', f'log_reference_efficiency_harvester_{process_name}')
content = content.replace('INPUT_FILES', ','.join(ref_efficiency_worker_outputs))
with open(os.path.join(submit_files_dir, 'ref_harvester.sub'), 'w') as f:
    f.write(content)

# Use the dag submission file to submit the jobs
ans = ask_user('This will submit the reference efficiency worker and harvester jobs to HTCondor. Do you want to proceed?')
if ans == 'yes':
    logging.info('Submitting jobs to HTCondor')
    result = subprocess.run(
        ['condor_submit_dag','workflow.dag'],
        cwd=submit_files_dir,
    )
    # Check the result of the command
    if result.returncode != 0:
        logging.critical('Command failed with return code {}'.format(result.returncode))
        sys.exit(1)
else:
    logging.warning('Stopping aggregation here')
    sys.exit(0)