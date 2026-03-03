## Tracking efficiency aggregation

### Setup

Setup with the following command, tuned to the desired CMSSW version (check that the corresponding rolling-calib branch exists)

```
cmssw_version=15_0_6
cmsrel CMSSW_$cmssw_version
cd CMSSW_$cmssw_version/src
cmsenv
scram-venv
cmsenv
git cms-merge-topic CTPPS:rolling-calib-tracking-efficiency-$cmssw_version # tracking-efficiency
scram b -j10

# Need SSH key uploaded to gitlab for the following command to work
pip3 install "git+https://gitlab.cern.ch/cms-ecal-dpg/ECALELFS/automation-control.git" 
```

### Aggregate results

You will first need to add the InfluxDB credentials to the environment, using
```
export ECAL_INFLUXDB_USER='user'
export ECAL_INFLUXDB_PWD='password'
```

If you don't have credentials, ask A. Bellora or a PPS Automation maintainer.

Results can be aggregated using the command below and providing the required inputs.

```
cd $CMSSW_BASE/src/RecoPPS/RPixEfficiencyTools/test
python3 aggregateEfficiency.py --campaign 2023_track_eff_firstbunches_v4 --begin-run 366727 --end-run 366801 --process-name 'test'
```

Using the `--help` command you can see the description of all the options
- The `--campaign` variable is the automation campaign used in the efficiency production.
- The `--begin-run` and `--end-run` arguments specifies a range of runs to include in the efficiency aggregation (the extremes are included).
- The `process-name` variable is used to tag the log/output files and directories. Use a different one for each aggregation!

The script will allow you to skip some steps (in case you already performed them previously and don't want to repeat them). 

It will first produce the average plane and multiRP efficiency file in the `output_files_${process-name}` directory (called `EfficiencyAnalysis_${process-name}.root`), and then submit a condor job that produces the average pixel unit efficiency (called `ReferenceAnalysis_${process-name}.root`).

The latter output will be available once the HTCondor jobs finish running. The `condor_q` command can be used to monitor the job execution status.