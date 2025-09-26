## CalculateDiamondEfficiency.C

This macro generates efficiencies per diamond (and plots with "squares") based on the timing efficiency root file produced by automation and JSON with diamonds' coordinates.

The macro takes following arguments:

- `infile`: path to the efficiency root file (this has to be a path relative to `/eos/cms/store/group/dpg_ctpps/comm_ctpps/pps-automation`)
- `coordsJson`: path to the JSON file with diamonds' coordinates
- `runNumber`: run number
- `fourthDiamondInCyl`: 1 if there are 4 diamonds in the cylindrical RPs, 0 otherwise
- `outputPath`: path to the folder where the output files will be generated (if you want the data to be generated where you are just leave blank)


### Example script that plots efficiency for whole 2024 (all bunches and picked bunches)

```sh
input_prefix=/eos/cms/store/group/dpg_ctpps/comm_ctpps/pps-automation
input_suffix=dev/pps-test-timing-efficiency/test_timing_efficiency_for_2024


for scheme in "all_bunches" "picked_bunches"; do
    for file in "${input_prefix}/${input_suffix}/"*; do 
        file_basename=$(basename $file) # this will be a run number
        if [[ $file_basename =~ ^[0-9]+$ ]]; then
            if (( file_basename >= $1 && file_basename <= $2 )); then
                if [ -d $file/$scheme ] && [ "$(ls -A $file/$scheme)" ]; then
                    (
                        mkdir -p $file_basename/$scheme
                        root_file=$(ls $file/$scheme/outputTimingEfficiencyWorker_*.root | head -n 1)
                        root_file=${root_file#$input_prefix} # input root file without the prefix
                        root -l -b -q "/eos/user/j/jwisniew/CMSSW_15_0_6/src/RecoPPS/TimingEfficiency/tools/CalculateDiamondEfficiency.C(\"$root_file\", \"coords.json\", \"${file_basename}\", \"1\", \"$file_basename/$scheme\")"
                    ) &
                fi
            fi
        fi
    done
    wait

    mkdir summary
    for file in *; do
        if [[ $file =~ ^[0-9]+$ ]]; then
            if (( file >= $1 && file <= $2 )); then
                if [ -d $file/$scheme ] && [ "$(ls -A $file/$scheme)" ]; then
                    (
                        cd $file/$scheme
                        python3 /eos/user/j/jwisniew/CMSSW_15_0_6/src/RecoPPS/TimingEfficiency/tools/merge_jsons.py --input_json efficiency.json --output_json ../../summary/efficiencies_so_far_${scheme}.json
                    )
                fi
            fi
        fi
    done

    mkdir -p summary/$scheme
    if [ -f summary/efficiencies_so_far_${scheme}.json ]; then 
        mv summary/efficiencies_so_far_${scheme}.json summary/$scheme
        (
            cd summary/$scheme
            python3 /eos/user/j/jwisniew/CMSSW_15_0_6/src/RecoPPS/TimingEfficiency/tools/ploter.py --json efficiencies_so_far_${scheme}.json --width 40
        )
    fi
done
```
