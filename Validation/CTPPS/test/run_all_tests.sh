#!/bin/bash

function RunOne()
{
	local tag="$1"

	cmsRun "${tag}_cfg.py" &> "${tag}.log" &
}


cd "optics_parametrisation_validation"
RunOne "plot_optical_functions"
cd - &> /dev/null


cd "fast_simu_validation"
RunOne "test_acceptance_shape"
RunOne "test_reco_simu_diff_with_det_sm"
RunOne "test_reco_simu_diff_without_det_sm"
RunOne "test_y_vs_x_profile"
RunOne "test_long_extrapolation"
cd - &> /dev/null


cd "fast_simu_with_phys_generator/qgsjet/global"
RunOne "test"
cd - &> /dev/null

cd "fast_simu_with_phys_generator/qgsjet/relative_to_beam"
RunOne "test"
cd - &> /dev/null

#cd "rp_positions/with_particle_gun"
#RunOne "test_with_margin_beg"
#RunOne "test_without_margin_beg"
#RunOne "test_without_margin_end"
#cd - &> /dev/null

cd "rp_positions/with_qgsjet"
RunOne "test_with_margin_beg"
RunOne "test_without_margin_beg"
RunOne "test_without_margin_end"
cd - &> /dev/null
