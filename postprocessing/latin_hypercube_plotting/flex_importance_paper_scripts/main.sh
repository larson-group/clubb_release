#!/bin/bash

TOTAL_EQ='$\left(\frac{\partial r_r}{\partial t}\right)$'
AUTOCONV_EQ='$\left(\frac{\partial r_r}{\partial t}\right)_\mathrm{auto}$'
ACCR_EQ='$\left(\frac{\partial r_r}{\partial t}\right)_\mathrm{accr}$'
EVAP_EQ='$\left(\frac{\partial r_r}{\partial t}\right)_\mathrm{evap}$'

TOTAL_STR="Total microphysical rain tendency $TOTAL_EQ"
AUTOCONV_STR="Rain autoconversion tendency $AUTOCONV_EQ"
ACCR_STR="Rain accretion tendency $ACCR_EQ"
EVAP_STR="Rain evaporation tendency $EVAP_EQ"

# Parameters
# 1: CLUBB directory
# 2: Case name
# 3: Output directory
run_silhs_sp()
{
    "$1"/postprocessing/latin_hypercube_plotting/variance_analysis_scripts/rms_vs_sample_points/silhs_varying_sp_output.sh --stats_file "$1"/postprocessing/latin_hypercube_plotting/flex_importance_paper_scripts/stats.in --case_name "$2" --output_dir "$3"
}

# Parameters
# 1: CLUBB directory
# 2: Case name
# 3: Simulation directory 1
# 4: Simulation directory 2
# 5: Output directory
plot_rms_two_dir_all()
{
    script="$1"/postprocessing/latin_hypercube_plotting/variance_analysis_scripts/rms_vs_sample_points/silhs_rms_plot_mult_sim.py
    if [[ "$2" == "rico_lh" ]]
    then
        time1=0
        time2=4320
    else
        time1=0
        time2=360
    fi
    # Total tendency
    $script --plot_title "$TOTAL_STR" --time1 $time1 --time2 $time2 --case_name "$2" --clubb_var_str "rrm_mc" --silhs_var_str "lh_rrm_mc" --output_file "$5"/rrm_mc_rms.svg "$3" "$4"
    # Autoconversion tendency
    $script --plot_title "$AUTOCONV_STR" --time1 $time1 --time2 $time2 --case_name "$2" --clubb_var_str "rrm_auto" --silhs_var_str "lh_rrm_auto" --output_file "$5"/rrm_auto_rms.svg "$3" "$4"
    # Accretion tendency
    $script --plot_title "$ACCR_STR" --time1 $time1 --time2 $time2 --case_name "$2" --clubb_var_str "rrm_accr" --silhs_var_str "lh_rrm_accr" --output_file "$5"/rrm_accr_rms.svg "$3" "$4"
    # Evaporation tendency
    $script --plot_title "$EVAP_STR" --time1 $time1 --time2 $time2 --case_name "$2" --clubb_var_str "rrm_cond" --silhs_var_str "lh_rrm_evap" --output_file "$5"/rrm_evap_rms.svg "$3" "$4"
}
