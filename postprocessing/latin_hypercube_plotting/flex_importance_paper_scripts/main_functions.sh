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
# 3: Plot output directory
# +: Simulation directories
plot_rms_n_dir_all()
{
    clubb_dir=$1; shift
    case_name=$1; shift
    plot_dir=$1; shift
    script="$clubb_dir"/postprocessing/latin_hypercube_plotting/variance_analysis_scripts/rms_vs_sample_points/silhs_rms_plot_mult_sim.py
    mkdir -p "$plot_dir"
    if [[ "$case_name" == "rico_lh" ]]
    then
        time1=0
        time2=4320
    else
        time1=0
        time2=360
    fi
    # Total tendency
    echo "Total tendency"
    $script --plot_title "$TOTAL_STR" --time1 $time1 --time2 $time2 --case_name "$case_name" --clubb_var_str "rrm_mc" --silhs_var_str "lh_rrm_mc" --output_file "$plot_dir"/rrm_mc_rms.svg "$@"
    # Autoconversion tendency
    echo "Autoconversion tendency"
    $script --plot_title "$AUTOCONV_STR" --time1 $time1 --time2 $time2 --case_name "$case_name" --clubb_var_str "rrm_auto" --silhs_var_str "lh_rrm_auto" --output_file "$plot_dir"/rrm_auto_rms.svg "$@"
    # Accretion tendency
    echo "Accretion tendency"
    $script --plot_title "$ACCR_STR" --time1 $time1 --time2 $time2 --case_name "$case_name" --clubb_var_str "rrm_accr" --silhs_var_str "lh_rrm_accr" --output_file "$plot_dir"/rrm_accr_rms.svg "$@"
    # Evaporation tendency
    echo "Evaporation tendency"
    $script --plot_title "$EVAP_STR" --time1 $time1 --time2 $time2 --case_name "$case_name" --clubb_var_str "rrm_cond" --silhs_var_str "lh_rrm_evap" --output_file "$plot_dir"/rrm_evap_rms.svg "$@"
}

# Parameters
# 1: CLUBB directory
# 2: Case name
# 3: Simulation directory 1
# 4: Simulation directory 2
# 5: Output directory
plot_timeseries_two_dir_all()
{
    clubb_dir="$1"
    case_name="$2"
    sim_dir_1="$3"
    sim_dir_2="$4"
    out_dir="$5"

    script="$clubb_dir/postprocessing/latin_hypercube_plotting/variance_analysis_scripts/timeseries_k_lh_start/silhs_kplot.py"
    mkdir -p "$out_dir"

    tmp_dir="tmp_$RANDOM"
    mkdir $tmp_dir
    ln -s "`readlink -f \"$sim_dir_1/silhs_20\"`" "$tmp_dir/`basename \"$sim_dir_1\"`_20"
    ln -s "`readlink -f \"$sim_dir_2/silhs_20\"`" "$tmp_dir/`basename \"$sim_dir_2\"`_20"

    if [[ "$case_name" == "rico_lh" ]]
    then
        time1=0
        time2=4320
    else
        time1=0
        time2=360
    fi

    # Total tendency
    echo "Total tendency"
    $script --title "$TOTAL_STR" --time1 $time1 --time2 $time2 --case_name "$case_name" --clubb_var "rrm_mc" --silhs_var "lh_rrm_mc" --ylabel "$TOTAL_EQ" --output_file "$out_dir/rrm_mc_timeseries.svg" "$tmp_dir/`basename \"$sim_dir_1\"`_20" "$tmp_dir/`basename \"$sim_dir_2\"`_20"
    echo "Autoconversion tendency"
    # Autoconversion tendency
    $script --title "$AUTOCONV_STR" --time1 $time1 --time2 $time2 --case_name "$case_name" --clubb_var "rrm_auto" --silhs_var "lh_rrm_auto" --ylabel "$AUTOCONV_EQ" --output_file "$out_dir/rrm_auto_timeseries.svg" "$tmp_dir/`basename \"$sim_dir_1\"`_20" "$tmp_dir/`basename \"$sim_dir_2\"`_20"
    # Accretion tendency
    echo "Accretion tendency"
    $script --title "$ACCR_STR" --time1 $time1 --time2 $time2 --case_name "$case_name" --clubb_var "rrm_accr" --silhs_var "lh_rrm_accr" --ylabel "$ACCR_EQ" --output_file "$out_dir/rrm_accr_timeseries.svg" "$tmp_dir/`basename \"$sim_dir_1\"`_20" "$tmp_dir/`basename \"$sim_dir_2\"`_20"
    # Evaporation tendency
    echo "Evaporation tendency"
    $script --title "$EVAP_STR" --time1 $time1 --time2 $time2 --case_name "$case_name" --clubb_var "rrm_cond" --silhs_var "lh_rrm_evap" --ylabel "$EVAP_EQ" --output_file "$out_dir/rrm_evap_timeseries.svg" "$tmp_dir/`basename \"$sim_dir_1\"`_20" "$tmp_dir/`basename \"$sim_dir_2\"`_20"

    # Clean-up
    rm -r $tmp_dir
}

# Parameters
# 1: Case file to write the result to
apply_rico_presc_probs()
{
    sed 's/^eight_cluster_presc_probs%cloud_precip_comp1\s*=.*$/eight_cluster_presc_probs%cloud_precip_comp1 = 0.3/g' -i "$1"
    sed 's/^eight_cluster_presc_probs%cloud_precip_comp2\s*=.*$/eight_cluster_presc_probs%cloud_precip_comp2 = 0.04/g' -i "$1"
    sed 's/^eight_cluster_presc_probs%nocloud_precip_comp1\s*=.*$/eight_cluster_presc_probs%nocloud_precip_comp1 = 0.2/g' -i "$1"
    sed 's/^eight_cluster_presc_probs%nocloud_precip_comp2\s*=.*$/eight_cluster_presc_probs%nocloud_precip_comp2 = 0.04/g' -i "$1"
    sed 's/^eight_cluster_presc_probs%cloud_noprecip_comp1\s*=.*$/eight_cluster_presc_probs%cloud_noprecip_comp1 = 0.3/g' -i "$1"
    sed 's/^eight_cluster_presc_probs%cloud_noprecip_comp2\s*=.*$/eight_cluster_presc_probs%cloud_noprecip_comp2 = 0.04/g' -i "$1"
    sed 's/^eight_cluster_presc_probs%nocloud_noprecip_comp1\s*=.*$/eight_cluster_presc_probs%nocloud_noprecip_comp1 = 0.04/g' -i "$1"
    sed 's/^eight_cluster_presc_probs%nocloud_noprecip_comp2\s*=.*$/eight_cluster_presc_probs%nocloud_noprecip_comp2 = 0.04/g' -i "$1"
}

# Parameters
# 1: Case file to write the result to
apply_dycoms_presc_probs()
{
    sed 's/^eight_cluster_presc_probs%cloud_precip_comp1\s*=.*$/eight_cluster_presc_probs%cloud_precip_comp1 = 0.49/g' -i "$1"
    sed 's/^eight_cluster_presc_probs%cloud_precip_comp2\s*=.*$/eight_cluster_presc_probs%cloud_precip_comp2 = 0.39/g' -i "$1"
    sed 's/^eight_cluster_presc_probs%nocloud_precip_comp1\s*=.*$/eight_cluster_presc_probs%nocloud_precip_comp1 = 0.02/g' -i "$1"
    sed 's/^eight_cluster_presc_probs%nocloud_precip_comp2\s*=.*$/eight_cluster_presc_probs%nocloud_precip_comp2 = 0.02/g' -i "$1"
    sed 's/^eight_cluster_presc_probs%cloud_noprecip_comp1\s*=.*$/eight_cluster_presc_probs%cloud_noprecip_comp1 = 0.02/g' -i "$1"
    sed 's/^eight_cluster_presc_probs%cloud_noprecip_comp2\s*=.*$/eight_cluster_presc_probs%cloud_noprecip_comp2 = 0.02/g' -i "$1"
    sed 's/^eight_cluster_presc_probs%nocloud_noprecip_comp1\s*=.*$/eight_cluster_presc_probs%nocloud_noprecip_comp1 = 0.02/g' -i "$1"
    sed 's/^eight_cluster_presc_probs%nocloud_noprecip_comp2\s*=.*$/eight_cluster_presc_probs%nocloud_noprecip_comp2 = 0.02/g' -i "$1"
}

# Parameters
# 1: script
# 2: file
run_script_in_place()
{
    tmpfile="tmpfile_$RANDOM"
    "$1" "$2" "$tmpfile"
    mv -f "$tmpfile" "$2"
}

# Parameters
# 1: Model file
turn_off_presc_probs()
{
    sed 's/^l_lh_clustered_sampling\s*=.*$/l_lh_clustered_sampling = .false./g' -i "$1"
}
