#!/bin/bash

# Parameters
# 1: CLUBB directory
# 2: Case name
# 3: Output directory
# 4: Starting seed
run_silhs_sp()
{
    "$1"/postprocessing/latin_hypercube_plotting/variance_analysis_scripts/rms_vs_sample_points/silhs_varying_sp_output.sh --stats_file "$1"/postprocessing/latin_hypercube_plotting/flex_importance_paper_scripts/stats.in --case_name "$2" --output_dir "$3" --starting_seed "$4" --run_each_with_seeds 12
}

# Parameters
# 1: CLUBB dir
# 2: Case name
# 3: Plot output file
# +: Input directories
plot_rms_n_dir_all()
{
    clubb_dir="$1"
    shift
    case_name="$1"
    shift
    plot_output_file="$1"
    shift

    if [[ $case_name == "rico_lh" ]]; then
        time1=0
        time2=864
    else
        time1=0
        time2=360
    fi

    plot_something_n_dir_all "$clubb_dir" "$case_name" $time1 $time2 "$plot_output_file" --rms "$@"
}

# Parameters
# 1: CLUBB dir
# 2: Case name
# 3: Plot output file
# +: Input directories
plot_timeseries_n_dir_all()
{
    clubb_dir="$1"
    shift
    case_name="$1"
    shift
    plot_output_file="$1"
    shift

    if [[ $case_name == "rico_lh" ]]; then
        time1=664
        time2=864
    else
        time1=160
        time2=360
    fi

    plot_something_n_dir_all "$clubb_dir" "$case_name" $time1 $time2 "$plot_output_file" --timeseries "$@"
}

# Parameters
# 1: CLUBB dir
# 2: Case name
# 3: time1
# 4: time2
# 5: Plot output file
# 6: Plot type argument
# +: Input directories
plot_something_n_dir_all()
{
    clubb_dir="$1"
    shift
    case_name="$1"
    shift
    time1="$1"
    shift
    time2="$1"
    shift
    plot_file="$1"
    shift
    plot_type_arg="$1"
    shift

    script="$clubb_dir/postprocessing/latin_hypercube_plotting/variance_analysis_scripts/rms_vs_sample_points/silhs_rms_timeseries_profiles_4pan_mult_sim.py"

    # Call the script!
    "$script" --case_name "$case_name" --time1 "$time1" --time2 "$time2" --output_file "$plot_file" "$plot_type_arg" "$@"
}

# Parameters
# 1: CLUBB dir
# 2: Case name
# 3: Plot output file
# +: Input directories
plot_profiles_n_dir_all()
{
    clubb_dir="$1"
    shift
    case_name="$1"
    shift
    plot_output_file="$1"
    shift

    if [[ $case_name == "rico_lh" ]]; then
        time1=0
        time2=864
    else
        time1=0
        time2=360
    fi

    plot_something_n_dir_all "$clubb_dir" "$case_name" $time1 $time2 "$plot_output_file" --profiles "$@"
}


# Parameters
# 1: Case file to write the result to
apply_rico_presc_probs()
{
    sed 's/^eight_cluster_presc_probs%cloud_precip_comp1\s*=.*$/eight_cluster_presc_probs%cloud_precip_comp1 = 0.35/g' -i "$1"
    sed 's/^eight_cluster_presc_probs%cloud_precip_comp2\s*=.*$/eight_cluster_presc_probs%cloud_precip_comp2 = 0.14/g' -i "$1"
    sed 's/^eight_cluster_presc_probs%nocloud_precip_comp1\s*=.*$/eight_cluster_presc_probs%nocloud_precip_comp1 = 0.24/g' -i "$1"
    sed 's/^eight_cluster_presc_probs%nocloud_precip_comp2\s*=.*$/eight_cluster_presc_probs%nocloud_precip_comp2 = 0.06/g' -i "$1"
    sed 's/^eight_cluster_presc_probs%cloud_noprecip_comp1\s*=.*$/eight_cluster_presc_probs%cloud_noprecip_comp1 = 0.14/g' -i "$1"
    sed 's/^eight_cluster_presc_probs%cloud_noprecip_comp2\s*=.*$/eight_cluster_presc_probs%cloud_noprecip_comp2 = 0.07/g' -i "$1"
    sed 's/^eight_cluster_presc_probs%nocloud_noprecip_comp1\s*=.*$/eight_cluster_presc_probs%nocloud_noprecip_comp1 = 0.0/g' -i "$1"
    sed 's/^eight_cluster_presc_probs%nocloud_noprecip_comp2\s*=.*$/eight_cluster_presc_probs%nocloud_noprecip_comp2 = 0.0/g' -i "$1"
}

# Parameters
# 1: Case file to write the result to
apply_dycoms_presc_probs()
{
    apply_rico_presc_probs "$1"
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

# Parameters
# 1: Model file
turn_on_presc_probs()
{
    sed 's/^l_lh_clustered_sampling\s*=.*$/l_lh_clustered_sampling = .true./g' -i "$1"
}

# Parameters
# 1: Model file
turn_off_importance()
{
    sed 's/^l_lh_importance_sampling\s*=.*$/l_lh_importance_sampling = .false./g' -i "$1"
}


# Parameters
# 1: Model file
set_cluster_strategy_to_1()
{
    sed 's/^cluster_allocation_strategy\s*=.*$/cluster_allocation_strategy = 1/g' -i "$1"
}

# Parameters
# 1: Model file
turn_off_weight_limiter()
{
    sed 's/^l_lh_limit_weights\s*=.*$/l_lh_limit_weights = .false.,/g' -i "$1"
}

# Parameters
# 1: Model file
turn_on_var_frac()
{
    sed 's/^l_lh_var_frac\s*=.*$/l_lh_var_frac = .true.,/g' -i "$1"
}
