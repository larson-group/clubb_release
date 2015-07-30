#!/bin/bash

# Get function definitions
source main_functions.sh

## Main code ##
SCRIPT_DIR=`readlink -m "\`dirname \"$0\"\`"`
CLUBB_DIR="$SCRIPT_DIR/../../.."
SIM_OUTPUT_DIR="$SCRIPT_DIR/sim_output"
PLOTS_OUTPUT_DIR="$SCRIPT_DIR/plots_output"

rm -rf "$SIM_OUTPUT_DIR"
rm -rf "$PLOTS_OUTPUT_DIR"
mkdir -p "$SIM_OUTPUT_DIR"
mkdir -p "$PLOTS_OUTPUT_DIR"
echo "Outputting simulation data to `readlink -m $SIM_OUTPUT_DIR`"
echo "Outputting plots to `readlink -m $PLOTS_OUTPUT_DIR`"

# Make backup copies of the model.in files
RICO_LH_MODEL_IN_BACKUP="rico_lh_model.in.backup_$RANDOM"
DYCOMS_MODEL_IN_BACKUP="dycoms_rf02_do_model.in.backup_$RANDOM"

cp "$CLUBB_DIR/input/case_setups/rico_lh_model.in" $RICO_LH_MODEL_IN_BACKUP
cp "$CLUBB_DIR/input/case_setups/dycoms2_rf02_do_model.in" $DYCOMS_MODEL_IN_BACKUP

run_script_in_place $SCRIPT_DIR/modify_model_in_nonint_silhs.py "$CLUBB_DIR/input/case_setups/rico_lh_model.in"
run_script_in_place $SCRIPT_DIR/modify_model_in_nonint_silhs.py "$CLUBB_DIR/input/case_setups/dycoms2_rf02_do_model.in"

run_script_in_place $SCRIPT_DIR/modify_model_in_presc_probs.py "$CLUBB_DIR/input/case_setups/rico_lh_model.in"
run_script_in_place $SCRIPT_DIR/modify_model_in_presc_probs.py "$CLUBB_DIR/input/case_setups/dycoms2_rf02_do_model.in"

# 2Cat-CldPcp
echo -e "\n-----Running rico_lh using 2Cat-CldPcp-----"
run_silhs_sp "$CLUBB_DIR" rico_lh "$SIM_OUTPUT_DIR/rico_lh/2Cat-CldPcp" 2000
echo -e "\n-----Running dycoms2_rf02_do using 2Cat-CldPcp-----"
run_silhs_sp "$CLUBB_DIR" dycoms2_rf02_do "$SIM_OUTPUT_DIR/dycoms2_rf02_do/2Cat-CldPcp" 2500

# 8Cat
turn_on_var_frac "$CLUBB_DIR/input/case_setups/rico_lh_model.in"
turn_on_var_frac "$CLUBB_DIR/input/case_setups/dycoms2_rf02_do_model.in"
set_cluster_strategy_to_1 "$CLUBB_DIR/input/case_setups/rico_lh_model.in"
set_cluster_strategy_to_1 "$CLUBB_DIR/input/case_setups/dycoms2_rf02_do_model.in"
apply_rico_presc_probs "$CLUBB_DIR/input/case_setups/rico_lh_model.in"
apply_dycoms_presc_probs "$CLUBB_DIR/input/case_setups/dycoms2_rf02_do_model.in"
echo -e "\n-----Running rico_lh using 8Cat-----"
run_silhs_sp "$CLUBB_DIR" rico_lh "$SIM_OUTPUT_DIR/rico_lh/8Cat" 3000
echo -e "\n-----Running dycoms2_rf02_do using 8Cat-----"
run_silhs_sp "$CLUBB_DIR" dycoms2_rf02_do "$SIM_OUTPUT_DIR/dycoms2_rf02_do/8Cat" 3500

# 2Cat-Cld
turn_off_presc_probs "$CLUBB_DIR/input/case_setups/rico_lh_model.in"
turn_off_presc_probs "$CLUBB_DIR/input/case_setups/dycoms2_rf02_do_model.in"
echo -e "\n-----Running rico_lh using 2Cat-Cld-----"
run_silhs_sp "$CLUBB_DIR" rico_lh "$SIM_OUTPUT_DIR/rico_lh/2Cat-Cld" 4000
echo -e "\n-----Running dycoms2_rf02_do using 2Cat-Cld-----"
run_silhs_sp "$CLUBB_DIR" dycoms2_rf02_do "$SIM_OUTPUT_DIR/dycoms2_rf02_do/2Cat-Cld" 4500

# LH-only
turn_off_importance "$CLUBB_DIR/input/case_setups/rico_lh_model.in"
turn_off_importance "$CLUBB_DIR/input/case_setups/dycoms2_rf02_do_model.in"
echo -e "\n-----Running rico_lh without importance sampling-----"
run_silhs_sp "$CLUBB_DIR" rico_lh "$SIM_OUTPUT_DIR/rico_lh/LH-only" 5000
echo -e "\n-----Running dycoms2_rf02_do without importance sampling-----"
run_silhs_sp "$CLUBB_DIR" dycoms2_rf02_do "$SIM_OUTPUT_DIR/dycoms2_rf02_do/LH-only" 5500


# Plots!!
echo -e "\n-----Generating rico_lh RMS plots-----"
plot_rms_n_dir_all "$CLUBB_DIR" rico_lh "$PLOTS_OUTPUT_DIR/rico_lh_rms.pdf" "$SIM_OUTPUT_DIR/rico_lh/LH-only"  "$SIM_OUTPUT_DIR/rico_lh/2Cat-Cld" "$SIM_OUTPUT_DIR/rico_lh/2Cat-CldPcp" "$SIM_OUTPUT_DIR/rico_lh/8Cat"

echo -e "\n-----Generating dycoms2_rf02_do RMS plots-----"
plot_rms_n_dir_all "$CLUBB_DIR" dycoms2_rf02_do "$PLOTS_OUTPUT_DIR/dycoms2_rf02_do_rms.pdf" "$SIM_OUTPUT_DIR/dycoms2_rf02_do/LH-only" "$SIM_OUTPUT_DIR/dycoms2_rf02_do/2Cat-Cld" "$SIM_OUTPUT_DIR/dycoms2_rf02_do/2Cat-CldPcp" "$SIM_OUTPUT_DIR/dycoms2_rf02_do/8Cat"

echo -e "\n-----Generating rico_lh timeseries plots-----"
plot_timeseries_n_dir_all "$CLUBB_DIR" rico_lh "$PLOTS_OUTPUT_DIR/rico_lh_timeseries.pdf" "$SIM_OUTPUT_DIR/rico_lh/2Cat-Cld" "$SIM_OUTPUT_DIR/rico_lh/2Cat-CldPcp"

echo -e "\n-----Generating dycoms2_rf02_do timeseries plots-----"
plot_timeseries_n_dir_all "$CLUBB_DIR" dycoms2_rf02_do "$PLOTS_OUTPUT_DIR/dycoms2_rf02_do_timeseries.pdf" "$SIM_OUTPUT_DIR/dycoms2_rf02_do/2Cat-Cld" "$SIM_OUTPUT_DIR/dycoms2_rf02_do/2Cat-CldPcp"

echo -e "\n-----Generating rico_lh profile plots-----"
plot_profiles_n_dir_all "$CLUBB_DIR" rico_lh "$PLOTS_OUTPUT_DIR/rico_lh_profiles.pdf" "$SIM_OUTPUT_DIR/rico_lh/2Cat-Cld" "$SIM_OUTPUT_DIR/rico_lh/2Cat-CldPcp" "$SIM_OUTPUT_DIR/rico_lh/8Cat"

echo -e "\n-----Generating dycoms2_rf02_do profile plots-----"
plot_profiles_n_dir_all "$CLUBB_DIR" dycoms2_rf02_do "$PLOTS_OUTPUT_DIR/dycoms2_rf02_do_profiles.pdf" "$SIM_OUTPUT_DIR/dycoms2_rf02_do/2Cat-Cld" "$SIM_OUTPUT_DIR/dycoms2_rf02_do/2Cat-CldPcp" "$SIM_OUTPUT_DIR/dycoms2_rf02_do/8Cat"

##### Restore model files #####
mv $RICO_LH_MODEL_IN_BACKUP "$CLUBB_DIR/input/case_setups/rico_lh_model.in"
mv $DYCOMS_MODEL_IN_BACKUP "$CLUBB_DIR/input/case_setups/dycoms2_rf02_do_model.in"
