#!/bin/bash

# Get function definitions
source main_functions.sh

## Main code ##
SCRIPT_DIR=`readlink -m "\`dirname \"$0\"\`"`
CLUBB_DIR="$SCRIPT_DIR/../../.."
SIM_OUTPUT_DIR="$SCRIPT_DIR/sim_output"
PLOTS_OUTPUT_DIR="$SCRIPT_DIR/plots_output"

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
echo -e "\n-----Running dycoms2_rf02_do using 2Cat-CldPcp-----"

# 8Cat
turn_on_var_frac "$CLUBB_DIR/input/case_setups/rico_lh_model.in"
turn_on_var_frac "$CLUBB_DIR/input/case_setups/dycoms2_rf02_do_model.in"
set_cluster_strategy_to_1 "$CLUBB_DIR/input/case_setups/rico_lh_model.in"
set_cluster_strategy_to_1 "$CLUBB_DIR/input/case_setups/dycoms2_rf02_do_model.in"
apply_rico_presc_probs "$CLUBB_DIR/input/case_setups/rico_lh_model.in"
apply_dycoms_presc_probs "$CLUBB_DIR/input/case_setups/dycoms2_rf02_do_model.in"
echo -e "\n-----Running rico_lh using 8Cat-----"
echo -e "\n-----Running dycoms2_rf02_do using 8Cat-----"

# 2Cat-Cld
turn_off_presc_probs "$CLUBB_DIR/input/case_setups/rico_lh_model.in"
turn_off_presc_probs "$CLUBB_DIR/input/case_setups/dycoms2_rf02_do_model.in"
echo -e "\n-----Running rico_lh using 2Cat-Cld-----"
echo -e "\n-----Running dycoms2_rf02_do using 2Cat-Cld-----"

# LH-only
turn_off_importance "$CLUBB_DIR/input/case_setups/rico_lh_model.in"
turn_off_importance "$CLUBB_DIR/input/case_setups/dycoms2_rf02_do_model.in"
echo -e "\n-----Running rico_lh without importance sampling-----"
echo -e "\n-----Running dycoms2_rf02_do without importance sampling-----"


# Plots!!
echo -e "\n-----Generating rico_lh RMS plots-----"

echo -e "\n-----Generating dycoms2_rf02_do RMS plots-----"

echo -e "\n-----Generating rico_lh timeseries plots-----"

echo -e "\n-----Generating dycoms2_rf02_do timeseries plots-----"

echo -e "\n-----Generating rico_lh profile plots-----"
plot_profiles_n_dir_all "$CLUBB_DIR" rico_lh "$PLOTS_OUTPUT_DIR/rico_lh_profiles_mean.pdf" --profile_mean "$SIM_OUTPUT_DIR/rico_lh/2Cat-CldPcp" "$SIM_OUTPUT_DIR/rico_lh/8Cat" "$SIM_OUTPUT_DIR/rico_lh/2Cat-Cld"

plot_profiles_n_dir_all "$CLUBB_DIR" rico_lh "$PLOTS_OUTPUT_DIR/rico_lh_profiles_rmse.pdf" --profile_rmse "$SIM_OUTPUT_DIR/rico_lh/2Cat-CldPcp" "$SIM_OUTPUT_DIR/rico_lh/8Cat" "$SIM_OUTPUT_DIR/rico_lh/2Cat-Cld"

echo -e "\n-----Generating dycoms2_rf02_do profile plots-----"
plot_profiles_n_dir_all "$CLUBB_DIR" dycoms2_rf02_do "$PLOTS_OUTPUT_DIR/dycoms2_rf02_do_profiles_mean.pdf" --profile_mean "$SIM_OUTPUT_DIR/dycoms2_rf02_do/2Cat-CldPcp" "$SIM_OUTPUT_DIR/dycoms2_rf02_do/8Cat" "$SIM_OUTPUT_DIR/dycoms2_rf02_do/2Cat-Cld"

plot_profiles_n_dir_all "$CLUBB_DIR" dycoms2_rf02_do "$PLOTS_OUTPUT_DIR/dycoms2_rf02_do_profiles_rmse.pdf" --profile_rmse "$SIM_OUTPUT_DIR/dycoms2_rf02_do/2Cat-CldPcp" "$SIM_OUTPUT_DIR/dycoms2_rf02_do/8Cat" "$SIM_OUTPUT_DIR/dycoms2_rf02_do/2Cat-Cld"

##### Restore model files #####
mv $RICO_LH_MODEL_IN_BACKUP "$CLUBB_DIR/input/case_setups/rico_lh_model.in"
mv $DYCOMS_MODEL_IN_BACKUP "$CLUBB_DIR/input/case_setups/dycoms2_rf02_do_model.in"
