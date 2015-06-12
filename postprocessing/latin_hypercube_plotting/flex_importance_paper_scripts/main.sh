#!/bin/bash

# Get function definitions
source main_functions.sh

## Main code ##
SCRIPT_DIR=`readlink -f "\`dirname \"$0\"\`"`
CLUBB_DIR="$SCRIPT_DIR/../../.."
SIM_OUTPUT_DIR="$SCRIPT_DIR/sim_output"
PLOTS_OUTPUT_DIR="$SCRIPT_DIR/plots_output"

rm -rf "$SIM_OUTPUT_DIR"
rm -rf "$PLOTS_OUTPUT_DIR"
mkdir -p "$SIM_OUTPUT_DIR"
mkdir -p "$PLOTS_OUTPUT_DIR"
echo "Outputting simulation data to `readlink -f $SIM_OUTPUT_DIR`"
echo "Outputting plots to `readlink -f $PLOTS_OUTPUT_DIR`"

# Make backup copies of the model.in files
RICO_LH_MODEL_IN_BACKUP="rico_lh_model.in.backup_$RANDOM"
DYCOMS_MODEL_IN_BACKUP="dycoms_rf02_do_model.in.backup_$RANDOM"

cp "$CLUBB_DIR/input/case_setups/rico_lh_model.in" $RICO_LH_MODEL_IN_BACKUP
cp "$CLUBB_DIR/input/case_setups/dycoms2_rf02_do_model.in" $DYCOMS_MODEL_IN_BACKUP

run_script_in_place $SCRIPT_DIR/modify_model_in_nonint_silhs.py "$CLUBB_DIR/input/case_setups/rico_lh_model.in"
run_script_in_place $SCRIPT_DIR/modify_model_in_nonint_silhs.py "$CLUBB_DIR/input/case_setups/dycoms2_rf02_do_model.in"

run_script_in_place $SCRIPT_DIR/modify_model_in_presc_probs.py "$CLUBB_DIR/input/case_setups/rico_lh_model.in"
run_script_in_place $SCRIPT_DIR/modify_model_in_presc_probs.py "$CLUBB_DIR/input/case_setups/dycoms2_rf02_do_model.in"

apply_rico_presc_probs "$CLUBB_DIR/input/case_setups/rico_lh_model.in"
apply_dycoms_presc_probs "$CLUBB_DIR/input/case_setups/dycoms2_rf02_do_model.in"

echo -e "\n-----Running rico_lh using prescribed probabilities-----"
run_silhs_sp "$CLUBB_DIR" rico_lh "$SIM_OUTPUT_DIR/rico_lh/prescribed"
echo -e "\n-----Running dycoms2_rf02_do using prescribed probabilities-----"
run_silhs_sp "$CLUBB_DIR" dycoms2_rf02_do "$SIM_OUTPUT_DIR/dycoms2_rf02_do/prescribed"

# Cloud weighted sampling
turn_off_presc_probs "$CLUBB_DIR/input/case_setups/rico_lh_model.in"
turn_off_presc_probs "$CLUBB_DIR/input/case_setups/dycoms2_rf02_do_model.in"

echo -e "\n-----Running rico_lh using cloud weighted sampling-----"
run_silhs_sp "$CLUBB_DIR" rico_lh "$SIM_OUTPUT_DIR/rico_lh/cloud_weighted"
echo -e "\n-----Running dycoms2_rf02_do using cloud weighted sampling-----"
run_silhs_sp "$CLUBB_DIR" dycoms2_rf02_do "$SIM_OUTPUT_DIR/dycoms2_rf02_do/cloud_weighted"

# Plots!!
echo -e "\n-----Generating rico_lh RMS plots-----"
mkdir "$PLOTS_OUTPUT_DIR/rico_lh_rms"
plot_rms_n_dir_all "$CLUBB_DIR" rico_lh "$SIM_OUTPUT_DIR/rico_lh/cloud_weighted" "$SIM_OUTPUT_DIR/rico_lh/prescribed" "$PLOTS_OUTPUT_DIR/rico_lh_rms"

echo -e "\n-----Generating dycoms2_rf02_do RMS plots-----"
mkdir "$PLOTS_OUTPUT_DIR/dycoms2_rf02_do_rms"
plot_rms_n_dir_all "$CLUBB_DIR" dycoms2_rf02_do "$SIM_OUTPUT_DIR/dycoms2_rf02_do/cloud_weighted" "$SIM_OUTPUT_DIR/dycoms2_rf02_do/prescribed" "$PLOTS_OUTPUT_DIR/dycoms2_rf02_do_rms"

echo -e "\n-----Generating rico_lh timeseries plots-----"
mkdir "$PLOTS_OUTPUT_DIR/rico_lh_timeseries"
plot_timeseries_two_dir_all "$CLUBB_DIR" rico_lh "$SIM_OUTPUT_DIR/rico_lh/cloud_weighted" "$SIM_OUTPUT_DIR/rico_lh/prescribed" "$PLOTS_OUTPUT_DIR/rico_lh_timeseries"

echo -e "\n-----Generating dycoms2_rf02_do timeseries plots-----"
mkdir "$PLOTS_OUTPUT_DIR/dycoms2_rf02_do_timeseries"
plot_timeseries_two_dir_all "$CLUBB_DIR" dycoms2_rf02_do "$SIM_OUTPUT_DIR/dycoms2_rf02_do/cloud_weighted" "$SIM_OUTPUT_DIR/dycoms2_rf02_do/prescribed" "$PLOTS_OUTPUT_DIR/dycoms2_rf02_do_timeseries"

apply_dycoms_presc_probs "$CLUBB_DIR/input/case_setups/rico_lh_model.in"
apply_rico_presc_probs "$CLUBB_DIR/input/case_setups/dycoms2_rf02_do_model.in"

echo -e "\n-----Running rico_lh using suboptimal prescribed probabilities-----"
run_silhs_sp "$CLUBB_DIR" rico_lh "$SIM_OUTPUT_DIR/rico_lh/prescribed_subopt"
echo -e "\n-----Running dycoms2_rf02_do using suboptimal prescribed probabilities-----"
run_silhs_sp "$CLUBB_DIR" dycoms2_rf02_do "$SIM_OUTPUT_DIR/dycoms2_rf02_do/prescribed_subopt"

echo -e "\n-----Generating rico_lh RMS plots-----"
mkdir "$PLOTS_OUTPUT_DIR/rico_lh_rms_subopt"
plot_rms_n_dir_all "$CLUBB_DIR" rico_lh "$SIM_OUTPUT_DIR/rico_lh/cloud_weighted" "$SIM_OUTPUT_DIR/rico_lh/prescribed" "$SIM_OUTPUT_DIR/rico_lh/prescribed_subopt" "$PLOTS_OUTPUT_DIR/rico_lh_rms_subopt"

echo -e "\n-----Generating dycoms2_rf02_do RMS plots-----"
mkdir "$PLOTS_OUTPUT_DIR/dycoms2_rf02_do_rms_subopt"
plot_rms_n_dir_all "$CLUBB_DIR" dycoms2_rf02_do "$SIM_OUTPUT_DIR/dycoms2_rf02_do/cloud_weighted" "$SIM_OUTPUT_DIR/dycoms2_rf02_do/prescribed" "$SIM_OUTPUT_DIR/dycoms2_rf02_do/prescribed_subopt" "$PLOTS_OUTPUT_DIR/dycoms2_rf02_do_rms_subopt"

##### Restore model files #####
mv $RICO_LH_MODEL_IN_BACKUP "$CLUBB_DIR/input/case_setups/rico_lh_model.in"
mv $DYCOMS_MODEL_IN_BACKUP "$CLUBB_DIR/input/case_setups/dycoms2_rf02_do_model.in"
