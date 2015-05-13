#!/bin/bash

# This script will run CLUBB-SILHS with several numbers of sample points, and produce
# one output folder in $OUTPUT_DIR for each number of sample points CLUBB is run with.
# For example, $OUTPUT_DIR/out_16 corresponds to output from SILHS run with 16 sample
# points.

CLUBB_DIR="../../../.."
OUTPUT_DIR="../../../../rms_plot_output"

# The different numbers of sample points to use, separated by spaces
SAMPLE_POINT_VALUES="2 10 20 30 40 50 60 70 80 90 100 110 120 130 140 150 160 170 180 190 200"

for num_samples in $SAMPLE_POINT_VALUES
do
  sed 's/lh_num_samples\s*=\s*[0-9]*/lh_num_samples = '"$num_samples"'/g' \
    -i $CLUBB_DIR/input/case_setups/rico_lh_model.in
  echo "Running with $num_samples samples"
  $CLUBB_DIR/run_scripts/run_scm.bash rico_lh -o $OUTPUT_DIR/silhs_$num_samples &>/dev/null
  if [ ! $? -eq 0 ]
  then
    echo 'A run failed!'
    exit 1
  fi
done
