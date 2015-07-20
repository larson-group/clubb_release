#!/bin/bash

# This script will run CLUBB-SILHS with several numbers of sample points, and produce
# one output folder in $OUTPUT_DIR for each number of sample points CLUBB is run with.
# For example, $OUTPUT_DIR/out_16 corresponds to output from SILHS run with 16 sample
# points.

CLUBB_DIR="`readlink -m \`dirname $0\``/../../../.."
OUTPUT_DIR="$CLUBB_DIR/rms_plot_output"
CASE_NAME="rico_lh"
STATS_FILE=""
STARTING_SEED=""

# The different numbers of sample points to use, separated by spaces
SAMPLE_POINT_VALUES="2 4 8 16 32 64 128 256"

# Specifices the number of seeds to use to run SILHS with each of the sample point
# values above.
RUN_EACH_WITH_SEEDS=1

#----------------------------------------------------------------------------------------------

# Read command line arguments
while [[ $# -gt 0 ]]
do
    case $1 in
    "--output_dir")
        shift
        OUTPUT_DIR="$1"
        ;;
    "--case_name")
        shift
        CASE_NAME=$1
        ;;
    "--stats_file")
        shift
        STATS_FILE=$1
        ;;
    "--starting_seed")
        shift
        STARTING_SEED=$1
        ;;
    "--run_each_with_seeds")
        shift
        RUN_EACH_WITH_SEEDS=$1
        ;;
    *)
        >&2 echo "Invalid argument: $1"
        exit 1
        ;;
    esac
    shift
done

MODEL_FILE="$CLUBB_DIR/input/case_setups/""$CASE_NAME""_model.in"

# We need to use absolute paths because of a bug in run_scm.bash.
OUTPUT_DIR=`readlink -m "$OUTPUT_DIR"`

# Sanity checks
if [[ ! -f $CLUBB_DIR/run_scripts/run_scm.bash ]]
then
  >&2 echo "The CLUBB run script was not found. Please check the CLUBB_DIR " \
           "variable in the script."
  exit 1
fi

if [[ ! -f $MODEL_FILE ]]
then
  >&2 echo "The specified model file was not found. Please check the MODEL_FILE " \
           "variable in the script."
  exit 1
fi

if [[ -e $OUTPUT_DIR ]]
then
  >&2 echo "The specified output directory already exists. Please check the " \
           "OUTPUT_DIR variable in the script and the corresponding directory " \
           "on the filesystem."
  exit 1
fi

mkdir -p "$OUTPUT_DIR"

# Statistics command line argument string
STATS_CMD_STRING=""
if [[ -n $STATS_FILE ]]
then
    STATS_CMD_STRING="-s $STATS_FILE"
fi

seed=$STARTING_SEED
for num_samples in $SAMPLE_POINT_VALUES
do
    sed 's/lh_num_samples\s*=\s*[0-9]*/lh_num_samples = '"$num_samples"'/g' \
        -i $MODEL_FILE
    for iseed in `seq 1 $RUN_EACH_WITH_SEEDS`
    do
        if [[ -n $seed ]]
        then
            sed 's/lh_seed\s*=\s*[0-9]*/lh_seed = '"$seed"'/g' -i $MODEL_FILE
        fi
        echo "Running with $num_samples samples, iteration $iseed"
        $CLUBB_DIR/run_scripts/run_scm.bash $CASE_NAME $STATS_CMD_STRING \
            -o $OUTPUT_DIR/silhs_"$num_samples"_"$iseed" --netcdf &>/dev/null
        if [[ ! $? -eq 0 ]]
        then
            >&2 echo 'A run failed!'
            exit 2
        fi
        if [[ -n $seed ]]
        then
            seed=$[seed+1]
        fi
    done
done
