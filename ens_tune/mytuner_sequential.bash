#!/bin/bash
# $Id$
################################################################################
#
# Description: 
#   This is a sequential version of Chris Golaz's ensemble tuning script,
#   which was created after the Messner cluster was dismantled.
#   Each ensemable member is run in sequence.
################################################################################

################################################################################
# Global variables
################################################################################

MEMBERMAX=2 # Total number of ensemble members
CLUBB=$HOME/clubb
ENSEMBLE_DIR=`dirname $0` # Directory with this script
CASE=combined_001 # Legacy controls in input_misc/tuner/error_combined_001.in.
CURRENT_PARAMS=""
ARCHIVE=$HOME/archive/ensemble # Directory for the output
EXPERIMENTS=( bomex dycoms2_rf01 ) # Add additional simulations here

################################################################################
# Function to copy files and run the clubb_tuner
################################################################################
function tune ( ) {
	local member=$1

	# Create random seed data
	$CLUBB/install/latest/int2txt /dev/urandom > $CLUBB/ens_tune/rand_seed.dat

	local params_args=()
	if [[ -n "$CURRENT_PARAMS" ]]; then
		params_args=(-params "$CURRENT_PARAMS")
	fi

	python3 "$CLUBB/run_scripts/run_tuner.py" -prepare_only -controls "$CASE" \
        -run_type multiple -model_mult "${EXPERIMENTS[@]}" "${params_args[@]}" \
        -stats_tune tuning_stats -override debug_level=0 -output_dir "$ENSEMBLE_DIR"
	# Execute tuning run
	echo Tuning ensemble member $member on `hostname` at `date +"%R %D"` | tee $CLUBB/ens_tune/tuner.log 
	cd $ENSEMBLE_DIR && echo no | ../install/latest/clubb_tuner &> tune.log

	# Delete GrADS files from tuning run
	rm -f $CLUBB/output/*_zt.??? $CLUBB/output/*_zm.??? $CLUBB/output/*_sfc.???

	TUNED_PARAMS=$(ls -t "$CLUBB"/input/tunable_parameters_*.json | head -1)
	for EXP in "${EXPERIMENTS[@]}"; do
        python3 "$CLUBB/run_scripts/run_scm.py" -params "$TUNED_PARAMS" "$EXP" >> "$ENSEMBLE_DIR/tune.log"
	done

	# Move results data
	mkdir -p $ARCHIVE/$CASE/ens_tune_$member
	mv "$TUNED_PARAMS" "$ARCHIVE/$CASE/ens_tune_$member/"
	mv $ENSEMBLE_DIR/tune.log  $ARCHIVE/$CASE/ens_tune_$member
	# for NetCDF output
	mv $CLUBB/output/*_*.nc  $ARCHIVE/$CASE/ens_tune_$member
}

################################################################################
# Main code
################################################################################

# Check that the archive directory exists
if [ ! -e $ARCHIVE/$CASE ]; then
	# if not, create it
	mkdir -p $ARCHIVE/$CASE
	cp sortresults.bash $ARCHIVE/$CASE
	cp analyze_results.py $ARCHIVE/$CASE
#	cp average.bash $ARCHIVE
fi

# Input conversion occurs immediately before each native tuner launch.

cd $ENSEMBLE_DIR

#Make the les_data folder if it doesn't exsist
mkdir $CLUBB/les_data &> $CLUBB/ens_tune/tuner.log

# Download LES data for experiments
for EXP in "${EXPERIMENTS[@]}"; do
	if [ ! -e "$CLUBB/les_data/${EXP}_coamps_sm.ctl" ]; then
		svn export http://carson.math.uwm.edu/repos/clubb_benchmark_runs/les_runs/${EXP}_coamps_sm.ctl $CLUBB/les_data/${EXP}_coamps_sm.ctl
		svn export http://carson.math.uwm.edu/repos/clubb_benchmark_runs/les_runs/${EXP}_coamps_sm.dat $CLUBB/les_data/${EXP}_coamps_sm.dat
	fi
done

# Loop over the number of ensemble points, incrementing by 1
for (( MEMBER=1; MEMBER <= $MEMBERMAX; MEMBER++ )); do 
	tune $MEMBER
done

