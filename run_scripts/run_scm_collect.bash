#!/bin/bash
#######################################################################
# $Id$
#
# Description:
# Script to run an instance of GCSS ARM or other case for profiling 
#
#######################################################################

if [ -z $1 ]; then
	echo "Usage: "$0" <EXPERIMENT NAME> [MODEL CASE]"
	echo "Default model case is GCSS ARM"
	exit
else
	EXPERIMENT_NAME=$1
fi

if [ -z $2 ]; then
	MODEL_FILE='../input/case_setups/arm_model.in'
else
	MODEL_FILE='../input/case_setups/'$2'_model.in'
fi

PARAMS_FILE="../input/tunable_parameters/tunable_parameters.in"
SILHS_PARAMS_FILE="../input/tunable_parameters/silhs_parameters.in"
FLAGS_FILE="../input/tunable_parameters/configurable_model_flags.in"

#######################################################################
# Check for necessary namelists.  If files exist, then
# copy them over to the general input files.
#######################################################################

if [ ! -e "$MODEL_FILE" ]; then
	echo "$MODEL_FILE does not exist"
	exit 1
fi

if [ ! -e "$PARAMS_FILE" ]; then
	echo "$PARAMS_FILE does not exist"
	exit 1
fi


cat $PARAMS_FILE $SILHS_PARAMS_FILE > 'clubb.in'
cat $FLAGS_FILE >> 'clubb.in'
#######################################################################
# The following sed commands set l_stats to .false., debug_level to 0,
# and dt = 1.0 sec for profiling purposes.
#######################################################################
cat $MODEL_FILE | sed 's/l_stats\s*=\s*.*/l_stats = \.false\./g' | sed 's/debug_level\s*=\s*.*/debug_level = 0/g' | sed 's/dt_main\s*=\s*.*/dt_main = 1\.0/g' | sed 's/dt_rad\s*=\s*.*/dt_rad = 1\.0/g'  >> 'clubb.in'

echo "Running benchmark"

#######################################################################
# Run collect on the CLUBB model
# Note that -p hi only works on Solaris right now
#######################################################################
collect -p hi -A copy -o $EXPERIMENT_NAME.er ../bin/clubb_standalone 

# Remove the namelists
rm -f 'clubb.in'

echo "Use analyzer "$EXPERIMENT_NAME".er to view the results"
