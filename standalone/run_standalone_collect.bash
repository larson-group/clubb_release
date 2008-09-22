#!/bin/bash
#######################################################################
# $Id: run_standalone.bash,v 1.29 2008-04-17 00:02:03 dschanen Exp $
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
	MODEL_FILE='../model/arm_model.in'
else
	MODEL_FILE='../model/'$2'_model.in'
fi

PARAMS_FILE="default_parameters.in"

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


cat $PARAMS_FILE > 'clubb.in'
#######################################################################
# The following sed commands set l_stats to .false., debug_level to 0,
# and dt = 1.0 sec for profiling purposes.
#######################################################################
cat $MODEL_FILE | sed 's/l_stats\s*=\s*.*/l_stats = \.false\./g' | sed 's/debug_level\s*=\s*.*/debug_level = 0/g' | sed 's/dtmain\s*=\s*.*/dtmain = 1\.0/g' | sed 's/dtclosure\s*=\s*.*/dtclosure = 1\.0/g'  >> 'clubb.in'

echo "Running benchmark"

#######################################################################
# Run collect on the CLUBB model
# Note that -p hi only works on Solaris right now
#######################################################################
collect -p hi -A copy -o $EXPERIMENT_NAME.er ../bin/clubb_standalone 

# Remove the namelists
rm -f 'clubb.in'

echo "Use analyzer " $EXPERIMENT_NAME".er to view the results"
