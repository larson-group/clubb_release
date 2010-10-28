#!/bin/bash
#######################################################################
# $Id$
#
# Description:
# Script to run the CLUBB jacobian program.  
# Tested with GNU Bash v2 & 3.  Might work with Ksh.
#
#######################################################################
# Useful variable on multiprocessor machines with OpenMP capable 
# Fortran, uncomment and set to use.  
# The hoc_standalone program is not parallel, but many LAPACK/BLAS 
# libraries are.
#
# export OMP_NUM_THREADS=2
#######################################################################

if [ -z $1 ]; then
	echo "Usage: "$0" <MODEL CASE> [PARAMETER FILE] [STATS FILE]"
	exit
else
	MODEL_FILE='../input/case_setups/'$1'_model.in'

	if [ -z $2 ]; then
		PARAMS_FILE="../input/tunable_parameters/tunable_parameters.in"
	else
		PARAMS_FILE=$2
	fi

	if [ -z $3 ]; then
		STATS_FILE="../input/stats/standard_stats.in"
	else
		STATS_FILE=$3
	fi
fi
#######################################################################
# Enable G95 Runtime option that sets uninitialized 
# memory to a NaN value
#######################################################################
G95_MEM_INIT="NAN"
export G95_MEM_INIT

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

if [ ! -e "$STATS_FILE" ]; then
	echo "$STATS_FILE does not exist"
	exit 1
fi

cat "../input_misc/jacobian.in" $PARAMS_FILE $MODEL_FILE $STATS_FILE > 'jacobian.in'
#######################################################################
# State which case is being run
#######################################################################
echo "Running" $RUN_CASE

# Run the CLUBB model
../bin/jacobian 

# Remove the namelists
rm -f 'jacobian.in'
