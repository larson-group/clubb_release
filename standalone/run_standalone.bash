#!/bin/bash
#######################################################################
# $Id: run_standalone.bash,v 1.29 2008-04-17 00:02:03 dschanen Exp $
#
# Description:
# Script to run the standalone hoc program.  
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
	echo "Usage: run_standalone.bash <MODEL CASE> [PARAMETER FILE] [STATS FILE]"
	exit
else
	MODEL_FILE='../model/'$1'_model.in'

	if [ -z $2 ]; then
		PARAMS_FILE="default_parameters.in"
	else
		PARAMS_FILE=$2
	fi

	if [ -z $3 ]; then
		STATS_FILE="../stats/all_stats.in"
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

cat $PARAMS_FILE $MODEL_FILE $STATS_FILE > 'clubb.in'

#######################################################################
# State which case is being run
#######################################################################
echo "Running" $RUN_CASE

# Run the HOC model
../bin/clubb_standalone 

# Remove the namelists
rm -f 'clubb.in'
