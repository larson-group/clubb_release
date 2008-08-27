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
	echo "Usage: run_standalone.bash <MODEL CASE> [STATS FILE]"
	exit
else
	RUN_CASE=$1

	if [ -z $2 ]; then
		STATS_FILE="all_stats.in"
	else
		STATS_FILE=$2
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

STANDALONE_IN='standalone_'$RUN_CASE'.in'
MODEL_IN='../model/'$RUN_CASE'_model.in'
STATS_IN='../stats/'$STATS_FILE

if [ ! -e "$STANDALONE_IN" ] ; then
	echo $STANDALONE_IN " does not exist"
	exit 1
fi

if [ -e 'standalone.in' ] ; then
	rm -f 'standalone.in'
fi

ln -s $STANDALONE_IN 'standalone.in'
cat $MODEL_IN $STATS_IN > $RUN_CASE'_hoc.in'

#######################################################################
# State which case is being run
#######################################################################
echo "Running" $RUN_CASE

# Run the HOC model
../bin/hoc_standalone 

# Remove the namelists
rm -f 'standalone.in'
rm -f $RUN_CASE'_hoc.in'
