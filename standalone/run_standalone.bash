#!/bin/bash
#######################################################################
# $Id: run_standalone.bash,v 1.28 2007-04-24 22:07:25 dschanen Exp $
#
# Description:
# Script to run the standalone hoc program.  
# Tested with bash v2.  Might work with Ksh.
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
	echo "Usage: run_standalone.bash <CASE>"
	exit
else
	RUN_CASE=$1
fi


#######################################################################
# Check for necessary namelists.  If files exist, then
# copy them over to the general input files.
#######################################################################

STANDALONE_IN='standalone_'$RUN_CASE'.in'
MODEL_IN='../model/'$RUN_CASE'_model.in'
STATS_IN='../stats/'$RUN_CASE'_stats.in'

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
