#!/bin/bash
#######################################################################
# $Id: run_inputfields.bash,v 1.9 2008-05-14 22:13:02 dschanen Exp $
#
# Script to run the hoc inputfields program.  
# Tested with BASH.  Not tested with Korn shell or Bourne(sh) shell.
# Edit to change run
#
#######################################################################
# Useful on compilers which have OpenMP
# OMP_NUM_THREADS=2
#######################################################################
# Check for necessary namelists.  If files exist, then
# copy them over to the general input files.
#######################################################################
if [ -z $1 ]; then
	echo "Usage: run_inputfields.bash <MODEL CASE>"
	exit
else
	RUN_CASE=$1
fi


INPUTFIELDS_IN=$RUN_CASE'_inputfields.in'
MODEL_IN='../model/'$RUN_CASE'_model.in'
STATS_IN='../stats/nobudgets_stats.in'

if [ ! -e "$INPUTFIELDS_IN" ] ; then
	echo $INPUTFIELDS_IN " does not exist"
	exit 1
fi

if [ -e 'inputfields.in' ] ; then
	rm -f 'inputfields.in'
fi

ln -s $INPUTFIELDS_IN 'inputfields.in'
cat $MODEL_IN $STATS_IN > $RUN_CASE'_hoc.in' 

#######################################################################
#
# State which case is being run
echo "Running" $RUN_CASE
# Run HOC
../bin/hoc_inputfields
 
# remove the namelists
rm -f 'inputfields.in'
rm -f $RUN_CASE'_hoc.in'
