#!/bin/bash
#######################################################################
# $Id: run_inputfields.bash,v 1.9 2008-05-14 22:13:02 dschanen Exp $
#
# Script to run the CLUBB inputfields program, which runs CLUBB with
# the prognostic fields specified in a GrADS file.
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


INPUTFIELDS_IN='../input_misc/inputfields/'$RUN_CASE'_inputfields.in'
MODEL_IN='../input/case_setups/'$RUN_CASE'_model.in'
STATS_IN='../input/stats/nobudgets_stats.in'
PARAMS_IN='../input/tunable_parameters.in'

if [ ! -e "$INPUTFIELDS_IN" ] ; then
	echo $INPUTFIELDS_IN " does not exist"
	exit 1
fi

if [ -e 'inputfields.in' ] ; then
	rm -f 'inputfields.in'
fi

cat $INPUTFIELDS_IN $PARAMS_IN > 'inputfields.in'
cat $MODEL_IN $STATS_IN > $RUN_CASE'_hoc.in' 

#######################################################################
#
# State which case is being run
echo "Running" $RUN_CASE
# Run CLUBB
../bin/clubb_inputfields
 
# remove the namelists
rm -f 'inputfields.in'
rm -f $RUN_CASE'_hoc.in'
