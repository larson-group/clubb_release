#!/bin/bash
#######################################################################
# $Id: run_standalone.bash,v 1.24 2006-12-01 14:13:40 griffinb Exp $
#
# Script to run the standalone hoc program.  
# Tested with bash v2.  Might work with Ksh.
#
#######################################################################
# Useful on multiprocessor machines with OpenMP capable Fortran
# OMP_NUM_THREADS=2
#######################################################################
# Select a run, comment out the rest
if [ -z $1 ]; then
  RUN_CASE=arm 
# RUN_CASE=atex 
# RUN_CASE=bomex 
# RUN_CASE=cobra
# RUN_CASE=dycoms2_rf01 
# RUN_CASE=dycoms2_rf02_do
# RUN_CASE=dycoms2_rf02_ds
# RUN_CASE=dycoms2_rf02_nd 
# RUN_CASE=dycoms2_rf02_so
# RUN_CASE=fire 
# RUN_CASE=jun25_altocu
# RUN_CASE=mpace
# RUN_CASE=nov11_altocu
# RUN_CASE=wangara  

#######################################################################
# Or alternatively specify it as a command line argument in the form
# run_standalone.bash <CASE>.
else
 RUN_CASE=$1
fi


#######################################################################
# Check for necessary namelists.  If files exist, then
# copy them over to the general input files.

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
#
# State which case is being run
 echo "Running" $RUN_CASE
# Run HOC
 ../bin/hoc_standalone

# remove the namelists
 rm -f 'standalone.in'
 rm -f $RUN_CASE'_hoc.in'
