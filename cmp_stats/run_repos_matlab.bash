#!/bin/bash
#######################################################################
# $Id: run_repos_matlab.bash,v 1.5 2008-01-03 20:48:24 nielsenb Exp $
#
# Script to run the standalone hoc program for all models, with 
# repository (non-static) constants.
# Tested with bash v2.  Might work with Ksh.
#
#######################################################################
# Useful on multiprocessor machines with OpenMP capable Fortran
#export OMP_NUM_THREADS=2
#######################################################################

EXIT_CODE=( [0]=0 [1]=0 [2]=0 [3]=0 [4]=0 [5]=0 [6]=0 [7]=0 [8]=0 [9]=0 [10]=0 [11]=0 [12]=0 [13]=0 [14]=0 [15]=0 )

RUN_CASE=( \
 arm atex bomex cobra dycoms2_rf01 dycoms2_rf02_do dycoms2_rf02_ds\
 dycoms2_rf02_nd dycoms2_rf02_so fire gabls2 jun25_altocu mpace_b nov11_altocu\
 rico wangara )

#start nielsenb's changes
#eliminate the previous HOC results
#this prevents spurious profile generation resulting from
#previous profiles not getting overwritten
rm -f HOC_previous/*
#end nielsenb's changes

mv HOC_current/*.ctl HOC_previous/
mv HOC_current/*.dat HOC_previous/

# This will loop over all runs in sequence 
for (( x=0; x < "${#RUN_CASE[@]}"; x++ )); do
#######################################################################
# Check for necessary namelists.  If files exist, then
# copy them over to the general input files.

 STANDALONE_IN='../standalone/standalone_'"${RUN_CASE[$x]}"'.in'
 MODEL_IN='../model/'"${RUN_CASE[$x]}"'_model.in'
 STATS_IN="${RUN_CASE[$x]}"'_stats.in'

 if [ ! -e "$STANDALONE_IN" ] ; then
	echo $STANDALONE_IN " does not exist"
	exit 1
 fi

 if [ -e 'standalone.in' ] ; then
	rm -f 'standalone.in'
 fi

 ln -s $STANDALONE_IN 'standalone.in'
 cat $MODEL_IN $STATS_IN > "${RUN_CASE[$x]}"'_hoc.in'

#######################################################################
#
# State which case is being run
 echo "Running ""${RUN_CASE[$x]}"
# Run HOC 
 RESULT=`../bin/hoc_standalone 2>&1 |grep 'normal'`
#echo $RESULT
 if [ -z "$RESULT" ]; then
	EXIT_CODE[$x]=-1
 fi

# remove the namelists
 rm -f 'standalone.in'
 rm -f "${RUN_CASE[$x]}"'_hoc.in'

done

# Print the results & copy files
for (( x=0; x < "${#RUN_CASE[@]}"; x++ )); do
  if [ "${EXIT_CODE[$x]}" != 0 ]; then
	echo "${RUN_CASE[$x]}"' failure'
	rm "${RUN_CASE[$x]}"_zt.ctl
	rm "${RUN_CASE[$x]}"_zt.dat
	rm "${RUN_CASE[$x]}"_zm.ctl
	rm "${RUN_CASE[$x]}"_zm.dat
  else
	mv "${RUN_CASE[$x]}"_zt.ctl HOC_current/
	mv "${RUN_CASE[$x]}"_zt.dat HOC_current/
	mv "${RUN_CASE[$x]}"_zm.ctl HOC_current/
	mv "${RUN_CASE[$x]}"_zm.dat HOC_current/
  fi
done
