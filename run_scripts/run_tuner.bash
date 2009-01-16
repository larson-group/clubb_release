#!/bin/bash
#######################################################################
# $Id: run_tuner.bash,v 1.15 2008-04-30 23:11:43 dschanen Exp $
#
# Script to run the tuner.  
# Tested with BASH.  Not tested with Korn shell or Bourne(sh) shell.
# Edit to change run
#
#######################################################################
# Useful on SMP machines with OpenMP
# The tuner is capable of tuning on case per processor if compiled with
# OpenMP compiler flags.
OMP_NUM_THREADS=2
#######################################################################

# Select 'single' for a single case tuning run, or select 'multiple'
# for a multiple case tuning run.
RUN_TYPE='single'
#RUN_TYPE='multiple'

if [ $RUN_TYPE = 'single' ] ; then # Single Case.

   # Note: For a single case tuning run, RUN_CASE should be the name
   #       of the case being tuned for and should match the input
   #       filename "error_RUN_CASE.in".

   # Select a run, comment out the rest.
   # RUN_CASE=arm
   # RUN_CASE=atex
   # RUN_CASE=bomex
   # RUN_CASE=dycoms2_rf01
   # RUN_CASE=dycoms2_rf02_do
   # RUN_CASE=dycoms2_rf02_ds
   # RUN_CASE=dycoms2_rf02_nd
     RUN_CASE=fire
   # RUN_CASE=gabls2 (Needs an error_gabls2.in file before running)
   # RUN_CASE=jun25_altocu (Needs an error_jun25_altocu.in file before running)
   # RUN_CASE=nov11_altocu
   # RUN_CASE=rico (Needs an error_rico.in file before running)
   # RUN_CASE=wangara 

elif [ $RUN_TYPE = 'multiple' ] ; then # Multiple Cases.

   # Note: For a multiple case tuning run, RUN_CASE should match the
   #       input filename "error_RUN_CASE.in".  The names of each of the 
   #       individual cases should be contained in the MODEL_MULT array.

   # Select a RUN_CASE and a MODEL_MULT, comment out the rest.

   # example: all
   # all includes all models with LES data, except for Nov. 11 Altocu
   # (and GABLS2, Jun. 25 Altocu, and RICO).
   # RUN_CASE=all
   # MODEL_MULT=(arm atex bomex dycoms2_rf01 dycoms2_rf02_do\
   #  dycoms2_rf02_ds dycoms2_rf02_nd fire wangara)

   # example: BOMEX and FIRE
   # RUN_CASE=bomex_fire
   # MODEL_MULT=(bomex fire)

   # example: four cases
     RUN_CASE=messner_001
     MODEL_MULT=(arm bomex dycoms2_rf01 dycoms2_rf02_do)

fi

#######################################################################
# Check for necessary namelists.  If files exist, then
# copy them over to the general input files.
#######################################################################
# Path to 'model' directory.
MODEL_DIR='../input/case_setups/'
# Path to 'stats' directory.
STATS_DIR='../input/stats/'
# Today's date in YYYY-MM-DD format.
DATE=`date +%F`

# The error_*.in file.
ERROR_IN=../input/tuner/'error_'$RUN_CASE'.in'

if [ ! -e "$ERROR_IN" ] ; then
	echo $ERROR_IN " does not exist"
	exit 1
fi

if [ $RUN_TYPE = 'single' ] ; then # Single Case.

   # The <CASE>_model.in file.
   MODEL_IN=$MODEL_DIR$RUN_CASE'_model.in'
   if [ ! -e "$MODEL_IN" ] ; then
	   echo $MODEL_IN " does not exist"
	   exit 1
   fi

   # The <STATS>_model.in file.
   # One stats file for tuning, another for the optimal result.
   # They may be the same file, if you wish.
   STATS_TUNE_IN=$STATS_DIR'nobudgets_stats.in'
   if [ ! -e "$STATS_TUNE_IN" ] ; then
	   echo $STATS_TUNE_IN " does not exist"
	   exit 1
   fi
   STATS_OPT_IN=$STATS_DIR'all_stats.in'
   if [ ! -e "$STATS_OPT_IN" ] ; then
	   echo $STATS_OPT_IN " does not exist"
	   exit 1
   fi

   # Concatenate *_model.in and *_stats.in into *_hoc.in
   cat $MODEL_IN $STATS_TUNE_IN > $RUN_CASE'_hoc.in'

elif [ $RUN_TYPE = 'multiple' ] ; then # Multiple Cases.

   for EACH_CASE in "${MODEL_MULT[@]}"; do

           # The <CASE>_model.in file.
           MODEL_IN=$MODEL_DIR$EACH_CASE'_model.in'
           if [ ! -e "$MODEL_IN" ] ; then
	           echo $MODEL_IN " does not exist"
	           exit 1
           fi

           # The <STATS>_model.in file.
           # One stats file for tuning, another for the optimal result
           # They may be the same file, if you wish.
           STATS_TUNE_IN=$STATS_DIR$EACH_CASE'_stats.in'
           if [ ! -e "$STATS_TUNE_IN" ] ; then
	           echo $STATS_TUNE_IN " does not exist"
	           exit 1
           fi
           STATS_OPT_IN=$STATS_DIR$EACH_CASE'_stats.in'
           if [ ! -e "$STATS_OPT_IN" ] ; then
	           echo $STATS_OPT_IN " does not exist"
	           exit 1
           fi

           # Concatenate *_model.in and *_stats.in into *_hoc.in
           cat $MODEL_IN $STATS_TUNE_IN > $EACH_CASE'_hoc.in'

   done

fi

# Copy error_*.in file to error.in
cat $ERROR_IN > 'error.in'

#######################################################################
#
# State which case is being tuned
 echo 'Tuning '$RUN_CASE

# Run tuner, keep a log of STDOUT & STDERR
#
#../bin/clubb_tuner 2>&1 | tee 'tuner_'$RUN_CASE'_'$DATE'.log'
#
# Run tuner, don't keep a log
 ../bin/clubb_tuner 

# Remove the temporary error.in file and the *_hoc.in file(s).
rm -f 'error.in'
rm -f *'_hoc.in'

#######################################################################
# do a run with the optimal constants
echo "Running with optimal constants"

# The newest parameter file should have the optimal set
PARAMS_IN=`ls -t ../input/tunable_parameters* | head -n 1` 

if [ $RUN_TYPE = 'single' ] ; then # Single Case.

   # Concatenate *_model.in and *_stats.in into hoc.in
   # Note:  The <CASE>_hoc.in file is not used because
   #        when the tuner writes standalone_<DATE>.in,
   #        it makes a generic reference to hoc.in
   cat $PARAMS_IN $MODEL_IN $STATS_OPT_IN > 'clubb.in'
   ../bin/clubb_standalone

elif [ $RUN_TYPE = 'multiple' ] ; then # Multiple Cases.

   for EACH_CASE in "${MODEL_MULT[@]}"; do
           MODEL_IN=$MODEL_DIR$EACH_CASE'_model.in'
           STATS_OPT_IN=$STATS_DIR$EACH_CASE'_stats.in'
           # Concatenate *_model.in and *_stats.in into hoc.in
           # Note:  The <CASE>_hoc.in file is not used because
           #        when the tuner writes standalone_<DATE>.in,
           #        it makes a generic reference to hoc.in
           cat $PARAMS_IN $MODEL_IN $STATS_OPT_IN > 'clubb.in'
           ../bin/clubb_standalone
   done

fi

rm -f 'clubb.in'
