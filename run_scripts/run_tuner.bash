#!/bin/bash
#######################################################################
# $Id: run_tuner.bash,v 1.15 2008-04-30 23:11:43 dschanen Exp $
#
# Script to run the clubb_tuner.  
# Tested with BASH.  Not tested with Korn shell or Bourne(sh) shell.
# Edit to change run
#
#######################################################################
# Useful on SMP machines with OpenMP
# The tuner is capable of tuning one case per processor if compiled with
# OpenMP compiler flags.
OMP_NUM_THREADS=2
#######################################################################

NIGHTLY=false
# Select 'single' for a single case tuning run, or select 'multiple'
# for a multiple case tuning run.
RUN_TYPE='single'
#RUN_TYPE='multiple'

# directory to save initial (pre-tuner) runs
INITIAL_OUTPUT_DIR='../initial_output'
INITIAL_OUTPUT=false

# The code below is borrowed from run_scm.bash to allow command line arguments to this script

# Note that we use `"$@"' to let each command-line parameter expand to a 
# separate word. The quotes around `$@' are essential!
# We need TEMP as the `eval set --' would nuke the return value of getopt.
TEMP=`getopt -o inh --long initial-output,nightly,help \
     -n 'run_tuner.bash' -- "$@"`

if [ $? != 0 ] ; then echo "Terminating..." >&2 ; exit 1 ; fi

# Note the quotes around `$TEMP': they are essential!
eval set -- "$TEMP"

while true ; do
	case "$1" in
	-i|--initial-output)
	# Before tuning, run CLUBB in standalone mode for all cases included in the tuning setup
	# and store the output in folder INITIAL_OUTPUT_DIR
	    INITIAL_OUTPUT=true
	    echo "Executing initial standalone run before tuning"
	    shift;;
	-n|--nightly) 
		NIGHTLY=true
		echo "nightly"
		shift;;
	-h|--help) # Print the help message
		echo -e "Usage: run_tuner.bash [OPTION]"
		echo -e "\t-n, --nightly\t\tRun in a nightly test configuration"
		echo -e "\t-i, --initial-output\tBefore tuning, run CLUBB standalone and store output in separate folder"
		echo -e "\t\t\t\tIf the -n option is given, -i will not work"
		echo -e "\t-h, --help\t\tPrints this help message"

		exit 1 ;;
		--) shift ; break ;;
		*) echo "Something bad happened!" ; exit 1 ;;
	esac
done

# end borrowed code

if [ $RUN_TYPE = 'single' ] ; then # Single Case.

   # Note: For a single case tuning run, RUN_CASE should be the name
   #       of the case being tuned for and should match the input
   #       filename "error_RUN_CASE.in".

   # Select a run, comment out the rest.
   # The selected case will be used in the nightly tuning run
   # RUN_CASE=arm
   # RUN_CASE=arm_97
   # RUN_CASE=atex
   # RUN_CASE=bomex
   # RUN_CASE=cgils_s6
   # RUN_CASE=cgils_s11
   # RUN_CASE=cgils_s12
   # RUN_CASE=dycoms2_rf01
   # RUN_CASE=dycoms2_rf02_do
   # RUN_CASE=dycoms2_rf02_ds
   # RUN_CASE=dycoms2_rf02_nd
    RUN_CASE=fire
   # RUN_CASE=gabls2
   # RUN_CASE=gabls3_night
   # RUN_CASE=jun25_altocu (Needs an error_jun25_altocu.in file before running)
   # RUN_CASE=lba
   # RUN_CASE=nov11_altocu
   # RUN_CASE=rico
   # RUN_CASE=wangara

elif [ $RUN_TYPE = 'multiple' ] ; then # Multiple Cases.

   # Note: For a multiple case tuning run, RUN_CASE should match the
   #       input filename "error_RUN_CASE.in".  The names of each of the 
   #       individual cases should be contained in the MODEL_MULT array.

   # Select a RUN_CASE and a MODEL_MULT, comment out the rest.

   # example: all
   # all includes all models with LES data, except for Nov. 11 Altocu
   # (and GABLS2, Jun. 25 Altocu, and RICO).
   RUN_CASE=fire_atex
   MODEL_MULT=(fire atex)

   # example: BOMEX and FIRE
   # RUN_CASE=bomex_fire
   # MODEL_MULT=(bomex fire)

   # example: four cases
   # RUN_CASE=messner_001
   # MODEL_MULT=(arm bomex dycoms2_rf01 dycoms2_rf02_do)

   # example: rico, cgils_s6/11/12, rf02_nd, and arm
   #RUN_CASE=ticket_756
   #MODEL_MULT=(arm rico dycoms2_rf02_nd)

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

# The tunable_parameters.in file
PARAMS_FILE="../input/tunable_parameters/tunable_parameters.in"
SILHS_PARAMS_FILE="../input/tunable_parameters/silhs_parameters.in"
FLAGS_FILE="../input/tunable_parameters/configurable_model_flags.in"

if [ ! -e "$PARAMS_FILE" ] ; then
	echo $PARAMS_FILE " does not exist"
	exit 1
fi
if [ ! -e "$FLAGS_FILE" ] ; then
	echo $FLAGSS_FILE " does not exist"
	exit 1
fi

# The error_*.in file.
ERROR_IN=../input_misc/tuner/'error_'$RUN_CASE'.in'

if [ ! -e "$ERROR_IN" ] ; then
	echo $ERROR_IN " does not exist"
	exit 1
fi

# The random seed file
RAND_SEED=../input_misc/tuner/rand_seed.dat

if [ ! -e "$RAND_SEED" ] ; then
	echo $RAND_SEED " does not exist"
	exit 1
fi

if [ $RUN_TYPE = 'single' ] ; then # Single Case.

   # The <CASE>_model.in file.
   MODEL_FILE=$MODEL_DIR$RUN_CASE'_model.in'
   if [ ! -e "$MODEL_FILE" ] ; then
	   echo $MODEL_FILE " does not exist"
	   exit 1
   fi

   # The <STATS>_model.in file.
   # One stats file for tuning, another for the optimal result.
   # They may be the same file, if you wish.
   # Note: if STATS_TUNE_IN = tuning_stats.in, the variables listed
   #   there should match the variables being tuned for in 
   #   the error_<RUN_CASE>.in file.
   STATS_TUNE_IN=$STATS_DIR'tuning_stats.in'
   if [ ! -e "$STATS_TUNE_IN" ] ; then
	   echo $STATS_TUNE_IN " does not exist"
	   exit 1
   fi
   STATS_OPT_IN=$STATS_DIR'standard_stats.in'
   if [ ! -e "$STATS_OPT_IN" ] ; then
	   echo $STATS_OPT_IN " does not exist"
	   exit 1
   fi

	# Concatenate *_model.in and *_stats.in into *_hoc.in
	cat $MODEL_FILE $STATS_TUNE_IN $FLAGS_FILE $SILHS_PARAMS_FILE > $RUN_CASE'_hoc.in'
	sed -i -e 's/\!.*//' $RUN_CASE'_hoc.in'

elif [ $RUN_TYPE = 'multiple' ] ; then # Multiple Cases.

   for EACH_CASE in "${MODEL_MULT[@]}"; do

           # The <CASE>_model.in file.
           MODEL_FILE=$MODEL_DIR$EACH_CASE'_model.in'
           if [ ! -e "$MODEL_FILE" ] ; then
	           echo $MODEL_FILE " does not exist"
	           exit 1
           fi

           # The <STATS>_model.in file.
           # One stats file for tuning, another for the optimal result
           # They may be the same file, if you wish.
           STATS_TUNE_IN=$STATS_DIR'tuning_stats.in'
           if [ ! -e "$STATS_TUNE_IN" ] ; then
	           echo $STATS_TUNE_IN " does not exist"
	           exit 1
           fi
           STATS_OPT_IN=$STATS_DIR'standard_stats.in'
           if [ ! -e "$STATS_OPT_IN" ] ; then
	           echo $STATS_OPT_IN " does not exist"
	           exit 1
           fi

        # Concatenate *_model.in and *_stats.in into *_hoc.in
	cat $MODEL_FILE $STATS_TUNE_IN $FLAGS_FILE $SILHS_PARAMS_FILE > $EACH_CASE'_hoc.in'
	sed -i -e 's/\!.*//' $EACH_CASE'_hoc.in'

   done

fi

# Copy error_*.in file and tunable_parameters to error.in
cat $ERROR_IN $PARAMS_FILE > 'error.in'
sed -i -e 's/\!.*//' 'error.in'

# Copy random seed
cp $RAND_SEED .

#######################################################################
#
# # run test cases with initial parameters & save output to separate directory
if ([ $NIGHTLY=false ] && [ $INITIAL_OUTPUT=true ]) ; then

  mkdir $INITIAL_OUTPUT_DIR

  if [ $RUN_TYPE = 'single' ] ; then # Single Case.

    # Concatenate *_model.in and *_stats.in into clubb.in
    cat $STATS_OPT_IN $PARAMS_FILE $SILHS_PARAMS_FILE $MODEL_FILE $FLAGS_FILE | sed -e 's/\!.*//' > 'clubb.in'
     ../bin/clubb_standalone
    mv ../output/$RUN_CASE'_zt'* ../output/$RUN_CASE'_zm'* ../output/$RUN_CASE'_sfc'* $INITIAL_OUTPUT_DIR
    mv ../output/$RUN_CASE'_lh_zt'* ../output/$RUN_CASE'_lh_sfc'* ../output/$RUN_CASE'_rad'* $INITIAL_OUTPUT_DIR
    mv ../output/$RUN_CASE'_setup'* $INITIAL_OUTPUT_DIR

  elif [ $RUN_TYPE = 'multiple' ] ; then # Multiple Cases.

    for EACH_CASE in "${MODEL_MULT[@]}"; do
      MODEL_FILE=$MODEL_DIR$EACH_CASE'_model.in'
      # Concatenate *_model.in and *_stats.in into clubb.in
      cat $STATS_OPT_IN $PARAMS_FILE $MODEL_FILE $FLAGS_FILE | sed -e 's/\!.*//' > 'clubb.in'
        ../bin/clubb_standalone
      mv ../output/$EACH_CASE'_zt'* ../output/$EACH_CASE'_zm'* ../output/$EACH_CASE'_sfc'* $INITIAL_OUTPUT_DIR
      mv ../output/$EACH_CASE'_lh_zt'* ../output/$EACH_CASE'_lh_sfc'* ../output/$EACH_CASE'_rad'* $INITIAL_OUTPUT_DIR
      mv ../output/$EACH_CASE'_setup'* $INITIAL_OUTPUT_DIR
    done

  fi

fi

#######################################################################
#
# State which case is being tuned
 echo 'Tuning '$RUN_CASE

# Run tuner, keep a log of STDOUT & STDERR
#
if ( $NIGHTLY ); then
	logFile='tuner_'$RUN_CASE'_'$DATE'.log'
else
	logFile='tmp.log'
fi

../bin/clubb_tuner 2>&1 | tee $logFile
exitStatus=${PIPESTATUS[0]}

#
# Run tuner, don't keep a log
# Note: Tuner code has option to save results to a file or not
# ../bin/clubb_tuner

# copy the output files to a directory where they will not be erased
# should only work for nightly tests -meyern
if ( $NIGHTLY ); then
	cp ../input/tuning_run_results* $HOME/tuner_output
	cp ../input/error_* $HOME/tuner_output
	cp ../input/tunable_parameters/tunable_parameters_* $HOME/tuner_output
	mv ./*.log $HOME/tuner_output
fi

# Remove the temporary error.in file and the *_hoc.in file(s).
rm -f 'error.in'
rm -f *'_hoc.in'

# Copy the contents of the error_*.in file to the tuner log so they will be saved
TUNER_LOG=`ls -t ../input/tuning_run_results_* | head -n 1`
echo "" >> $TUNER_LOG
echo "********************************" >> $TUNER_LOG
echo "contents of " $ERROR_IN >> $TUNER_LOG
cat $ERROR_IN >> $TUNER_LOG

#######################################################################
# do a run with the optimal constants
echo "Running with the optimal parameter set"

# The newest parameter file should have the optimal set
PARAMS_FILE=`ls -t ../input/tunable_parameters/tunable_parameters_* | head -n 1` 

# For the model flag tuning runs
FLAGS_FILE=`ls -t ../input/tunable_parameters/configurable_model_flags* | head -n 1` 

if [ $RUN_TYPE = 'single' ] ; then # Single Case.

   # Concatenate *_model.in and *_stats.in into clubb.in
   cat $STATS_OPT_IN $PARAMS_FILE $SILHS_PARAMS_FILE $MODEL_FILE $FLAGS_FILE | sed -e 's/\!.*//' > 'clubb.in'
    ../bin/clubb_standalone 2>&1 | tee -a $logFile

elif [ $RUN_TYPE = 'multiple' ] ; then # Multiple Cases.

   for EACH_CASE in "${MODEL_MULT[@]}"; do
		MODEL_FILE=$MODEL_DIR$EACH_CASE'_model.in'
		# Concatenate *_model.in and *_stats.in into clubb.in
        cat $STATS_OPT_IN $PARAMS_FILE $MODEL_FILE $FLAGS_FILE | sed -e 's/\!.*//' > 'clubb.in'
		../bin/clubb_standalone 2>&1 | tee -a $logFile
   done

fi

rm -f 'clubb.in' 'rand_seed.dat'

if [ `grep -c "Fatal error" $logFile` -gt 0 ]; then
	exitStatus=1
fi

if ( ! $NIGHTLY ); then
	rm $logFile
fi

exit $exitStatus
