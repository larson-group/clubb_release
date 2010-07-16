#!/bin/bash
##########################################################################
# $Id:$
#
# Description:
# Script that runs a predefined set of test cases for creating an overall
#   timing profile.
#
# This script will run BOMEX, COBRA, nov11_altocu, RICO, and
#   cloud_feedback_s12 5 times (by default) each in order to ensure that
#   consistent testing data is created. This is to give a good sampling of
#   different types of cases that CLUBB can handle:
#     COBRA uses neither radiation nor microphysics
#     BOMEX uses simplified_bomex radiation
#     nov11_altocu uses COAMPS microphysics and simplified radiation
#     RICO uses KK microphysics
#     cloud_feedback_s12 uses Morrison microphysics and bugsrad radiation
#
# Additional cases can be run by adding these cases to the RUN_CASE array
#   declared at the top of the script. If you wish to extend the time_final
#   for better benchmarks, specify the time_final in the list where the
#   variable TIME_FINAL is declared. The typical recommendation depends on
#   case, but a good starting point is generally the standard value
#   multiplied by 10 or 20.
#
# NOTE: When running this script, it is highly recommended that the
#   computer is doing nothing else--having any other programs running,
#   even an instance of GVim, will likely increase the inconsistency
#   of the results. Ideally, try to ensure that the only program running
#   on that computer is your SSH session/terminal window that started
#   the script.
#
##########################################################################

# List of the cases that will run in the profile testing suite
RUN_CASE=( bomex cobra nov11_altocu rico cloud_feedback_s12 )

# Figure out the directory where the script is located
scriptPath=`dirname $0`

# Store the current working directory so it can be restored
restoreDir=`pwd`

# Change directories to the script's directory
cd $scriptPath

if [ -z $1 ]; then
	echo "Usage: "$0" <OUTPUT DIRECTORY> [NUMBER OF TESTS]"
	echo "The output directory will contain the output from"
	echo "all of the profile tests."
	echo ""
	echo "The default number of tests is 5."
	exit
else
	OUTPUT_DIR=$1
fi

if [ -z $2 ]; then
	NUM_TESTS=5
else
	NUM_TESTS=$2
fi

echo "NOTE: It is highly recommended that no other programs are running"
echo "on this computer. If other programs are running, the results may"
echo "be inconsistent and difficult to interpret, and may be misleading."

# Create the output directory if it does not exist
if [ ! -d "./$OUTPUT_DIR" ]; then
	mkdir $OUTPUT_DIR
else
	echo "ERROR: Output directory already exists. Exiting."
	exit 1
fi

echo "Running test suite--each case will run "$NUM_TESTS" times."

##########################################################################
# Loop through the cases in the suite, setting up the parameters so testing
# runs can be completed
##########################################################################
for (( x=0; x < "${#RUN_CASE[@]}"; x++ )); do
	echo "Performing tests for case: "${RUN_CASE[$x]}"."
	echo "---------------------------------------------"

	# Ensure necessary namelists exist and copy them to the general
	# input files.

	MODEL_FILE="../input/case_setups/"${RUN_CASE[$x]}"_model.in"
	PARAMS_FILE="../input/tunable_parameters.in"

	if [ ! -e "$MODEL_FILE" ]; then
		echo "ERROR: $MODEL_FILE does not exist."
		exit 1
	fi

	if [ ! -e "$PARAMS_FILE" ]; then
		echo "ERROR: $PARAMS_FILE does not exist."
		exit 1
	fi

	cat $PARAMS_FILE > "clubb.in"

	# Determine proper end-time parameter for various cases
	case "${RUN_CASE[$x]}" in
		"bomex" ) FINAL_TIME="6480000\.0";;
		"cobra" ) FINAL_TIME="2304000\.0";;
		"nov11_altocu" ) FINAL_TIME="3780000\.0";;
		"rico" ) FINAL_TIME="5184000\.0";;
		"cloud_feedback_s12" ) FINAL_TIME="1296000\.0";;
		* )
		echo "WARNING: We should not have gotten to this point."
		echo "Potential error in script at end-time parameters."
		echo "The current case will run with a default time_final."
		;;
	esac

	# Use sed to set l_stats to .false., debug_level to 0,
	# dt = 60.0 sec, and ending time of case to FINAL_TIME
	# which is set in the case statement above.
	cat $MODEL_FILE | sed 's/l_stats\s*=\s*.*/l_stats = \.false\./g' | sed 's/debug_level\s*=\s*.*/debug_level = 0/g' | sed 's/dtmain\s*=\s*.*/dtmain = 60\.0/g' | sed 's/dtclosure\s*=\s*.*/dtclosure = 60\.0/g' | sed 's/time_final\s*=\s*.*/time_final = '$FINAL_TIME'/g' >> 'clubb.in'

	echo "Running "$NUM_TESTS" runs of case: "${RUN_CASE[$x]}"."

	# Run collect on the CLUBB model repeatedly
	# Produce the output within the desired output directory.

	for run in `seq $NUM_TESTS`; do
		COLLECT_OUT=$OUTPUT_DIR"/"${RUN_CASE[$x]}"-"$run".er"
		collect -p hi -A copy -o $COLLECT_OUT ../bin/clubb_standalone
	done

	# Remove the namelists
	rm -f 'clubb.in'
done

echo "Testing is complete. The output from the tests is in:"
echo $OUTPUT_DIR"/<case>-<run number>.er/ ."
echo "Run analyzer "$OUTPUT_DIR"/<case>-<run number>.er to view the results."

# Return to previous working directory.
cd $restoreDir
