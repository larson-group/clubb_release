#!/bin/bash
##########################################################################
# $Id$
#
# Description:
# Script that runs a predefined set of test cases for creating an overall
#   timing profile.
#
# This script will run six different cases multiple times to produce a timing
#   profile for CLUBB. Each of the cases is run with the statistics output
#   disabled, debug level set to 0, and otherwise default settings. The
#   multiple runs can then be aggregated together to provide a more robust
#   statistical sample of the cases' timings.
#
# Additional cases can be run by adding these cases to the RUN_CASE array
#   declared at the top of the script. Upon adding a case to the RUN_CASE
#   array, a value for the number of times to run that case must be added to
#   the NUM_TESTS array. Ideally, the number of times a case is run should
#   be a large enough number so that the total runtime for that case is
#   approximately equal to 7-8 minutes. This is a large enough sample time
#   to ensure reasonably accurate data, but won't take a large amount of time
#   to run the tests.
#
# WARNING: A successful run of this script will likely consume approximately
#   25-30 GB of hard drive space! Please plan accordingly to ensure that the
#   output does not fill up the hard drive on the host computer!
#
# NOTE: When running this script, it is highly recommended that the
#   computer is doing nothing else--having any other programs running,
#   even an instance of GVim, will likely increase the inconsistency
#   of the results. Ideally, try to ensure that the only program running
#   on that computer besides the script is your SSH session/terminal window
#   that started the script, at the most. This includes programs that are run
#   by other users.
#
##########################################################################

RUN_CASE=( dycoms2_rf02_do mpace_a cobra gabls3_night jun25_altocu rico )

NUM_TESTS=( "105" "16" "46" "48" "155" "215" )

# Figure out the directory where the script is located
scriptPath=`dirname $0`

# Store the current working directory so it can be restored
restoreDir=`pwd`

# Change directories to the script's directory
cd $scriptPath

if [ -z $1 ]; then
	echo "Usage: "$0" <OUTPUT DIRECTORY>"
	echo "The output directory will contain the output from"
	echo "all of the profile tests."
	exit
else
	OUTPUT_DIR=$1
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

echo "Running test suite..."

##########################################################################
# Loop through the cases in the suite, setting up the parameters so testing
# runs can be completed
##########################################################################
for (( x=0; x < "${#RUN_CASE[@]}"; x++ )); do
	echo "---------------------------------------------"
	echo "Performing tests for case: "${RUN_CASE[$x]}"."
	echo "Tests will run "${NUM_TESTS[$x]}" times."
	echo "---------------------------------------------"

	# Create testing output subdirectory
	mkdir $OUTPUT_DIR"/${RUN_CASE[$x]}"

	# Ensure necessary namelists exist and copy them to the general
	# input files.

	MODEL_FILE="../input/case_setups/"${RUN_CASE[$x]}"_model.in"
	PARAMS_FILE="../input/tunable_parameters/tunable_parameters.in"
    SILHS_PARAMS_FILE="../input/tunable_parameters/silhs_parameters.in"
	FLAGS_FILE="../input/tunable_parameters/configurable_flags.in"

	if [ ! -e "$MODEL_FILE" ]; then
		echo "ERROR: $MODEL_FILE does not exist."
		exit 1
	fi

	if [ ! -e "$PARAMS_FILE" ]; then
		echo "ERROR: $PARAMS_FILE does not exist."
		exit 1
	fi

	if [ ! -e "$FLAGS_FILE" ]; then
		echo "ERROR: $FLAGS_FILE does not exist."
		exit 1
	fi

	cat $PARAMS_FILE $SILHS_PARAMS_FILE > "clubb.in"
	cat $FLAGS_FILE >> "clubb.in"

	# Use sed to set l_stats to .false. and debug_level to 0.
	cat $MODEL_FILE | sed 's/l_stats\s*=\s*.*/l_stats = \.false\./g' | sed 's/debug_level\s*=\s*.*/debug_level = 0/g' >> 'clubb.in'

	# Run collect on the CLUBB model repeatedly
	# Produce the output within the desired output directory.

	for run in `seq ${NUM_TESTS[$x]}`; do
		COLLECT_OUT=$OUTPUT_DIR"/"${RUN_CASE[$x]}"/"$run".er"
		collect -p hi -A copy -o $COLLECT_OUT ../bin/clubb_standalone
	done

	# Remove the namelists
	rm -f 'clubb.in'
done

echo "Testing is complete. The output from the tests is in:"
echo $OUTPUT_DIR"/<case>/<run number>.er/ ."
echo "Run analyzer "$OUTPUT_DIR"/<case>-<run number>.er to view the results."
echo "For best results, it is recommended that all runs for a particular"
echo "case are opened simulaneously in analyzer. Use the Open Case command"
echo "in the analyzer window, and highlight all runs for a particular case,"
echo "and the analyzer tool will aggregate all of the results together."

# Return to previous working directory.
cd $restoreDir
