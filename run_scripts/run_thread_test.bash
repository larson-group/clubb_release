#!/bin/bash
###############################################################################
# $Id$
#
# Desciption:
#   Script to run the test whether CLUBB is still threadsafe.
#   We need this because not all versions of Fortran accept command line
#   arguments (a Fortran 2003 feature), and GNU Fortran has trouble with
#   comments in namelists
# Notes:
#   This has been tested with GNU Bash v4.1.2 on GNU/Linux.  The sed command
#   used only works with GNU sed.
###############################################################################

# Set the number of threads to
export OMP_NUM_THREADS=8

RUN_CASES=( fire bomex rico_lh gabls3 )
NAMELISTS=( 'clubb_1.in' 'clubb_2.in' 'clubb_3.in' 'clubb_4.in' )
FLAGS_FILE="../input/tunable_parameters/configurable_model_flags.in"
PARAMS_FILE="../input/tunable_parameters/tunable_parameters.in"
SILHS_PARAMS_FILE="../input/tunable_parameters/silhs_parameters.in"
STATS_FILE="../input/stats/standard_stats.in"
SERIAL="../output/serial"
PARALLEL="../output/parallel"

run_serial()
{
	echo -n "Running CLUBB in serial... "
	for (( x=0;  x < "${#RUN_CASES[@]}"; x++ )); do
		MODEL_FILE='../input/case_setups/'"${RUN_CASES[$x]}"'_model.in'
		cat $MODEL_FILE $FLAGS_FILE $PARAMS_FILE $SILHS_PARAMS_FILE $STATS_FILE > 'clubb.in'
 

		# Set debug level to 0
#		sed -i 's/debug_level\s*=\s*.*/debug_level = 0/g' 'clubb.in'

		# This is a kluge for Fortran compilers that the can't handle comments in 
		# a namelist by using the sed command to remove them.
		# Since -i is a GNU sed extension the command might be 'gsed' on some systems.
		sed -i -e 's/\!.*//' 'clubb.in'

		if [ -e ../bin/clubb_standalone ]; then
			../bin/clubb_standalone &> /dev/null
		else
			echo "clubb_standalone not found (did you re-compile?)"
			exit
		fi

		# Remove the namelist
		rm -f 'clubb.in'
	done
	echo "Done"'!'
}

run_parallel()
{
	echo -n "Running CLUBB in parallel... "
	for (( x=0;  x < "${#RUN_CASES[@]}"; x++ )); do
		MODEL_FILE='../input/case_setups/'"${RUN_CASES[$x]}"'_model.in'
		cat $MODEL_FILE $FLAGS_FILE $PARAMS_FILE $SILHS_PARAMS_FILE $STATS_FILE > "${NAMELISTS[$x]}"


		# Set debug level to 0
#		sed -i 's/debug_level\s*=\s*.*/debug_level = 0/g' "${NAMELISTS[$x]}"

		# This is a kluge for Fortran compilers that the can't handle comments in 
		# a namelist by using the sed command to remove them.
		# Since -i is a GNU sed extension the command might be 'gsed' on some systems.
		sed -i -e 's/\!.*//' "${NAMELISTS[$x]}"
	done

	# Run the CLUBB thread test
	if [ -e ../bin/clubb_thread_test ]; then
		../bin/clubb_thread_test &> /dev/null
	else
		echo "clubb_thread_test not found (did you re-compile?)"
		exit
	fi

	# Remove the namelists
	for (( x=0;  x < "${#RUN_CASES[@]}"; x++ )); do
		rm -f "${NAMELISTS[$x]}"
	done
	echo 'Done!'
}

# Figure out the directory where the script is located
scriptPath=`dirname $0`

# Store the current directory location so it can be restored
restoreDir=`pwd`

# Change directories to the one the script is located in
cd $scriptPath

# Run in serial mode
run_serial

mkdir $SERIAL
mv ../output/*.??? $SERIAL

# Run in parallel mode
run_parallel

mkdir $PARALLEL
mv ../output/*.??? $PARALLEL

echo -n "Diffing the output... "
for (( x=0;  x < "${#RUN_CASES[@]}"; x++ )); do
	diff $SERIAL/"${RUN_CASES[$x]}"_zt.dat $PARALLEL/"${RUN_CASES[$x]}"_zt.dat &>> 'diff.txt'
	diff $SERIAL/"${RUN_CASES[$x]}"_zm.dat $PARALLEL/"${RUN_CASES[$x]}"_zm.dat &>> 'diff.txt'
	diff $SERIAL/"${RUN_CASES[$x]}"_sfc.dat $PARALLEL/"${RUN_CASES[$x]}"_sfc.dat &>> 'diff.txt'
done
echo -n "Done"
echo '!'

if [[ -s 'diff.txt' ]]; then
	cat 'diff.txt'
	RESULT=1
else
	echo "No differences found"
	RESULT=0
fi

setup_files=`ls $SERIAL/*_setup.txt |wc -l`
if [ "${#RUN_CASES[@]}" -ne "$setup_files" ]; then
	echo "One or more simulations failed to run in serial"
	RESULT=1
fi

setup_files=`ls $PARALLEL/*_setup.txt |wc -l`
if [ "${#RUN_CASES[@]}" -ne "$setup_files" ]; then
	echo "One or more simulations failed to run in parallel"
	RESULT=1
fi

rm -f 'diff.txt'
rm -rf $SERIAL
rm -rf $PARALLEL

cd $restoreDir

exit $RESULT
