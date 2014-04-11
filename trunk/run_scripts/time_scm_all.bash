#!/bin/bash
#####################################################################
# $Id$
#
# Script that times a run of the standalone CLUBB program for all models.
# Utilizes the -e (--performance_test) option of run_scm(_all).bash, and
#   the time system command.
#
# Accepts one optional argument, which pipes the output of the time
#   command to a file rather than output to standard out. This is because
#   the time command's output behaves oddly, and needs to be piped in a
#   special manner.
#####################################################################

RUN_CASE=( \
	arm atex bomex dycoms2_rf01 dycoms2_rf02_nd dycoms2_rf02_so \
	fire gabls2 mpace_b wangara cobra gabls3_night jun25_altocu )

# Back up the old RUN_CASES file...
mv RUN_CASES RUN_CASES.bak

# Populate the "new" RUN_CASES file...
for (( x=0; x < "${#RUN_CASE[@]}"; x++ )); do
	echo "${RUN_CASE[$x]}" >> "RUN_CASES"
done

# As the time for a run is now much shorter, this script aggregates multiple
# runs together. As such, it is required that the number of runs is specified.
if [ -z $1 ]; then
	echo "Usage: time_scm_all.bash <NUMBER OF RUNS> [OUTPUT FILE]"
	echo "Note: Number of runs is a required argument."
	echo "The recommended number of runs depends on the compiler:"
	echo -e "\tg95: approximately 8 runs"
	echo -e "\tIntel Fortran: approximately 15 runs"
	echo -e "\tSun Studio Fortran: approximately 60 runs"
	exit 1
else
	NUM_RUNS=$1
fi

echo "Running timing test on a modified run_scm_all.bash..."
echo
echo "The only cases that will be run are those that not use"
echo "microphysics or non-simplified radiation schemes."

if [ -z $2 ]; then
	echo
	echo "NOTE: If you wish to have the output piped to a file, perform:"
	echo -e "\t"$0" <NUMBER OF RUNS> <OUTPUT FILE>"
	echo
else
	OUTPUT_FILE=$2

	echo
	echo "Output will been piped to "$OUTPUT_FILE"."
	echo "If "$OUTPUT_FILE" exists, it will be overwritten."
	echo
	echo "Timing test, with "$NUM_RUNS" runs:" > $OUTPUT_FILE
fi

totaltime=0

for (( x=0; x < $NUM_RUNS; x++ )) do
	exec 3>&1 4>&2
	cputime=$( TIMEFORMAT="%U"; { time ./run_scm_all.bash -e 1>&3 2>&4; } 2>&1)
	exec 3>&- 4>&-
	if [ -z $2 ]; then
		echo "Run "$x" time: "$cputime"s"
	else
		echo "Run "$x" time: "$cputime"s" >> $OUTPUT_FILE
	fi
	cputime=${cputime/./}
	let totaltime=$totaltime+$cputime
done

#Perform calculations to get the relevant pieces of information
#from the total time in milliseconds, for display to the user.
let totalseconds=$(($totaltime/1000))
let totalmillis=$(($totaltime%1000))
let totaltenths=$(($totalmillis/100))
let totalhundths=$(( ( $totalmillis/10 ) % 10 ))
let totalmills=$(($totalmillis%10))

if [ -z $2 ]; then
	echo "Timing test, with "$NUM_RUNS" runs:"
	echo "Total User CPU Time: "$totalseconds"."$totaltenths""$totalhundths""$totalmills"s"
else
	echo "Total User CPU Time: "$totalseconds"."$totaltenths""$totalhundths""$totalmills"s" >> $OUTPUT_FILE
fi

# Restore the backup of the RUN_CASES file
mv RUN_CASES.bak RUN_CASES
