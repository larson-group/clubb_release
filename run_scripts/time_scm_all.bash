#!/bin/bash
#####################################################################
# $Id:$
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

echo "Running timing test on run_scm_all.bash..."

if [ -z $1 ]; then
	echo "NOTE: If you wish to have the output piped to a file, perform:"
	echo -e "\t"$0" <OUTPUT FILE>"

	{ time ./run_scm_all.bash -e ; }
else
	OUTPUT_FILE=$1

	echo "Output will be piped to "$OUTPUT_FILE"."
	echo "If "$OUTPUT_FILE" exists, it will be overwritten."

	{ time ./run_scm_all.bash -e ; } 2> $OUTPUT_FILE
fi
