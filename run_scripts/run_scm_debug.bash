#!/bin/bash
#######################################################################
# $Id$
#
# Description:
# Script to run the standalone CLUBB program with a debugger.  
# Tested with GNU Bash v2 & 3.  Might work with Ksh.
#
#######################################################################

# Figure out the directory where the script is located
scriptPath=`dirname $0`

# Store the current directory location so it can be restored
restoreDir=`pwd`

# Executable for the debugging program
#DEBUG=dbx
DEBUG=pgdbg

# Change directories to the one the script is located in
cd $scriptPath

if [ -z $1 ]; then
	echo "Usage: "$0" <MODEL CASE> [PARAMETER FILE] [STATS FILE]"
	exit
else
	MODEL_FILE='../input/case_setups/'$1'_model.in'

	if [ -z $2 ]; then
		PARAMS_FILE="../input/tunable_parameters.in"
	else
		PARAMS_FILE=$2
	fi

	if [ -z $3 ]; then
		STATS_FILE="../input/stats/all_stats.in"
	else
		STATS_FILE=$3
	fi
fi

#######################################################################
# Check for necessary namelists.  If files exist, then
# copy them over to the general input files.
#######################################################################

if [ ! -e "$MODEL_FILE" ]; then
	echo "$MODEL_FILE does not exist"
	exit 1
fi

if [ ! -e "$PARAMS_FILE" ]; then
	echo "$PARAMS_FILE does not exist"
	exit 1
fi

if [ ! -e "$STATS_FILE" ]; then
	echo "$STATS_FILE does not exist"
	exit 1
fi

cat $PARAMS_FILE $MODEL_FILE $STATS_FILE > 'clubb.in'

#######################################################################
# State which case is being run
#######################################################################
echo "Running" $RUN_CASE

# Run the CLUBB model
$DEBUG ../bin/clubb_standalone 

# Remove the namelists
rm -f 'clubb.in'

cd $restoreDir
