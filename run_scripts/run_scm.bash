#!/bin/bash
#######################################################################
# $Id: run_standalone.bash,v 1.29 2008-04-17 00:02:03 dschanen Exp $
#
# Description:
# Script to run the standalone CLUBB program.  
# Tested with GNU Bash v2 & 3.  Might work with Ksh.
#
#######################################################################
# Useful variable on multiprocessor machines with OpenMP capable 
# Fortran, uncomment and set to use.  
# The hoc_standalone program is not parallel, but many LAPACK/BLAS 
# libraries are.
#
# export OMP_NUM_THREADS=2
#######################################################################

ZT_GRID_TEST=false
ZM_GRID_TEST=false

# Figure out the directory where the script is located
scriptPath=`dirname $0`

# Store the current directory location so it can be restored
restoreDir=`pwd`

# Change directories to the one the script is located in
cd $scriptPath

if [ -z $1 ]; then
	echo "Usage: "$0" [--zm_grid_test OR --zt_grid_test grid_levels grid_file]" \
			" <MODEL CASE> [PARAMETER FILE] [STATS FILE]"
	exit
else
	# Loop through to pull out options passed into the run.
	# This is done as a loop to facilitate additional options later on,
	# if necessary.
	while [ -n "$(echo $1 | grep "-")" ]; do
		case $1 in
			--zt_grid_test ) ZT_GRID_TEST=true
				if [ "$2" == "" ]; then
					echo "Option '--zt_grid_test': The number of levels in the grid" \
					"needs to be entered following '--zt_grid_test'."
					exit 1
				elif [ -n "$(echo $2 | grep "-")" ]; then
					echo "Option '--zt_grid_test': The number of grid levels" \
					"needs to follow '--zt_grid_test', not another option."
					exit 1
				elif [ "$3" == "" ]; then
					echo "Option '--zt_grid_test': The path to the grid file" \
					"needs to be entered following '--zt_grid_test'."
					exit 1
				elif [ -n "$(echo $3 | grep "-")" ]; then
					echo "Option '--zt_grid_test': The path to the grid file" \
					"needs to be entered following '--zt_grid_test'," \
					"not another option."
					exit 1
				elif [ $ZM_GRID_TEST == true ]; then
					echo "Only --zt_grid_test or --zm_grid_test may be used, not both."
					exit 1
				else
					test_grid_nz=$2
					test_grid_name=$3
					test_grid_name=`echo $test_grid_name | sed 's/\//\\\\\//g'`
					echo "Running case using specified zt grid."
					shift
					shift
				fi;;

			--zm_grid_test ) ZM_GRID_TEST=true
				if [ "$2" == "" ]; then
					echo "Option '--zm_grid_test': The number of levels in the grid" \
					"needs to be entered following '--zm_grid_test'."
					exit 1
				elif [ -n "$(echo $2 | grep "-")" ]; then
					echo "Option '--zm_grid_test': The number of grid levels" \
					"needs to follow '--zm_grid_test', not another option."
					exit 1
				elif [ "$3" == "" ]; then
					echo "Option '--zm_grid_test': The path to the grid file" \
					"needs to be entered following '--zm_grid_test'."
					exit 1
				elif [ -n "$(echo $3 | grep "-")" ]; then
					echo "Option '--zm_grid_test': The path to the grid file" \
					"needs to be entered following '--zm_grid_test'," \
					"not another option."
					exit 1
				elif [ $ZT_GRID_TEST == true ]; then
					echo "Only --zt_grid_test or --zm_grid_test may be used, not both."
					exit 1
				else
					test_grid_nz=$2
					test_grid_name=$3
					test_grid_name=`echo $test_grid_name | sed 's/\//\\\\\//g'`
					echo "Running case using specified zm grid."
					shift
					shift
				fi;;
			esac
		# Shift moves the parameters up one. Example: $2 -> $1, etc.
		# This is so we only have to check $1 for every iteration.
		shift
	done

	MODEL_FILE='../input/case_setups/'$1'_model.in'
	RUN_CASE=$1

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
# Enable G95 Runtime option that sets uninitialized 
# memory to a NaN value
#######################################################################
G95_MEM_INIT="NAN"
export G95_MEM_INIT

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

##########################################################################
# Set up the namelists for the CLUBB model.
# If we are using the --zt_grid_test or --zm_grid_test options,
# we need to modify the input files according to the options passed in.
##########################################################################

if [ $ZT_GRID_TEST == true ]; then
	cat $PARAMS_FILE > 'clubb.in'
	cat $MODEL_FILE | sed -e 's/^nzmax\s*=\s*.*//g' \
		-e 's/^grid_type\s*=\s*.*//g' \
		-e 's/^zm_grid_fname\s*=\s*.*//g' \
		-e "s/^zt_grid_fname\s*=\s*.*//g" \
		-e 's/^\&model_setting/\&model_setting\n \
		nzmax = '$test_grid_nz'\n \
		zt_grid_fname ='\'$test_grid_name\''\n \
		grid_type = 2\n/g' >> 'clubb.in'
	cat $STATS_FILE >> 'clubb.in'
elif [ $ZM_GRID_TEST == true ]; then
	cat $PARAMS_FILE > 'clubb.in'
	cat $MODEL_FILE | sed -e 's/nzmax\s*=\s*.*//g' \
		-e 's/grid_type\s*=\s*.*//g' \
		-e 's/zt_grid_fname\s*=\s*.*//g' \
		-e 's/zm_grid_fname\s*=\s*.*//g'
		-e 's/^\&model_setting/\&model_setting\n \
		nzmax = '$test_grid_nz'\n \
		zm_grid_fname ='\'$test_grid_name\''\n \
		grid_type = 3\n/g' >> 'clubb.in'
	cat $STATS_FILE >> 'clubb.in'
else
	cat $PARAMS_FILE $MODEL_FILE $STATS_FILE > 'clubb.in'
fi

#######################################################################
# State which case is being run
#######################################################################
echo "Running" $RUN_CASE

# Run the CLUBB model
../bin/clubb_standalone 

# Remove the namelists
rm -f 'clubb.in'

cd $restoreDir
