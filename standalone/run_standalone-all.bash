#!/bin/bash
#######################################################################
# $Id: run_standalone-all.bash,v 1.35 2008-07-25 18:03:45 senkbeir Exp $
#
# Script to run the standalone hoc program for all models.
# Tested with bash v2.  Might work with Ksh.
#
#######################################################################
# Useful on multiprocessor machines with OpenMP capable Fortran
#export OMP_NUM_THREADS=2
#######################################################################

NIGHTLY=false

# This function reads all the arguments and sets variables that will be used
# later in the script.
set_args()
{
	# Loop through the list of arguments ($1, $2...). This loop ignores
	# anything not starting with '-'.
	while [ -n "$(echo $1 | grep "-")" ]; do
		case $1 in
			# '--nightly' sets the script to run the nightly version
			--nightly ) NIGHTLY=true;;
			--help | -h | -? | * ) echo -e "Usage:\n  run_standalone-all.bash [OPTION]..."
					       echo "Options:"
					       echo -e "  --nightly\tPerforms the nightly run."
					       echo -e "  -h, --help\tShows help (this)."
			                       exit 1;;
		esac
		# Shift moves the parameters up one. Ex: $2 -> $1 and so on.
		# This is so we only have to check $1 on each iteration.
		shift
	done
}

set_args $*

if [ $NIGHTLY == true ] ; then
	echo -e "\nPerforming nightly run...\n"
else
	echo -e "\nPerforming standard run (all models)\n"
fi

EXIT_CODE=( [0]=0 [1]=0 [2]=0 [3]=0 [4]=0 [5]=0 [6]=0 [7]=0 [8]=0 [9]=0 \
	    [10]=0 [11]=0 [12]=0 [13]=0 [14]=0 [15]=0 [16]=0 [17]=0 [18]=0 )

RUN_CASE=( \
	arm atex bomex clex9_nov02 clex9_oct14 cobra dycoms2_rf01
        dycoms2_rf02_do dycoms2_rf02_ds	dycoms2_rf02_nd dycoms2_rf02_so \
        fire gabls2 jun25_altocu mpace_a mpace_b nov11_altocu rico wangara )

# Since everyone seems to like to add new cases without adding exit codes,
# we try and catch that error here...
if [ "${#RUN_CASE[@]}" -ne "${#EXIT_CODE[@]}" ] ; then
	echo "RUN_CASE: ${#RUN_CASE[@]}" "EXIT_CODE: ${#EXIT_CODE[@]}" 
	echo "EXIT_CODE is not equal in size to RUN_CASE"
	exit 1
fi

if [ $NIGHTLY == true ] ; then
	# Eliminate the previous CLUBB results.
	# This prevents spurious profile generation resulting from
	# previous profiles not getting overwritten
	rm -f CLUBB_previous/*

	mv CLUBB_current/*.ctl CLUBB_previous/
	mv CLUBB_current/*.dat CLUBB_previous/
fi
# This will loop over all runs in sequence 
for (( x=0; x < "${#RUN_CASE[@]}"; x++ )); do
#######################################################################
# Check for necessary namelists.  If files exist, then
# copy them over to the general input files.

	STANDALONE_IN='standalone_'"${RUN_CASE[$x]}"'.in'
	MODEL_IN='../model/'"${RUN_CASE[$x]}"'_model.in'
	if [ $NIGHTLY == true ] ; then
		STATS_IN='../stats/nightly_stats.in'
#		STATS_IN='../stats/nobudgets_stats.in'
	else
		STATS_IN='../stats/nobudgets_stats.in'
	fi

	if [ ! -e "$STANDALONE_IN" ] ; then
		echo $STANDALONE_IN " does not exist"
		exit 1
	fi

	if [ -e 'standalone.in' ] ; then
		rm -f 'standalone.in'
	fi

	ln -s $STANDALONE_IN 'standalone.in'

	if [ $NIGHTLY == true ] ; then
		# This is needed because the model file now contains stats_tout
		# Here we replace the repository version of stats_tout with an hour output
		# The regular expression use here matches:
		# 'stats_tout' (0 or > whitespaces) '=' (0 or > whitespaces) (0 or > characters)
		# and replaces it with 'stats_tout = 3600.'
		cat $MODEL_IN | sed 's/stats_tout\s*=\s*.*/stats_tout = 3600\./g' > "${RUN_CASE[$x]}"'_hoc.in'
		cat $STATS_IN >> "${RUN_CASE[$x]}"'_hoc.in'
	else
		cat $MODEL_IN $STATS_IN > "${RUN_CASE[$x]}"'_hoc.in'
	fi

	#######################################################################
	#
	# State which case is being run
	echo "Running ""${RUN_CASE[$x]}"
	# Run HOC 
	RESULT=`../bin/hoc_standalone 2>&1 |grep 'normal'`
	if [ -z "$RESULT" ]; then
		EXIT_CODE[$x]=-1
	fi

	# remove the namelists
	rm -f 'standalone.in'
	rm -f "${RUN_CASE[$x]}"'_hoc.in'
done

# Print the results
for (( x=0; x < "${#RUN_CASE[@]}"; x++ )); do
	if [ "${EXIT_CODE[$x]}" != 0 ]; then
		echo "${RUN_CASE[$x]}"' failure'

		if [ $NIGHTLY == true ] ; then
			rm "${RUN_CASE[$x]}"_zt.ctl
			rm "${RUN_CASE[$x]}"_zt.dat
			rm "${RUN_CASE[$x]}"_zm.ctl
			rm "${RUN_CASE[$x]}"_zm.dat
		fi
	else
		if [ $NIGHTLY == true ] ; then
			mv "${RUN_CASE[$x]}"_zt.ctl CLUBB_current/
			mv "${RUN_CASE[$x]}"_zt.dat CLUBB_current/
			mv "${RUN_CASE[$x]}"_zm.ctl CLUBB_current/
			mv "${RUN_CASE[$x]}"_zm.dat CLUBB_current/
		fi
 	fi
done
