#!/bin/bash
#######################################################################
# $Id$
#
# Script to run the standalone hoc program for all models.
# Tested with bash v2.  Might work with Ksh.
#
#######################################################################
# Useful on multiprocessor machines with OpenMP capable Fortran
#export OMP_NUM_THREADS=2
#######################################################################

NIGHTLY=false
TIMESTEP_TEST=false

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
                        # '--timestep_test' runs all cases at specified timestep.
                        # The specified timestep is listed right after '--timestep_test'
                        # and is entered into the code as "test_ts".
                        --timestep_test ) TIMESTEP_TEST=true
                                          if [ "$2" == "" ]; then
                                             # If '--timestep_test' is used, then a time step
                                             # length needs to be declared immediately.
                                             echo "Option '--timestep_test':  the time step length"\
                                                  "needs to be entered following '--timestep_test'."
                                             exit 1
                                          elif [ -n "$(echo $2 | grep "-")" ]; then
                                             # The time step length, and not another option, needs
                                             # to follow the '--timestep_test' option.
                                             echo "Option '--timestep_test':  the time step length"\
                                                  "needs to follow '--timestep_test', not another"\
                                                  "option."
                                             exit 1
                                          else
                                             # The time step length is entered as variable 'test_ts'.
                                             test_ts=$2
                                             echo ""
                                             echo "Running all cases at a(n) "$test_ts" second"\
                                                  "time step."
                                             shift
                                          fi;;
			--help | -h | -? | * ) echo -e "Usage:\n  run_standalone-all.bash [OPTION]..."
					       echo "Options:"
					       echo -e "  --nightly\t\t\t\tPerforms the nightly run."
					       echo -e "  --timestep_test  time_step_length\tRuns all"\
                                                          "cases at the specified time step length (sec)."
					       echo -e "  -h, --help\t\t\t\tShows help (this)."
			                       exit 1;;
		esac
		# Shift moves the parameters up one. Ex: $2 -> $1 and so on.
		# This is so we only have to check $1 on each iteration.
		shift
	done
}

run_case()
{
        #######################################################################
        # Enable G95 Runtime option that sets uninitialized 
        # memory to a NaN value
        #######################################################################
        G95_MEM_INIT="NAN"
        export G95_MEM_INIT

	#######################################################################
	#
	# State which case is being run
	echo "Running ""${RUN_CASE[$x]}"
	# Run HOC 
	#RESULT=`../bin/clubb_standalone 2>&1 |grep 'normal'`
	RESULT=`../bin/clubb_standalone 2>&1`

	if [ $NIGHTLY == true ]; then  
		echo -e "$RESULT";
	fi

	RESULT=`echo "$RESULT" | grep 'normal'`
	if [ -z "$RESULT" ]; then
		EXIT_CODE[$x]=-1
	fi

	# remove the namelists
	rm -f 'clubb.in'
}

set_args $*

if [ $NIGHTLY == true ] ; then
	echo -e "\nPerforming nightly run...\n"
else
	echo -e "\nPerforming standard run (all cases)\n"
fi

EXIT_CODE=( [0]=0 [1]=0 [2]=0 [3]=0 [4]=0 [5]=0 [6]=0 [7]=0 [8]=0 [9]=0 \
	    [10]=0 [11]=0 [12]=0 [13]=0 [14]=0 [15]=0 [16]=0 [17]=0 [18]=0 [19]=0 [20]=0 )

RUN_CASE=( \
	arm arm_97 atex bomex clex9_nov02 clex9_oct14 cobra dycoms2_rf01
        dycoms2_rf02_do dycoms2_rf02_ds	dycoms2_rf02_nd dycoms2_rf02_so \
        fire gabls2 jun25_altocu lba mpace_a mpace_b nov11_altocu rico wangara )

# Since everyone seems to like to add new cases without adding exit codes,
# we try and catch that error here...
if [ "${#RUN_CASE[@]}" -ne "${#EXIT_CODE[@]}" ] ; then
	echo "RUN_CASE: ${#RUN_CASE[@]}" "EXIT_CODE: ${#EXIT_CODE[@]}" 
	echo "EXIT_CODE is not equal in size to RUN_CASE"
	exit 1
fi

if [ $NIGHTLY == true ] ; then
	# Make the CLUBB_previous and CLUBB_current directories if they don't exist
	mkdir CLUBB_current
	mkdir CLUBB_previous
	
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

	#STANDALONE_IN='standalone_'"${RUN_CASE[$x]}"'.in'
	PARAMS_IN='default_parameters.in'
	MODEL_IN='../model/'"${RUN_CASE[$x]}"'_model.in'
	if [ $NIGHTLY == true ] ; then
		STATS_IN='../stats/nightly_stats.in'
#		STATS_IN='../stats/nobudgets_stats.in'
	else
		STATS_IN='../stats/nobudgets_stats.in'
	fi

	if [ ! -e $PARAMS_IN ] ; then
		echo $PARAMS_IN " does not exist"
		exit 1
	fi


	if [ $NIGHTLY == true ] ; then
		cat $PARAMS_IN > 'clubb.in'
		# This is needed because the model file now contains stats_tout
		# Here we replace the repository version of stats_tout with an hour output
		# The regular expression use here matches:
		# 'stats_tout' (0 or > whitespaces) '=' (0 or > whitespaces) (0 or > characters)
		# and replaces it with 'stats_tout = 3600.'
		cat $MODEL_IN | sed 's/stats_tout\s*=\s*.*/stats_tout = 3600\./g' >> 'clubb.in'
		cat $STATS_IN >> 'clubb.in'
		run_case

		#Move the ZT and SM files out of the way
		if [ "${EXIT_CODE[$x]}" != 0 ]; then
			rm "${RUN_CASE[$x]}"_zt.ctl
			rm "${RUN_CASE[$x]}"_zt.dat
			rm "${RUN_CASE[$x]}"_zm.ctl
			rm "${RUN_CASE[$x]}"_zm.dat
			rm "${RUN_CASE[$x]}"_sfc.ctl
			rm "${RUN_CASE[$x]}"_sfc.dat
		else
			mv "${RUN_CASE[$x]}"_zt.ctl CLUBB_current/
			mv "${RUN_CASE[$x]}"_zt.dat CLUBB_current/
			mv "${RUN_CASE[$x]}"_zm.ctl CLUBB_current/
			mv "${RUN_CASE[$x]}"_zm.dat CLUBB_current/
			rm "${RUN_CASE[$x]}"_sfc.ctl
			rm "${RUN_CASE[$x]}"_sfc.dat
		fi
		
		#Run again with a finer time step
		cat $PARAMS_IN > 'clubb.in'
		cat $MODEL_IN | sed 's/stats_tout\s*=\s*.*/stats_tout = 60\./g' >> 'clubb.in'
		cat $STATS_IN >> 'clubb.in'
		run_case

		#Now move the SFC file
		if [ "${EXIT_CODE[$x]}" != 0 ]; then
			rm "${RUN_CASE[$x]}"_zt.ctl
			rm "${RUN_CASE[$x]}"_zt.dat
			rm "${RUN_CASE[$x]}"_zm.ctl
			rm "${RUN_CASE[$x]}"_zm.dat
			rm "${RUN_CASE[$x]}"_sfc.ctl
			rm "${RUN_CASE[$x]}"_sfc.dat
		else
			rm "${RUN_CASE[$x]}"_zt.ctl
			rm "${RUN_CASE[$x]}"_zt.dat
			rm "${RUN_CASE[$x]}"_zm.ctl
			rm "${RUN_CASE[$x]}"_zm.dat
			mv "${RUN_CASE[$x]}"_sfc.ctl CLUBB_current/
			mv "${RUN_CASE[$x]}"_sfc.dat CLUBB_current/
		fi

	elif [ $TIMESTEP_TEST == true ]; then

                # Set the model timestep for all cases (and the stats output timestep
                # unless l_stats is overwritten to .false.) to timestep test_ts.
                cat $PARAMS_IN > 'clubb.in'
                cat $MODEL_IN | sed -e 's/dtmain\s*=\s*.*/dtmain = '$test_ts'/g' \
                                    -e 's/dtclosure\s*=\s*.*/dtclosure = '$test_ts'/g' \
                                    -e 's/stats_tsamp\s*=\s*.*/stats_tsamp = '$test_ts'/g' \
                                    -e 's/stats_tout\s*=\s*.*/stats_tout = '$test_ts'/g' >> 'clubb.in'
                cat $STATS_IN >> 'clubb.in'
                run_case

        else

		cat $PARAMS_IN $MODEL_IN $STATS_IN > 'clubb.in'
		run_case

	fi

done

# Print the results and copy files for a nightly run
for (( x=0; x < "${#RUN_CASE[@]}"; x++ )); do
	if [ "${EXIT_CODE[$x]}" != 0 ]; then
		echo "${RUN_CASE[$x]}"' failure'
 	fi
done
