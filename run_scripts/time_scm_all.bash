#!/bin/bash
#####################################################################
# $Id$
#
# Script that times a run of the standalone CLUBB program.
# Utilizes the run_scm.bash command with a -e flag to turn off stats
#   and debugging output, and uses the time command to collect data.
#
#####################################################################

CASES=( arm atex bomex fire wangara )

CASE_RUNS=( 1 1 4 8 4 )

#Declaration of default values
RUNS=3
TIMESTEP="NOT SPECIFIED"
PARAMETER_OPTIONS=b #b=r8029 backwards compatible parameters, c=current parameters

echo "" # Just for newline after command, I prefer how it looks

options=$( getopt -o d:p:h --long no-run,cpu_time: -- "$@" )

if [ $? != 0 ] ; then echo "Run with -h for help." >&2 ; exit 1 ; fi

eval set -- "$options"

while true; do
	case "$1" in
		-d )
            shift;
            if [[ -z "$( echo $1 | sed 's:[0-9]\+\.[0-9]\+::g')" ]] ; then
                 TIMESTEP=$1
                 sed -i "s:dt_main\s*=\s*[0-9.]*:dt_main = $1:g" ../input/case_setups/*_model.in
                 sed -i "s:dt_rad\s*=\s*[0-9.]*:dt_rad = $1:g" ../input/case_setups/*_model.in
            else
                 echo " The timestep specified by -d [timesetp] must be in the form of a number with a decimal. For example -d 1.0 is valid, but -d 1 is invalid."
                 exit 1
            fi
			;;
            
        -p )
            shift;
            if [ "$1" = c ]; then
                PARAMETER_OPTIONS=c
            elif [ "$1" = b ]; then
                PARAMETER_OPTIONS=b
            else
                echo " -p accepts either c or b"
                echo "    OPTION=c   Use current settings"
                echo "    OPTION=b   Use r8029 backwards compatible settings"
                exit 1
            fi
            ;;
            
        --cpu_time )
            shift;
            if [ "$1" = y ] || [ "$1" = yes ]; then
                sed -i 's: !call cpu_time(: call cpu_time(:g' ../src/clubb_driver.F90
                echo -e " Uncommenting all calls to cpu_time within clubb_driver"
                echo -e "   ----- THIS REQUIRES A RECOMPILE TO TAKE EFFECT ----- \n"
            elif [ "$1" = n ] || [ "$1" = no ]; then
                sed -i 's: call cpu_time(: !call cpu_time(:g' ../src/clubb_driver.F90
                echo -e " Commenting out all calls to cpu_time within clubb_driver"
                echo -e "   ----- THIS REQUIRES A RECOMPILE TO TAKE EFFECT ----- \n"
            else
                echo " --cpu_time accepts either y/yes or n/no"
                echo "    OPTION=y/yes   Uncomment all calls to cpu_time"
                echo "    OPTION=n/no    Comment out all calls to cpu_time"
                exit 1
            fi
            ;;
            
        --no-run )
            RUNS=0
            ;;
	    
    	-h )
		    echo -e " usage: time_scm_all.bash [RUNS] [-c] [-d TIMESTEP] [-p OPTION] "
            echo -e "                          [--no-run] [--list-defaults] [-h] "
            echo
		    echo -e "  RUNS               Specifies number of timing runs. Default = 3"
            echo -e "  -d TIMESTEP        Changes timestep of model files to be run"
            echo -e "  -p OPTION          Changes the parameters and model flags files"
            echo -e "                        OPTION=c   Use current settings"
            echo -e "                        OPTION=b   Use r8029 backwards compatible settings"
            echo -e "  --cpu_time OPTION  Calls to cpu_time are expensive themselves, you may want"
            echo -e "                     comment them out to get more meaningful timing results."
            echo -e "                     This will require compilation to take effect."
            echo -e "                           OPTION=y/yes   Uncomment all calls to cpu_time"
            echo -e "                           OPTION=n/no    Comment out all calls to cpu_time"
            echo -e "  --no-run           Do not run a timing test. Essentially turning this script"
            echo -e "                        into a configuration tool"
		    echo -e "  -h                 Prints this help message"
            exit 0
		    ;;
            
         -- )
            shift
            break
            ;;
	esac
    shift
done


#Modify the model flags and parameters
CURRENT_PARAMS_DIR=../input/tunable_parameters
BACKWARDS_PARAMS_DIR=../input/tunable_parameters_compatible_r8029
SAVED_PARAMS_DIR=../input/tunable_parameters/save

if [ "$PARAMETER_OPTIONS" = c ]; then
    
    PARAMETER_DESCRIPTION="Using current parameters and model flags."
    echo -e " Setting parameters and model flags to use current settings\n"
    
    if [ -d $SAVED_PARAMS_DIR ]; then
        cp $SAVED_PARAMS_DIR/* $CURRENT_PARAMS_DIR
    fi
        
elif [ "$PARAMETER_OPTIONS" = b ]; then
    
    PARAMETER_DESCRIPTION="Using r8029 backwards compatible parameters and model flags."
    echo -e " Setting parameters and model flags to use backwards compatible settings\n"
    
    if [ ! -d $SAVED_PARAMS_DIR ]; then
        mkdir $SAVED_PARAMS_DIR
        cp $CURRENT_PARAMS_DIR/* $SAVED_PARAMS_DIR
    fi
    
    cp $BACKWARDS_PARAMS_DIR/* $CURRENT_PARAMS_DIR
    
fi


#Check if any runs will be completed
if [ $RUNS = 0 ]; then
    echo " Number of runs has been specified to be 0 with either --no-run or inputting 0."
    echo " There will be no timming runs."
    exit 0
fi


#Set number of runs
if [ -z $1 ]; then
    echo -e " No input for number of runs. Using the default: "$RUNS
	NUM_RUNS=$RUNS
else
    if [ -z "$( echo $1 | sed 's:[0-9]::g')" ];  then
	    NUM_RUNS=$1
    else 
        echo -e " Number of runs must be an integer. '$1' not accepted. Exiting...\n"
        exit 1
    fi
fi


echo -e "\n Running the case group with '"${CASES[*]}"' "$NUM_RUNS" times."
echo -e " Each case will be run '"${CASE_RUNS[*]}"' times respectively within each group."
echo -e " $PARAMETER_DESCRIPTION"
echo -e " Timestep = "$TIMESTEP"\n"

TIMES_GROUP=()
TIMES_CUMULATIVE=()
for (( i=0; i < $NUM_RUNS; i++ )) do
    
    TIMES_INDIVIDUAL=()

    for (( j=0; j < ${#CASES[@]}; j++ )) do
        
        if [ "$i" = 0 ]; then TIMES_CUMULATIVE+=( 0 ); fi
    
        echo -n " "${CASES[$j]}" x"${CASE_RUNS[$j]}" -- "
        
    	REAL_TIME=$( TIMEFORMAT="%U"; { time ( \
            for (( k=0; k < ${CASE_RUNS[$j]}; k++)) do
                ./run_scm.bash -e ${CASES[$j]} &>/dev/null 
            done ); } 2>&1 )
            
        echo $REAL_TIME"s"
        TIMES_INDIVIDUAL+=( $REAL_TIME )
        TIMES_CUMULATIVE[$j]=$( bc <<< "${TIMES_CUMULATIVE[$j]} + $REAL_TIME" )
        
    done
    
    TOTAL=$( IFS="+"; bc <<< "${TIMES_INDIVIDUAL[*]}" )
    TIMES_GROUP+=( $TOTAL )
    
    echo -e " Group "$((i+1))"/"$NUM_RUNS" Total -- "$TOTAL"s\n"
    
done

TOTAL=$( IFS="+"; bc <<< "${TIMES_GROUP[*]}" )

AVG=$( bc -l <<< "scale=3; $TOTAL/$NUM_RUNS")

ST_DEV=$(
            echo ${TIMES_GROUP[*]} | 
            awk "{  
                    for( i=1; i<=NF; i++ ) {
                        sum += ( \$i - $AVG )^2
                    }  
                    print sqrt( sum / NF )
                 }"
        )
echo    " ------------------------------- "
echo -e " Total Runtime -- "$TOTAL"s\n"

echo -e " Average times by case: "
for (( j=0; j < ${#CASES[@]}; j++ )) do

    echo "  "${CASES[$j]}" x"${CASE_RUNS[$j]}" -- "$( bc -l <<< "scale=3; ${TIMES_CUMULATIVE[$j]}/$NUM_RUNS" )"s"

done

echo -e "\n Average Group Runtime -- "$AVG"s"
echo -e " Group Standard Deviation -- "$ST_DEV"s\n"
