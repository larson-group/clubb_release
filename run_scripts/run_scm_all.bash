#!/bin/bash
#######################################################################
# $Id$
#
# Script to run the standalone CLUBB program for all models.
# Tested with bash v2.  Might work with Ksh.
#
#######################################################################
# Useful on multiprocessor machines with OpenMP capable Fortran
#export OMP_NUM_THREADS=2
#######################################################################

NIGHTLY=false
SHORT_CASES=false
PRIORITY_CASES=false
MIN_CASES=false

OUTPUT_DIR=$nightlyOut

# Figure out the directory where the script is located
scriptPath=`dirname $0`

# Store the current directory location so it can be restored
restoreDir=`pwd`

# Change directories to the one the script is located in
cd $scriptPath

OPTIONS=$*

# Note that we use `"$@"' to let each command-line parameter expand to a
# separate word. The quotes around `$@' are essential!
# We need TEMP as the `eval set --' would nuke the return value of getopt.
TEMP=`getopt -o :nhcij --long nightly,help,short-cases,priority-cases,min-cases -n 'run_scm_all.bash' -- "$@"`

# Note the quotes around `$TEMP': they are essential!
eval set -- "$TEMP"

while true ; do
  case "$1" in
  -n|--nightly) # Use nightly mode
            NIGHTLY=true
            shift ;;
  -h|--help) # Print the help message
            echo -e "Usage: run_scm_all.bash [OPTION]..."
            echo -e "\t-n, --nightly\t\t\tRun in nightly mode"
            echo -e "\t-c, --short-cases\t\tRun short cases. This will omit\n\t\tthe gabls2, cloud_feedback_s6, cloud_feedback_s11,\n\t\tcloud_feedback_s12, and twp_ice cases"
            echo -e "\t-i, --priority-cases\t\tRun priority cases. This will include\n\t\tonly the following cases if they are in RUN_CASES: arm, atex,\n\t\tbomex, dycoms2_rf01, dycoms2_rf01_fixed_sst, dycoms2_rf02_ds,\n\t\tdycoms2_rf02_nd, mpace_b, rico, wangara, arm_97,\n\t\tcloud_feedback_s6, cloud_feedback_s11, cloud_feedback_s12,\n\t\tgabls3_night, lba, and twp_ice."
            echo -e "\t-j, --min-cases\t\t\tRun a minimal set of cases (e.g. to\n\t\teconomize on output). This will include only the following\n\t\tcases if they are in RUN_CASES: arm, atex, bomex, dycoms2_rf01,\n\t\tdycoms2_rf02_ds, rico, wangara, arm_97, gabls3_night, lba,\n\t\tand twp_ice."
            echo -e "\t-h, --help\t\t\tPrints this help message"

            # Since the options for run_scm.bash are also valid, print those too:
            ./run_scm.bash --help | grep -v "Usage: run_scm.bash" | grep -v "\-\-help"

            exit 1 ;;
  -c|--short-cases) # Omit the longest cases to perform a shorter run
      SHORT_CASES=true
      echo "SHORT"
      shift ;;
  -i|--priority-cases) # Include only priority cases to filter unneeded data
      PRIORITY_CASES=true
      shift ;;
  -j|--min-cases) # A subset of priority cases
      MIN_CASES=true
      shift ;;
  --) shift ; break ;;
  *) echo "Something bad happened!" ; exit 1 ;;
  esac
done

# Declare arrays that will be used later on in this script.
declare -a RUN_CASE
declare -a EXIT_CODES
if [ $SHORT_CASES == true ] ; then # Run only short cases
    # Remove -c and --short-cases from the options so they aren't
    # passed to the run_scm.bash script
    OPTIONS=${OPTIONS#-c}
    OPTIONS=${OPTIONS#--short-cases}
    SHORT_CASES=true

    declare -a IGNORE_CASES

    # Populate the list of cases to ingore here. These ignored cases are the
    # cases that take the longest to run.
    IGNORE_CASES[0]=gabls2
    IGNORE_CASES[1]=cgils_s6
    IGNORE_CASES[2]=cgils_s11
    IGNORE_CASES[3]=cgils_s12
    IGNORE_CASES[4]=cloud_feedback_s6
    IGNORE_CASES[5]=cloud_feedback_s11
    IGNORE_CASES[6]=cloud_feedback_s12
    IGNORE_CASES[7]=twp_ice

    a=0
    while read line
    do  # If the line is not commented out (does not start with '!')
      if [[ $line != !* ]] && [[ ! -z $line ]] ; then
        # Check to see if this case should be ignored
        ignore=0
        for (( x=0; x<${#IGNORE_CASES[@]}; x++ )); do
            if [ $line == ${IGNORE_CASES[$x]} ] ; then
                ignore=1
            fi
        done

        # If the case was found in IGNORE_CASES, don't add it.
        if [ $ignore == 0 ]; then
          RUN_CASE[$a]=$line
          a=$(($a+1));
        fi
      fi
    done < "RUN_CASES"
elif [ $PRIORITY_CASES == true ] ; then
    # Remove -i and --priority-cases from the options so they aren't
    # passed to the run_scm.bash script
    OPTIONS=${OPTIONS#-i}
    OPTIONS=${OPTIONS#--priority-cases}
    PRIORITY_CASES=true

    declare -a PRIORITY_CASE_ARRAY

    PRIORITY_CASE_ARRAY[0]="arm"
    PRIORITY_CASE_ARRAY[1]="atex"
    PRIORITY_CASE_ARRAY[2]="bomex"
    PRIORITY_CASE_ARRAY[3]="dycoms2_rf01"
    PRIORITY_CASE_ARRAY[4]="dycoms2_rf01_fixed_sst"
    PRIORITY_CASE_ARRAY[5]="dycoms2_rf02_ds"
    PRIORITY_CASE_ARRAY[6]="dycoms2_rf02_nd"
    PRIORITY_CASE_ARRAY[7]="mpace_b"
    PRIORITY_CASE_ARRAY[8]="rico"
    PRIORITY_CASE_ARRAY[9]="wangara"
    PRIORITY_CASE_ARRAY[10]="arm_97"
    PRIORITY_CASE_ARRAY[11]="cgils_s6"
    PRIORITY_CASE_ARRAY[12]="cgils_s11"
    PRIORITY_CASE_ARRAY[13]="cgils_s12"
    PRIORITY_CASE_ARRAY[14]="gabls3_night"
    PRIORITY_CASE_ARRAY[15]="lba"
    PRIORITY_CASE_ARRAY[16]="twp_ice"

    a=0
    while read line
    do  # If the line is not commented out (does not start with '!')
      if [[ $line != !* ]] && [[ ! -z $line ]] ; then
        # Check to see if this case should be included
        include=0
        for (( x=0; x<${#PRIORITY_CASE_ARRAY[@]}; x++ )); do
            if [ $line == ${PRIORITY_CASE_ARRAY[$x]} ] ; then
                include=1
            fi
        done

        # If the case was found in PRIORITY_CASE_ARRAY, add it to the run cases.
        if [ $include == 1 ]; then
          RUN_CASE[$a]=$line
          a=$(($a+1));
        fi
      fi
    done < "RUN_CASES"
elif [ $MIN_CASES == true ] ; then
    OPTIONS=${OPTIONS#-j}
    OPTIONS=${OPTIONS#--min-cases}
    MIN_CASES=true
   
    declare -a MIN_CASE_ARRAY

    MIN_CASE_ARRAY[0]="arm"
    MIN_CASE_ARRAY[1]="atex"
    MIN_CASE_ARRAY[2]="bomex"
    MIN_CASE_ARRAY[3]="dycoms2_rf01"
    MIN_CASE_ARRAY[4]="dycoms2_rf02_ds"
    MIN_CASE_ARRAY[5]="rico"
    MIN_CASE_ARRAY[6]="wangara"
    MIN_CASE_ARRAY[7]="arm_97"
    MIN_CASE_ARRAY[8]="gabls3_night"
    MIN_CASE_ARRAY[9]="lba"
    MIN_CASE_ARRAY[10]="twp_ice"

    a=0
    while read line
    do  # If the line is not commented out (does not start with '!')
      if [[ $line != !* ]] && [[ ! -z $line ]] ; then
        # Check to see if this case should be included
        include=0
        for (( x=0; x<${#MIN_CASE_ARRAY[@]}; x++ )); do
            if [ $line == ${MIN_CASE_ARRAY[$x]} ] ; then
                include=1
            fi
        done

        # If the case was found in MIN_CASE_ARRAY, add it to the run cases.
        if [ $include == 1 ]; then
          RUN_CASE[$a]=$line
          a=$(($a+1));
        fi
      fi
    done < "RUN_CASES"
else # Populate the RUN_CASE array normally
    a=0
    while read line
    do
      # If the line is not commented out (does not start with '!')
      if [[ $line != !* ]] && [[ ! -z $line ]] ; then
        RUN_CASE[$a]=$line
        a=$(($a+1));
      fi
     done < "RUN_CASES"
fi

# Initialize all elements in EXIT_CODES to 0
# There will be one EXIT_CODE element for each RUN_CASE
for (( x=0; x < "${#RUN_CASE[@]}"; x++ )); do
    EXIT_CODES[$x]=0
done


if [ $NIGHTLY == true ] ; then
    echo -e "\nPerforming nightly run...\n"

    # Make the CLUBB_previous and CLUBB_current directories if they don't exist
    mkdir -p $OUTPUT_DIR"/CLUBB_current"
    mkdir -p $OUTPUT_DIR"/CLUBB_previous"

    # Eliminate the previous CLUBB results.
    # This prevents spurious profile generation resulting from
    # previous profiles not getting overwritten
    rm -f $OUTPUT_DIR"/CLUBB_previous/*"

    mv $OUTPUT_DIR/CLUBB_current/*.ctl $OUTPUT_DIR/CLUBB_previous/
    mv $OUTPUT_DIR/CLUBB_current/*.dat $OUTPUT_DIR/CLUBB_previous/
    mv $OUTPUT_DIR/CLUBB_current/*.nc $OUTPUT_DIR/CLUBB_previous/
    mv $OUTPUT_DIR/CLUBB_current/*.txt $OUTPUT_DIR/CLUBB_previous/
elif [ $SHORT_CASES == true ] ; then
    echo -e "\nPerforming short-cases run\n"
elif [ $PRIORITY_CASES == true ] ; then
    echo -e "\nPerforming priority-cases run\n"
elif [ $MIN_CASES == true ] ; then
    echo -e "\nPerforming min-cases run\n"
else
    echo -e "\nPerforming standard run\n"
fi

# This will loop over all runs in sequence
for (( x=0; x < "${#RUN_CASE[@]}"; x++ ));
do
    echo -e "Running ${RUN_CASE[$x]}"
    TMP_OUT=tmp_out.log

    if [ $NIGHTLY == true ] ; then
      ./run_scm.bash $OPTIONS ${RUN_CASE[$x]}
			EXIT_CODES[$x]=$?
    else
      # Send standard output to the bit bucket so only error output is recorded
      #RESULT=`./run_scm.bash $OPTIONS ${RUN_CASE[$x]} 2>&1 >/dev/null`
      ./run_scm.bash $OPTIONS ${RUN_CASE[$x]} 2>$TMP_OUT >/dev/null
			EXIT_CODES[$x]=$?
    fi

    #RESULT_STATUS=`echo "$RESULT" | grep 'normally'`
    RESULT_STATUS=`grep 'normally' $TMP_OUT`

    if [ ${EXIT_CODES[$x]} -ne 0 ]; then
        # If there was an error, and this is not running in nightly mode,
        # it will not be displayed. So, display the error here.
        if [ $NIGHTLY != true ] ; then
            #echo -e "$RESULT" | tail
            tail $TMP_OUT
        fi
    fi
    rm -f $TMP_OUT
done

EXIT_STATUS=0

# Print the results and copy files for a nightly run
for (( x=0; x < "${#RUN_CASE[@]}"; x++ )); do
  if [ "${EXIT_CODES[$x]}" != 0 ]; then
    echo "${RUN_CASE[$x]}"' failure'
        EXIT_STATUS=1
   fi
done

# If no cases failed, print an all good message
if [[ $EXIT_STATUS -eq 0 ]]; then
   echo "" 
   echo "All cases ran to completion."
fi

exit $EXIT_STATUS
