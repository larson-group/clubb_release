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
OUTPUT_DIR="/home/`whoami`/nightly_tests/output"

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
TEMP=`getopt -o :nh --long nightly,help -n 'run_scm_all.bash' -- "$@"`

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
            echo -e "\t-h, --help\t\t\tPrints this help message"

            # Since the options for run_scm.bash are also valid, print those too:
            ./run_scm.bash --help | grep -v "Usage: run_scm.bash" | grep -v "\-\-help"

            exit 1 ;;
		--) shift ; break ;;
		*) echo "Something bad happened!" ; exit 1 ;;
	esac
done

# Declare arrays that will be used later on in this script.
declare -a RUN_CASE
declare -a EXIT_CODES

a=0
while read line
do
    # If the line is not commented out (does not start with '!')
    if [[ $line != !* ]] && [[ ! -z $line ]] ; then
        RUN_CASE[$a]=$line
        a=$(($a+1));
    fi
done < "RUN_CASES"

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
else
    echo -e "\nPerforming standard run\n"
fi

# This will loop over all runs in sequence 
for (( x=0; x < "${#RUN_CASE[@]}"; x++ )); 
do
    echo -e "Running ${RUN_CASE[$x]}"
    RESULT=`./run_scm.bash $OPTIONS ${RUN_CASE[$x]} 2>&1`

    if [ $NIGHTLY == true ] ; then
        echo -e "$RESULT"
    fi

    RESULT_STATUS=`echo "$RESULT" | grep 'normally'`

    if [ -z "$RESULT_STATUS" ]; then
        EXIT_CODES[$x]=-1
       
        # If there was an error, and this is not running in nightly mode,
        # it will not be displayed. So, display the error here.
#        if [ $NIGHTLY != true ] ; then
#            echo -e "$RESULT"
#        fi
    fi
done

EXIT_STATUS=0

# Print the results and copy files for a nightly run
for (( x=0; x < "${#RUN_CASE[@]}"; x++ )); do
	if [ "${EXIT_CODES[$x]}" != 0 ]; then
		echo "${RUN_CASE[$x]}"' failure'
        EXIT_STATUS=1
 	fi
done

exit $EXIT_STATUS
