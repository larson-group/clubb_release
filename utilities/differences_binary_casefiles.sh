#/bin/bash
#
# This is a simple script for diffing GrADS files of the same name in two different directories.
# Written by: Joshua Fasching
#

# Default Directories

# Path to GrADS files generated from current repository source code.
BASELINE=~/test/hoc_v2.2_tuner/standalone/

# Path to GrADS files generated by user modified source code.
WORKING=~/hoc_v2.2_tuner/standalone/

# Flag used determine if multiple diff processes are to be used.
MULTIPROCESS=0

set_args()
{
         # Loop through the list of arguments ($1, $2...). This loop ignores
         # anything not starting with '-'.
         while [ -n "$(echo $1 | grep "-")" ]; do
                 case $1 in
                         # '--nightly' sets the script to run the nightly version
                         --working | -w ) shift
			 		  WORKING=$1;;
                         --baseline | -b ) shift
			 	           BASELINE=$1;;
			 -mp ) MULTIPROCESS=1;;
			 -h | --help )	echo "                                            "
			 	        echo "Usage: ./differences_binary_casefiles.sh (-w|-b|-mp)"
					echo "                                            "
					echo " This script performs diff on the GrADS control and data files in two directories."
					echo "                                            "
					echo " By default those directories are:          "
                                        echo "   Baseline => "$BASELINE
					echo "   Working  +> "$WORKING
					echo "                                            "
					echo " Directories can be specified with these command-line options:"
			     		echo "   -w|working      Specify working directory"
			     		echo "   -b|baseline     Specify baseline directory"
					echo "                                            "
					echo " This script can be run such that it spawms multiple diff processes."
					echo " To activate, use the -mp command-line argument."
					echo "                                            "
					echo " To do a multi-process run in which both directories are specified:"
					echo "   ./differences_binary_casefiles.sh -w hoc_v2.2_tuner/standalone/ -b prev/hoc_v2.2_tuner/standalone/ -mp "
					echo "                                            "
			     		exit;;
                 esac
                 # Shift moves the parameters up one. Ex: $2 -> $1 and so on.
                 # This is so we only have to check $1 on each iteration.
                 shift
         done
}
 
set_args $*

# If -mp is one of the command-line arguments. 
# Alert the user that the script has been initiated
# in multiprocess mode.

if [ $MULTIPROCESS -gt 0 ]; then
	echo "############################################################"
	echo "######### Multi-Process Mode. Terminate with CTRL-Z ########"
	echo "############################################################"
fi

cd $WORKING

# Diff all GrADS files in the working directory with their counterpart in the baseline directory.
for FILE in *.ctl
do
	if [ $MULTIPROCESS -gt 0 ]; then
		diff $WORKING$FILE $BASELINE$FILE &
	else
		diff $WORKING$FILE $BASELINE$FILE
	fi
done

for FILE in *.dat
do
	if [ $MULTIPROCESS -gt 0 ]; then
		diff $WORKING$FILE $BASELINE$FILE &
	else
		diff $WORKING$FILE $BASELINE$FILE
	fi
done

wait 
