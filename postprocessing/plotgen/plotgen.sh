#!/bin/bash

# senkbeir changed the PLOTGEN_DIR variable to automatically figure out where it is
# This variable holds the path to the directory where plotgen.sh is located
# The readlink -f is necessary so if running from the symlink, it gets the 
# full path to plotgen.sh, and dirname $0 gets the directory.
#PLOTGEN_DIR=$(readlink -f $(dirname "$0"))
PLOTGEN_DIR=`readlink -f $0`
PLOTGEN_DIR=`dirname $PLOTGEN_DIR`
#External Use: Comment above 2 lines, uncomment line below
#PLOTGEN_DIR=`pwd`

#We need to unset the term type for this to work with Matlab 2008a
SessionType=$DISPLAY
unset DISPLAY

#Catch any termination signals
trap `DISPLAY=$SessionType;export DISPLAY; exit` INT TERM

#Provide help if the user asks for it
if [ "$1" == "" ]; then
	echo "plotgen: missing necessary arguments"
	echo "Try 'plotgen --help' for more information." 
	exit 1
elif [ "$1" == "--help" ]; then
	echo ""
        echo "Usage: 'plotgen Output_Arg HOC_dir1 [HOC_dir2] Output_Directory Plot_LES Plot_Best Plot_HOCDec17'"
	echo ""
	echo "This generates plots comparing the specified simulations, e.g." 
        echo "  plotgen -c /home/vlarson/HOC_dir1 /home/vlarson/HOC_dir2 /home/vlarson/HOC_output 1 0 0"
        echo ""
	echo "The results will be saved in 'Output_Directory' and"
        echo "  can be viewed with the following command:"
	echo "  'firefox Output_Directory/index.html'"
	echo ""
        echo "'Output_Arg' must be ONE of the following:"
	echo "  '-r' : Overwrite the files in the output directory"
	echo "  '-c' : Create the output directory before outputting"
	echo ""
        echo "'HOC_dir1' is the directory (including path) containing"
        echo "  GrADS output files for all cases you wish to plot."
        echo "All possible cases that plotgen recognizes"
        echo "  are listed at the end of generate_plots.sh." 
        echo "If a case is listed but is not contained in 'HOC_dir1',"
        echo "  then that case will not be plotted."
        echo ""
        echo "'HOC_dir2' is optional but follows the same rules"
        echo "  as 'HOC_dir1'."
        echo ""
       	echo "Plot_LES, Plot_Best, and Plot_HOCDec17 are all boolean," 
        echo " with 1 being plot and 0 being do not plot."
       	exit
fi

#Remember where we came from
working_directory=`pwd`

#Move to a directory that its easier to work in
cd $PLOTGEN_DIR

#These arguments are always in the same spot
output_arg="$1"
HOC_dir1="$2"

#See if we need to adjust the variables to take in to account directory relationships
HOC_dir1_rel=${HOC_dir1:0:1}

if [ "$HOC_dir1_rel" != "/"  ]; then
	HOC_dir1="$working_directory/$HOC_dir1"
fi

#Parse out the arguments
if [ "$#" == 7 ]; then
	HOC_dir2="$3"
	output_dir="$4"
	compare_LES="$5"
	compare_best="$6"
	compare_HOC="$7"

	#See if we need to adjust these arguments for directory relationships
	HOC_dir2_rel=${HOC_dir2:0:1}
	if [ "$HOC_dir2_rel" != "/"  ]; then
		HOC_dir2="$working_directory/$HOC_dir2"
	fi

	ln -s $HOC_dir1
	ln -s $HOC_dir2
elif [ "$#" == 6 ]; then
	HOC_dir2=0
	output_dir="$3"
	compare_LES="$4"
	compare_best="$5"
	compare_HOC="$6"

	ln -s $HOC_dir1
else
	echo "plotgen: missing necessary arguments"
	echo "Try 'plotgen --help' for more information." 
	cd $working_directory
	exit 1
fi

# senkbeir added checks to see if data directories exist if needed
# If we are plotting LES, check for LES_files directory
if [ $compare_LES == 1 ]
then
	# If the "LES_files" directory does not exist
	if [ ! -d "LES_files" ]
	then
		# Try to create a symbolic link to the directory.
		# We will assume it is being run from the repository, but
		# still check to see if the directory exists first
		if [ -d "../../les_data" ]
		then
			# It probably is being run from the repository... Let's create a symbolic link
			ln -s ../../les_data LES_files
		else
			# We cannot find the LES data, exit
			echo "No LES data available to plot! Create the LES_files directory first,"
			echo "or don't plot the LES data."
			exit 1
		fi
	fi
fi

# If we are plotting the best ever, check for Chris_Golaz_best_ever directory
if [ $compare_best == 1 ]
then
	if [ ! -d "Chris_Golaz_best_ever" ]
	then
		# The directory doesn't exist, warn and exit
		echo "No best ever data available to plot! Create the Chris_Golaz_best_ever directory first,"
		echo "or don't plot the best ever data."
		exit 1
	fi
fi

# If we are plotting HOC Dec 17, check for HOC_20051217 directory
if [ $compare_HOC == 1 ]
then
	if [ ! -d "HOC_20051217" ]
	then
		# The directory doesn't exist, warn and exit
		echo "No HOC Dec 17 data available to plot! Create the HOC_20051217 directory first,"
		echo "or don't plot the HOC Dec 17 data"
		exit 1
	fi
fi
# end senkbeir's changes

#Strip the directory information from the sim name
HOC_dir1=${HOC_dir1##/*/}
HOC_dir2=${HOC_dir2##/*/}

echo $HOC_dir1
echo $HOC_dir2		

#Before copying the plots, lets see if the desired directory exists
output_dir_rel=${output_dir:0:1}
if [ "$output_dir_rel" != "/"  ]; then
	output_dir="$working_directory/$output_dir"
fi

if [ "$output_arg" == "-r"  ]; then
	#External Use: Remove "sudo - u matlabuser" from command below
	sudo -u matlabuser $PLOTGEN_DIR/generate_plots.sh $HOC_dir1 $HOC_dir2 $compare_LES $compare_best $compare_HOC && \
	rm -rf "$output_dir"
	mkdir "$output_dir"
	cp -rf $PLOTGEN_DIR/profiles/* "$output_dir"
elif [ "$output_arg" == "-c" ]; then
	if [ -e "$output_dir" ]; then
		echo "Directory exists, use '-r' to overwrite." 
		cd $working_directory
		exit 1
	fi
	mkdir -p "$output_dir" && \
	#External Use: Remove "sudo - u matlabuser" from command below
	sudo -u matlabuser $PLOTGEN_DIR/generate_plots.sh $HOC_dir1 $HOC_dir2 $compare_LES $compare_best $compare_HOC && \
	
	#External Use: Comment out the line immediately below, uncomment the line below it
	cp -rf $PLOTGEN_DIR/profiles/* "$output_dir"
	#cp -rf $PLOTGEN_DIR/output/* "$output_dir"
else
	#External Use: Remove "sudo - u matlabuser" from below
	sudo -u matlabuser $PLOTGEN_DIR/generate_plots.sh $HOC_dir1 $HOC_dir2 $compare_LES $compare_best $compare_HOC
	echo "Invalid output argument."
	echo "Results stored in $PLOTGEN_DIR/profiles"
fi

#Remove the symlinks, but not if we're doing nightly plots
if [ "$2" != "nightly" ]; then
	rm -f $HOC_dir1 $HOC_dir2
fi

#Take us back to where we started
cd $working_directory

#Reset the session type
DISPLAY=$SessionType
export DISPLAY

exit
