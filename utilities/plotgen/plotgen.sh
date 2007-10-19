#!/bin/bash

#Provide help if the user asks for it
if [ "$1" == "" ]; then
	echo "plotgen: missing necessary arguments"
	echo "Try 'plotgen --help' for more information." 
	exit 1
elif [ "$1" == "--help" ]; then
	echo ""
        echo "Usage: 'plotgen Output_Arg HOC_sim1 [HOC_sim2] Output_Directory Plot_LES Plot_Best Plot_HOCDec17'"
	echo "Or:    'plotgen Output_Arg nightly Output_Directory'"
	echo ""
	echo "The first generates plots comparing the specified simulations, e.g." 
        echo "  plotgen -c /home/vlarson/HOC_dir1 /home/vlarson/HOC_dir2 /home/vlarson/HOC_output 1 0 0"
        echo ""
        echo "The latter generates the nightly plots, e.g."
        echo "  plotgen -r nightly /home/vlarson/HOC_output"
        echo ""
	echo "The results will be saved in 'Output_Directory' and"
        echo "  can be viewed with the following command:"
	echo "  'firefox Output_Directory/index.html'"
	echo ""
        echo "'Output_Arg' must be ONE of the following:"
	echo "  '-r' : Overwrite the files in the output directory"
	echo "  '-c' : Create the output directory before outputting"
	echo ""
        echo "'HOC_sim1' is the directory (including path) containing"
        echo "  GrADS output files for all cases you wish to plot."
        echo "All possible cases that plotgen recognizes"
        echo "  are listed at the end of generate_plots.sh." 
        echo "If a case is listed but is not contained in 'HOC_sim1',"
        echo "  then that case will not be plotted."
        echo ""
        echo "'HOC_sim2' is optional but follows the same rules"
        echo "  as 'HOC_sim1'."
        echo ""
        echo "If the user chooses the 'nightly' option, then"
        echo "  the new GrADS data to be plotted are assumed to be"
        echo "  in local sub-directories named" 
        echo "  HOC_current/ and HOC_previous/."
        echo ""
	echo "Plot_LES, Plot_Best, and Plot_HOCDec17 are all boolean," 
        echo " with 1 being plot and 0 being do not plot."
        echo "The 'nightly' option assumes all 3 flags are set to 1."
	exit
fi

#Remember where we came from
working_directory=`pwd`

#Move to a directory that its easier to work in
cd /home/matlabuser/plotgen

#These arguments are always in the same spot
output_arg="$1"
HOC_sim1="$2"

#See if we need to adjust the variables to take in to account directory relationships
HOC_sim1_dir=${HOC_sim1%/*}
HOC_sim1_dir_rel=${HOC_sim1_dir:0:1}

if [ "$HOC_sim1_dir_rel" != "/"  ]; then
	HOC_sim1="$working_directory/$HOC_sim1"
fi

echo $HOC_sim1

#Parse out the arguments
if [ "$2" == "nightly" ]; then
	if [ "$3" == "" ]; then
		echo "Please specify an output directory."
		cd $working_directory
		exit 1
	fi	
	HOC_sim1='HOC_previous'
	HOC_sim2='HOC_current'
	output_dir="$3"
	compare_LES=1
	compare_best=1
	compare_HOC=1
elif [ "$#" == 7 ]; then
	HOC_sim2="$3"
	output_dir="$4"
	compare_LES="$5"
	compare_best="$6"
	compare_HOC="$7"

	#See if we need to adjust these arguments for directory relationships
	HOC_sim2_dir=${HOC_sim2%/*}
	HOC_sim2_dir_rel=${HOC_sim2_dir:0:1}
	if [ "$HOC_sim2_dir_rel" != "/"  ]; then
		HOC_sim2="$working_directory/$HOC_sim2"
	fi

	ln -s $HOC_sim1
	ln -s $HOC_sim2
elif [ "$#" == 6 ]; then
	HOC_sim2=0
	output_dir="$3"
	compare_LES="$4"
	compare_best="$5"
	compare_HOC="$6"

	ln -s $HOC_sim1
else
	echo "plotgen: missing necessary arguments"
	echo "Try 'plotgen --help' for more information." 
	cd $working_directory
	exit 1
fi

#Strip the directory information from the sim name
HOC_sim1=${HOC_sim1##/*/}
HOC_sim2=${HOC_sim2##/*/}

echo $HOC_sim1
echo $HOC_sim2		

#Before copying the plots, lets see if the desired directory exists
output_dir_rel=${output_dir:0:1}
if [ "$output_dir_rel" != "/"  ]; then
	output_dir="$working_directory/$output_dir"
fi

if [ "$output_arg" == "-r"  ]; then
	sudo -u matlabuser /home/matlabuser/plotgen/generate_plots.sh $HOC_sim1 $HOC_sim2 $compare_LES $compare_best $compare_HOC && \
	rm -rf "$output_dir"
	mkdir "$output_dir"
	cp -rf /home/matlabuser/plotgen/profiles/* "$output_dir"
elif [ "$output_arg" == "-c" ]; then
	if [ -e "$output_dir" ]; then
		echo "Directory exists, use '-r' to overwrite." 
		cd $working_directory
		exit 1
	fi
	mkdir -p "$output_dir" && \
	sudo -u matlabuser /home/matlabuser/plotgen/generate_plots.sh $HOC_sim1 $HOC_sim2 $compare_LES $compare_best $compare_HOC && \
	cp -rf /home/matlabuser/plotgen/profiles/* "$output_dir"
else
	sudo -u matlabuser /home/matlabuser/plotgen/generate_plots.sh $HOC_sim1 $HOC_sim2 $compare_LES $compare_best $compare_HOC
	echo "Invalid output argument."
	echo "Results stored in /home/matlabuser/plotgen/profiles"
fi

#Remove the symlinks, but not if we're doing nightly plots
if [ "$2" != "nightly" ]; then
	rm -f $HOC_sim1 $HOC_sim2
fi

#Take us back to where we started
cd $working_directory
exit
