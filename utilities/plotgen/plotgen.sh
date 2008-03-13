#!/bin/bash

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
cd /home/matlabuser/plotgen

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
	sudo -u matlabuser /home/matlabuser/plotgen/generate_plots.sh $HOC_dir1 $HOC_dir2 $compare_LES $compare_best $compare_HOC && \
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
	sudo -u matlabuser /home/matlabuser/plotgen/generate_plots.sh $HOC_dir1 $HOC_dir2 $compare_LES $compare_best $compare_HOC && \
	cp -rf /home/matlabuser/plotgen/profiles/* "$output_dir"
else
	sudo -u matlabuser /home/matlabuser/plotgen/generate_plots.sh $HOC_dir1 $HOC_dir2 $compare_LES $compare_best $compare_HOC
	echo "Invalid output argument."
	echo "Results stored in /home/matlabuser/plotgen/profiles"
fi

#Remove the symlinks, but not if we're doing nightly plots
if [ "$2" != "nightly" ]; then
	rm -f $HOC_dir1 $HOC_dir2
fi

#Take us back to where we started
cd $working_directory
exit
