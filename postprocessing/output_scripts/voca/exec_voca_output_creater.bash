#!/bin/bash

#created by janhft 05/24/2012

myhelp()
{ 
echo "------------------------------------------------------------------------------------"
echo "This bash script was created to run the matlab script voca_output_creater.m."
echo "The Matlab script has to be located in the same directory as this script. To"
echo "run jobs on HD1 use the -H option. If you have large input files the -b "
echo "option might be very usefull. Use the -l <lat> option to set the latitude "
echo "<lat> for the cross section plots. The matlab script will be running in the "
echo "background then and you don't have to stay logged in. For more Information"
echo "have a look at ticket 47 of the WRF trac."
echo ""
echo "Note: If you are using this script to create plots, you can put additional"
echo "information on the plots by creating a file named run.nfo which contains one"
echo "line of text. This file has to be located in <indir>. You can also use the"
echo "-n option to assign a <name> to jobs. If the -m option is used <name> will"
echo "show up in the name of the pdf-output. If used on HD1, <name> will appear"
echo "in the bjobs table and in the job output filename. "
echo ""
echo "Attention: If you are using the -d option, all wrfout files have to be located"
echo "directly in <indir>. No other files should be in that directory except files"
echo "with the following extensions: .pdf, .png, .eps, .nfo."
echo ""
echo "Dependencies: matlab, [atd, bsub]"
echo "" 
echo "Usage: ./exec_voca_output_creater.bash [-bmH] [-n <name>] [-l <lat>] -a <action> -f <infile> | -d <indir>"
echo ""
echo "Options:"
echo "-h -- show this page"
echo "-f -- use input file <infile>"
echo "-d -- use input directory <indir>"
echo "-a -- perform action <action>"
echo "-l -- set latitude for cross section plots"
echo "-b -- run the script in the background"
echo "-m -- merge the eps-graphics to a pdf-file"
echo "-n -- Jobname <name>"
echo "-H -- run the script on HD1"
echo "------------------------------------------------------------------------------------"
}

# parse options
infile=''
indir=''
action=''
lat=''
background=false
merge=false
hd1=false
hd1name=MatlabJob

while getopts hHbmn:a:f:d:l: opt
do
  case "$opt" in
    f) infile="$OPTARG";;
    a) action="$OPTARG";;
    d) indir="$OPTARG";;
    l) lat="$OPTARG";;
    b) background=true;;
    m) merge=true;;
    H) hd1=true;;
    n) hd1name="$OPTARG";;
    h) myhelp; exit;;
    \?) echo "Error: Unknown option."; myhelp;;
  esac
done

# check parameters
if [ "$action" != "plot" ] && [ "$action" != "write" ] && [ "$action" != "keyboard" ]; then
	echo "Error: You have to specify an action (plot, write or keyboard). Use the -a option."
	exit
fi

if [ "$action" = "keyboard" ] && [ $background = true ]; then
	echo "Error: You can't use the -b option with the keyboard action."
	exit
fi

if [ "$infile" != "" ] && [ "$indir" != "" ]; then
	echo "Error: You cannot use the -f and -d option at the same time!"
	exit

elif [ "$infile" = "" ] && [ "$indir" = "" ]; then
	echo "Error: You have to specify a valid input file or an input directory! Use the -f or -d option respectively."
	exit
fi

if [ $background = true ] && [ $hd1 = true ]; then
	echo "Error: You cannot use the -b and -H option at the same time!"
	exit
fi

# if an existing input file was specified
if [ "$infile" != "" ] && [ -f $infile ]; then

	#run in background
	if [ $background = true ]; then

		if [ "$lat" != "" ]; then
			echo "matlab -r voca_output_creater\(\'$infile\'\,\'$action\'\,$lat\)" >tmp
		else
			echo "matlab -r voca_output_creater\(\'$infile\'\,\'$action\'\)" >tmp
		fi

		# merge eps to pdf
		if [ $merge = true ]; then	
			echo "./mergeEPS.bash -f ${infile%/*}/ ${infile%/*}/plots.pdf">>tmp
		fi 

		at now -f tmp
		rm tmp

	#run in foreground
	else
		# hd1 run
		if [ $hd1 = true ]; then
			if [ -f run_matlab_tmp.job ]; then
				rm run_matlab_tmp.job
			fi

			touch run_matlab_tmp.job

			echo '#BSUB -J '$hd1name>>run_matlab_tmp.job
			echo '#BSUB -o MatlabJob_'$hd1name'_%J'>>run_matlab_tmp.job
			echo '#BSUB -K'>>run_matlab_tmp.job

			if [ "$lat" != "" ]; then
				echo '/sharedapps/LS/matlab/bin/matlab -nodisplay -nosplash -nojvm -r '"voca_output_creater\(\'$infile\'\,\'$action\'\,$lat\)">>run_matlab_tmp.job
			else
				echo '/sharedapps/LS/matlab/bin/matlab -nodisplay -nosplash -nojvm -r '"voca_output_creater\(\'$infile\'\,\'$action\'\)">>run_matlab_tmp.job
			fi

			bsub < run_matlab_tmp.job

			while [ ! -e 'MatlabJob_'$hd1name'_'* ]; do
				sleep 1
			done

			mv 'MatlabJob_'$hd1name'_'* ${infile%/*}/

		# usual run		
		else
			if [ "$lat" != "" ]; then
				matlab -r voca_output_creater\(\'$infile\'\,\'$action\'\,$lat\)
			else
				matlab -r voca_output_creater\(\'$infile\'\,\'$action\'\)
			fi
		fi		
			
		# merge eps to pdf
		if [ $merge = true ]; then	
			./mergeEPS.bash -f ${infile%/*}/ plots_"$hd1name".pdf
		fi 
	fi	

# if an existing input directory was specified
elif [ "$indir" != "" ] && [ -d $indir ]; then
	
	k=0	
	# get files from indir; exclude all png and eps
	for curFile in $indir*; do
		if [ "${curFile/*./}" != "png" ] && [ "${curFile/*./}" != "eps" ] && [ "${curFile/*./}" != "pdf" ] && [ "${curFile/*./}" != "nfo" ]; then
			files[$k]=$curFile;
			let k=k+1;
		fi	
	done

	# create string for the matlab call 
	let n=k-1;	
	fileList="\{\'"${files[0]}"\'"
	for i in `seq 1 $n`; do
		 fileList=$fileList"\;\'"${files[$i]}"\'"
	done
	fileList=$fileList'\}'

	# run in background
	if [ $background = true ]; then

		if [ "$lat" != "" ]; then
			echo "matlab -r voca_output_creater\($fileList\,\'$action\'\,$lat\)" >tmp
		else
			echo "matlab -r voca_output_creater\($fileList\,\'$action\'\)" >tmp
		fi

		# merge eps to pdf
		if [ $merge = true ]; then	
			echo "./mergeEPS.bash -f ${infile%/*}/ ${infile%/*}/plots.pdf">>tmp
		fi
 
		at now -f tmp
		rm tmp

	# run in foreground
	else		
		# hd1 run
		if [ $hd1 = true ]; then

			if [ -f run_matlab_tmp.job ]; then
				rm run_matlab_tmp.job
			fi

			touch run_matlab.job

			echo '#BSUB -J '$hd1name>>run_matlab_tmp.job
			echo '#BSUB -o MatlabJob_'$hd1name'_output_%J'>>run_matlab_tmp.job
			echo '#BSUB -K'>>run_matlab_tmp.job

			if [ "$lat" != "" ]; then
				echo '/sharedapps/LS/matlab/bin/matlab -nodisplay -nosplash -nojvm -r '"voca_output_creater\($fileList\,\'$action\'\,$lat\)">>run_matlab_tmp.job
			else
				echo '/sharedapps/LS/matlab/bin/matlab -nodisplay -nosplash -nojvm -r '"voca_output_creater\($fileList\,\'$action\'\)">>run_matlab_tmp.job
			fi

			bsub < run_matlab_tmp.job

			while [ ! -e 'MatlabJob_'$hd1name'_'* ]; do
				sleep 1
			done

			echo "moving to ${indir%*/}/"
			mv 'MatlabJob_'$hd1name'_'* ${indir%*/}/

		# usual run
		else 
			if [ "$lat" != "" ]; then
				echo "matlab -r voca_output_creater\($fileList\,\'$action\'\,$lat\)"
				matlab -r voca_output_creater\($fileList\,\'$action\'\,$lat\)
			else
				echo "matlab -r voca_output_creater\($fileList\,\'$action\'\)"
				matlab -r voca_output_creater\($fileList\,\'$action\'\)
			fi
		fi		

		# merge eps to pdf
		if [ $merge = true ]; then	
			./mergeEPS.bash -f ${indir%*/}/ plots_"$hd1name".pdf
		fi 
	fi
	
else
	echo "Error: The specified file or directory does not exist!"
	exit
fi


