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
echo "If you are using the -t option, you have to use the following format for"
echo "the <timeinterval>: hh:mm:ss-hh:mm:ss."
echo ""
echo "Attention: If you are using the -d option, all wrfout files have to be located"
echo "directly in <indir>. No other files should be in that directory except files"
echo "with the following extensions: .pdf, .png, .eps, .nfo."
echo ""
echo "Dependencies: matlab, [atd, bsub]"
echo "" 
echo "Usage: ./exec_voca_output_creater.bash [-bmH] [-n <name>] [-l <lat>] [-t <timeinterval>] [-o <outfile_prefix>] -a <action> -f <infile> | -d <indir>"
echo ""
echo "Options:"
echo "-h -- show this page"
echo "-f -- use input file <infile>"
echo "-d -- use input directory <indir>"
echo "-o -- set output file prefix"
echo "-a -- perform action <action>"
echo "-t -- average over <timeinterval>"
echo "-l -- set latitude for cross section plots"
echo "-b -- run the script in the background"
echo "-m -- merge the eps-graphics to a pdf-file"
echo "-n -- job name of the hd1 job  <name>"
echo "-H -- run the script on HD1"
echo "-Y -- run the script on Yellowstone"
echo "------------------------------------------------------------------------------------"
}

# default values
lat_default=20

# parse options
infile=''
indir=''
action=''
outfile_prefix=''
lat=''
background=false
merge=false
hd1=false
yellowstone=false
hd1name=MatlabJob
timestr=''
l_valid_time=false

while getopts hHYbmn:a:f:d:o:l:t: opt
do
  case "$opt" in
    f) infile="$OPTARG";;
    a) action="$OPTARG";;
    d) indir="$OPTARG";;
    o) outfile_prefix="$OPTARG";;
    l) lat="$OPTARG";;
    b) background=true;;
    m) merge=true;;
    H) hd1=true;;
    Y) yellowstone=true;;
    n) hd1name="$OPTARG";;
    t) timestr="$OPTARG";;
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

# check timestr
if [ "$timestr" != "" ]; then
    
  tstart=${timestr:0:8}
  tend=${timestr:9:17}

  if [[ ! $tstart =~ [0-1][0-9]:[0-5][0-9]:[0-5][0-9]|2[0-3]:[0-5][0-9]:[0-5][0-9] ]] || [[ ! $tend =~ [0-1][0-9]:[0-5][0-9]:[0-5][0-9]|2[0-3]:[0-5][0-9]:[0-5][0-9] ]]; then 
    echo "Error: The entered time interval is invalid. It has to be of the format \"hh:mm:ss-hh:mm:ss\"."    
    exit
  else
    l_valid_time=true
  fi 

fi

# determine input files
if [ "$infile" != "" ] && [ -f $infile ]; then

  fileList="\'$infile\'"
  path=${infile%/*}"/"

	# set outfile_prefix if not specified
  if [ "$outfile_prefix" = "" ]; then
		outfile_prefix=${infile/*\//}
	fi

elif [ "$indir" != "" ] && [ -d $indir ]; then

  # make sure indir has a '/' at the end
	indir=${indir/%\//}/

  path=$indir

	k=0	
	# get files from indir
	for curFile in $indir*; do
		
		# substring to check if the file is a MatlabJob oputput file
		tmp=${curFile/*\//}

		# exclude all non-wrfout files
		if [[ $tmp == wrfout_*_[0-9][0-9]:[0-9][0-9]:[0-9][0-9] ]]; then
			files[$k]=$curFile;
			
			# set outfile_prefix if not specified
      if [ $k = 0 ] && [ "$outfile_prefix" = "" ]; then
				outfile_prefix=$tmp
			fi

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

else
	echo "Error: The specified file or directory does not exist!"
	exit
fi

# generate argument list for matlab call
arg_list="$fileList\,\'$action\'\,\'$outfile_prefix\'"

if [ "$lat" != "" ]; then
	arg_list=$arg_list"\,$lat"
else
	arg_list=$arg_list"\,$lat_default"
fi

if [ $l_valid_time = true ]; then
  arg_list=$arg_list"\,\'$tstart\'\,\'$tend\'"
fi

# matlab call (either background or foreground)
if [ $background = true ]; then # background run

  echo "matlab -r voca_output_creater\($arg_list\)" >tmp

	# merge eps to pdf
	if [ $merge = true ]; then	
		echo "./mergeEPS.bash -f $path "$path"plots.pdf">>tmp
	fi 

	at now -f tmp
	rm tmp

else # foreground run

	jobfilename='run_'$hd1name'_tmp.job'

	if [ $hd1 = true ]; then # hd1 run
		if [ -f $jobfilename ]; then
			rm $jobfilename
		fi

		touch $jobfilename

		echo '#BSUB -J '$hd1name>>$jobfilename
		echo '#BSUB -o MatlabJob_'$hd1name'_%J'>>$jobfilename
		echo '#BSUB -K'>>$jobfilename

		echo '/sharedapps/LS/matlab/bin/matlab -nodisplay -nosplash -nojvm -r '"voca_output_creater\($arg_list\)">>$jobfilename

		bsub < $jobfilename

		#while [ ! -e 'MatlabJob_'$hd1name'_'* ]; do
		#	sleep 1
		#done

		#echo "Moving output to $path"
  	#mv 'MatlabJob_'$hd1name'_'* $path

	# usual run
        elif [ $yellowstone = true ]; then # yellowstone run

                echo '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                echo '!!! Do not forget to load the Matlab module !!!'
                echo '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'

		if [ -f $jobfilename ]; then
			rm $jobfilename
		fi

		touch $jobfilename

                echo '#BSUB -J '$hd1name>>$jobfilename
                echo '#BSUB -o MatlabJob_'$hd1name'_%J'>>$jobfilename
                echo '#BSUB -K'>>$jobfilename
                echo '#BSUB -P P36741010'>>$jobfilename
                echo '#BSUB -n 1'>>$jobfilename
                echo '#BSUB -W 1:00'>>$jobfilename
                echo '#BSUB -q regular'>>$jobfilename

		echo 'matlab -nodisplay -nosplash -nojvm -r '"voca_output_creater\($arg_list\)">>$jobfilename

		bsub < $jobfilename
	else
			matlab -r "voca_output_creater\($arg_list\)"
      echo "If matlab did not run the script, use this command: "
      echo "matlab -r voca_output_creater\($arg_list\)"
  fi	

	# merge eps to pdf
	if [ $merge = true ]; then	
		./mergeEPS.bash -f $path plots_"$hd1name".pdf
	fi 
	
fi

exit
