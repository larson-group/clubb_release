#!/bin/bash

#created by janhft 05/24/2012

myhelp()
{ 
echo "------------------------------------------------------------------------------------"
echo "This bash script was created to run the matlab script voca_output_creater.m."
echo "The Matlab script has to be located in the same directory as this script. If"
echo "you have large input files the -b option might be very usefull. The matlab"
echo "script will be running in the background then and you don't have to stay"
echo "logged in. For more Information have a look at ticket 47 of the WRF trac."
echo "" 
echo "Usage: ./exec_voca_output_creater.bash [-bm] -a <action> -f <infile> | -d <indir>"
echo ""
echo "Options:"
echo "-h -- show this page"
echo "-f -- use input file <infile>"
echo "-d -- use input directory <indir>"
echo "-a -- perform action <action>"
echo "-b -- run the script in the background"
echo "-m -- merge the eps-graphics to a pdf-file"
echo "------------------------------------------------------------------------------------"
}

# parse options
infile=''
indir=''
action=''
background=false
merge=false

while getopts hbma:f:d: opt
do
  case "$opt" in
    f) infile="$OPTARG";;
    a) action="$OPTARG";;
    d) indir="$OPTARG";;
    b) background=true;;
    m) merge=true;;
    h) myhelp; exit;;
    \?) echo "Error: Unknown option $opt."; myhelp;;
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

# if an existing input file was specified
if [ "$infile" != "" ] && [ -f $infile ]; then

	if [ $background = true ]; then
		echo "matlab -r voca_output_creater\(\'$infile\'\,\'$action\'\)" >tmp
		if [ $merge = true ]; then	
			echo "./mergeEPS.bash -f ${infile%/*}/ ${infile%/*}/plots.pdf">>tmp
		fi 
		at now -f tmp
		rm tmp
	else
		matlab -r voca_output_creater\(\'$infile\'\,\'$action\'\)
		if [ $merge = true ]; then	
			./mergeEPS.bash -f ${infile%/*}/ ${infile%/*}/plots.pdf
		fi 
	fi	

# if an existing input directory was specified
elif [ "$indir" != "" ] && [ -d $indir ]; then
	k=0	
	#get files from indir; exclude all png and eps
	for curFile in $indir*; do
		if [ "${curFile/*./}" != "png" ] && [ "${curFile/*./}" != "eps" ]; then
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

	if [ $background = true ]; then
		echo "matlab -r voca_output_creater\($fileList\,\'$action\'\)" >tmp
		if [ $merge = true ]; then	
			echo "./mergeEPS.bash -f ${infile%/*}/ ${infile%/*}/plots.pdf">>tmp
		fi 
		at now -f tmp
		rm tmp
	else		
		echo "matlab -r voca_output_creater\($fileList\,\'$action\'\)"
		matlab -r voca_output_creater\($fileList\,\'$action\'\)
		if [ $merge = true ]; then	
			./mergeEPS.bash -f ${infile%/*}/ ${infile%/*}/plots.pdf
		fi 
	fi
	
else
	echo "Error: The specified file or directory does not exist!"
	exit
fi


