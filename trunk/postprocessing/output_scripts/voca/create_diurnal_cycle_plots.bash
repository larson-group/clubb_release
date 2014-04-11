#!/bin/bash

#created by janhft 02/05/2013

myhelp()
{ 
echo "-------------------------------------------------------------------------------------"
echo "This script uses the exec_voca_output_creater.bash to create plots for the diurnal   " 
echo "cycle (daytime, nighttime, overall time). It runs the exec script three times and    "
echo "stores the output of the runs in different directories (e.g. plots_daytime) in "
echo "<indir>. Using the -m option creates a pdf which contains the plots of all three runs"
echo "(using the mergePDF.bash). The pdf will be stored under <indir>/<pdf_name>.pdf"
echo ""
echo "Usage: ./create_diurnal_cycle_plots.bash [-m <pdf_name>] [-d <daytime>] [-n <nighttime>] <indir>"
echo ""
echo "Options:"
echo "-h -- show this help page"
echo "-m -- merge the plots to <pdf_name>.pdf"
echo "-d -- set daytime to the time interval <daytime>"
echo "-n -- set daytime to the time interval <nighttime>"
echo "-H -- run the script on HD1"
echo "-j -- name of the job"
echo "-i -- path to input directory"
echo ""
echo "Using the -d or -n option requires a the following format for the time intervals:   "
echo "hh:mm:ss-hh:mm:ss"
echo "------------------------------------------------------------------------------------"
}

# parse options
daytime='05:00:00-08:00:00'
nighttime='17:00:00-20:00:00'
pdf_name=''
hd1=false
hd1opt=''
merge=false
jobname="MatlabJob"
indir=''


while getopts Hd:n:m:j:i: opt
do
  echo $opt
  case "$opt" in
    H) hd1=true;;
    d) daytime="$OPTARG";;
    n) nighttime="$OPTARG";;
    m) merge=true; pdf_name="$OPTARG";;
    j) jobname="$OPTARG";;
    i) indir="$OPTARG";;
    h) myhelp; exit;;
    \?) echo "Error: Unknown option $opt."; myhelp;;
  esac
done

if [ "$indir" = "" ]; then
	echo "Error: You have to specify a valid input directory! Use the -i option."
	exit
fi

# running on HD1
if [ $hd1 = true ]; then
        hd1opt="-H"
fi

# daytime run
./exec_voca_output_creater.bash $hd1opt -d $indir -a plot -t $daytime -o daytime -n $jobname"_day"
	mkdir $indir/plots_daytime 
	mv $indir/daytime_* $indir/plots_daytime 

# nighttime run
./exec_voca_output_creater.bash $hd1opt -d $indir -a plot -t $nighttime -o nighttime -n $jobname"_night"
	mkdir $indir/plots_nighttime
	mv $indir/nighttime_* $indir/plots_nighttime

# overall run
./exec_voca_output_creater.bash $hd1opt -d $indir -a plot -o fulltime -n $jobname"_full"
	mkdir $indir/plots_fulltime
	mv $indir/fulltime_* $indir/plots_fulltime


if [ $merge = true ]; then
	
	sleep 10

	./mergeEPS.bash -d $indir $pdf_name

fi

