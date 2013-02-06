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
merge=false

while getopts Hd:n:m: opt
do
  case "$opt" in
    H) hd1=true; shift;;
    d) daytime="$OPTARG";;
    n) nighttime="$OPTARG";;
    m) merge=true; pdf_name="$OPTARG";;
    h) myhelp; exit;;
    \?) echo "Error: Unknown option $opt."; myhelp;;
  esac
done

# check arguments
indir=$1

# running on HD1
if [ $hd1 = true ]; then

  # daytime run
  ./exec_voca_output_creater.bash -H -d $indir -a plot -t $daytime -o daytime
	mkdir $indir/plots_daytime 
	mv $indir/daytime_* $indir/plots_daytime 

	# nighttime run
  ./exec_voca_output_creater.bash -H -d $indir -a plot -t $nighttime -o nighttime
	mkdir $indir/plots_nighttime
	mv $indir/nighttime_* $indir/plots_nighttime

	# overall run
  ./exec_voca_output_creater.bash -H -d $indir -a plot -o fulltime
	mkdir $indir/plots_fulltime
	mv $indir/fulltime_* $indir/plots_fulltime

# running elsewhere
else

  # daytime run
  ./exec_voca_output_creater.bash -d $indir -a plot -t $daytime -o daytime
	mkdir $indir/plots_daytime 
	mv $indir/daytime_* $indir/plots_daytime 

	# nighttime run
  ./exec_voca_output_creater.bash -d $indir -a plot -t $nighttime -o nighttime
	mkdir $indir/plots_nighttime
	mv $indir/nighttime_* $indir/plots_nighttime

	# overall run
  ./exec_voca_output_creater.bash -d $indir -a plot -o fulltime
	mkdir $indir/plots_fulltime
	mv $indir/fulltime_* $indir/plots_fulltime

fi

if [ $merge = true ]; then
	
	sleep 10

	./mergeEPS.bash -d $indir $pdf_name

fi

