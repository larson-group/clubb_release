#!/bin/bash

#created by janhft 07/17/2012

myhelp()
{ 
echo "------------------------------------------------------------------------------------"
echo "This script combines all eps-files contained in the input directory to one pdf file." 
echo "The pdf file will be stored in the input directory."
echo ""
echo "Usage: ./mergeEPS.bash [-fpd] <indir> <outfile>"
echo ""
echo "Options:"
echo "-h -- show this page"
echo "-p -- use png graphics"
echo "-f -- force overwrite if <outfile> already exists"
echo "-d -- create pdf for diurnal cycle"
echo "------------------------------------------------------------------------------------"
}

# parse options
force=false
eps=true
diurnal=false

while getopts pfdh opt
do
  case "$opt" in
    f) force=true; shift;;
    p) eps=false; shift;;
		d) diurnal=true; shift;;
    h) myhelp; exit;;
    \?) echo "Error: Unknown option $opt."; myhelp;;
  esac
done

# check arguments
indir=$1
outfile=$2
overwrite='no'

# make sure there is a '/' at the end of the path
indir=${indir%*/}'/'

if [ "$outfile" = "" ]; then
	echo "ERROR: No output file specified!"
	exit	
fi

if [ ! -d $indir ]; then
	echo "ERROR: Input directory does not exist!"
	exit
fi

if [ -f $indir''$outfile ] && [ $force = false ]; then
	echo "Outputfile already exists. Do you want to overwrite it? (yes/no)"
	read overwrite
	
	if [ "$overwrite" != "yes" ] && [ "$overwrite" != "y" ]; then
		echo "Operation canceled by user. No pdf file has been written."		
		exit
	fi
fi

# create normal pdf (for one set of plots)
if [ $diurnal = false ]; then

	#get eps graphics from indir
	k=0
	echo "Including: "
	for curFile in $indir*; do
		if [ $eps = true ] && [ "${curFile/*./}" = "eps" ]; then
			files[$k]=$curFile;
			echo $curFile 
			let k=k+1;
		elif [ $eps = false ] && [ "${curFile/*./}" = "png" ]; then
			files[$k]=$curFile;
			echo $curFile		
			let k=k+1;
		fi	
	done

	#write tex-file

	#header
	echo '\documentclass{beamer}'>tmp.tex
	echo '\usepackage{graphicx}'>>tmp.tex
	echo '\begin{document}'>>tmp.tex
	echo '\thispagestyle{empty}'>>tmp.tex

	#insert figures from eps-graphics 
	#figure 0
	echo '\begin{figure}'>>tmp.tex
	echo '\centering'>>tmp.tex
	echo '\fbox{'>>tmp.tex
	echo '	\includegraphics[width=0.9\textwidth]{'${files[0]}'}'>>tmp.tex
	echo '}'>>tmp.tex
	echo '\end{figure}'>>tmp.tex
	#echo '\clearpage'>>tmp.tex

	#remaining figures
	let n=k-1
	for i in `seq 1 $n`; do
		echo '\begin{figure}'>>tmp.tex
		echo '\centering'>>tmp.tex
		echo '\fbox{'>>tmp.tex
		echo '	\includegraphics[width=0.9\textwidth]{'${files[$i]}'}'>>tmp.tex
		echo '}'>>tmp.tex
		echo '\end{figure}'>>tmp.tex
		#echo '\clearpage'>>tmp.tex
	done
	echo '\end{document}'>>tmp.tex

# create pdf for diurnal cycle comparisson (3 sets of plots)
else

	dir_daytime=plots_daytime/
	dir_nighttime=plots_nighttime/
	dir_fulltime=plots_fulltime/

	k=0
	echo "Including: "
	for curFile in $indir$dir_daytime*; do
		if [ $eps = true ] && [ "${curFile/*./}" = "eps" ]; then
			files_dt[$k]=$curFile;
			echo $curFile 
			let k=k+1;
		elif [ $eps = false ] && [ "${curFile/*./}" = "png" ]; then
			files_dt[$k]=$curFile;
			echo $curFile		
			let k=k+1;
		fi	
	done

	k=0
	echo "Including: "
	for curFile in $indir$dir_nighttime*; do
		if [ $eps = true ] && [ "${curFile/*./}" = "eps" ]; then
			files_nt[$k]=$curFile;
			echo $curFile 
			let k=k+1;
		elif [ $eps = false ] && [ "${curFile/*./}" = "png" ]; then
			files_nt[$k]=$curFile;
			echo $curFile		
			let k=k+1;
		fi	
	done

	k=0
	echo "Including: "
	for curFile in $indir$dir_fulltime*; do
		if [ $eps = true ] && [ "${curFile/*./}" = "eps" ]; then
			files_ft[$k]=$curFile;
			echo $curFile 
			let k=k+1;
		elif [ $eps = false ] && [ "${curFile/*./}" = "png" ]; then
			files_ft[$k]=$curFile;
			echo $curFile		
			let k=k+1;
		fi	
	done

	#create tex file
	#header
	echo '\documentclass{beamer}'>tmp.tex
	echo '\usepackage{graphicx}'>>tmp.tex
	echo '\begin{document}'>>tmp.tex
	echo '\thispagestyle{empty}'>>tmp.tex

	#insert figures from eps-graphics 
	#figure 0
	# daytime and nighttime in one figure
	echo '\begin{figure}'>>tmp.tex
	echo '\centering'>>tmp.tex
	echo '\begin{tabular}{c c}'>>tmp.tex
	echo '\fbox{'>>tmp.tex
	echo '	\includegraphics[width=0.45\textwidth]{'${files_dt[0]}'}'>>tmp.tex
	echo '}'>>tmp.tex
	echo '&'>>tmp.tex
	echo '\fbox{'>>tmp.tex
	echo '	\includegraphics[width=0.45\textwidth]{'${files_nt[0]}'}'>>tmp.tex
	echo '}'>>tmp.tex
	echo '\end{tabular}'>>tmp.tex
	echo '\end{figure}'>>tmp.tex
	#echo '\clearpage'>>tmp.tex

	# fulltime seperate
	echo '\begin{figure}'>>tmp.tex
	echo '\centering'>>tmp.tex
	echo '\fbox{'>>tmp.tex
	echo '	\includegraphics[width=0.9\textwidth]{'${files_ft[0]}'}'>>tmp.tex
	echo '}'>>tmp.tex
	echo '\end{figure}'>>tmp.tex
  
	#remaining figures
	let n=k-1
	for i in `seq 1 $n`; do

		# daytime and nighttime in one figure
		echo '\begin{figure}'>>tmp.tex
		echo '\centering'>>tmp.tex
	  echo '\begin{tabular}{c c}'>>tmp.tex
	  echo '\fbox{'>>tmp.tex
	  echo '	\includegraphics[width=0.45\textwidth]{'${files_dt[$i]}'}'>>tmp.tex
	  echo '}'>>tmp.tex
	  echo '&'>>tmp.tex
	  echo '\fbox{'>>tmp.tex
	  echo '	\includegraphics[width=0.45\textwidth]{'${files_nt[$i]}'}'>>tmp.tex
	  echo '}'>>tmp.tex
	  echo '\end{tabular}'>>tmp.tex
		echo '\end{figure}'>>tmp.tex
		#echo '\clearpage'>>tmp.tex

		# fulltime seperate
		echo '\begin{figure}'>>tmp.tex
		echo '\centering'>>tmp.tex
		echo '\fbox{'>>tmp.tex
		echo '	\includegraphics[width=0.9\textwidth]{'${files_ft[$i]}'}'>>tmp.tex
		echo '}'>>tmp.tex
		echo '\end{figure}'>>tmp.tex
	done
	echo '\end{document}'>>tmp.tex

fi

#convert tex->pdf
if [ $eps = true ]; then
	#convert
	latex tmp.tex >> /dev/null
	dvipdf tmp.dvi >> /dev/null

	#clean up
	rm tmp.tex tmp.aux tmp.dvi tmp.toc tmp.nav tmp.out tmp.snm tmp.log  
else
	#convert
	pdflatex tmp.tex >> /dev/null

	#clean up
	rm tmp.tex tmp.aux tmp.toc tmp.nav tmp.out tmp.snm tmp.log
fi

mv tmp.pdf $indir''$outfile

if [ -f $indir''$outfile ]; then
	echo "Operation successful. Output has been written to $outfile."
else
	echo "ERROR: Something went wrong while compiling the tex-document."
fi
