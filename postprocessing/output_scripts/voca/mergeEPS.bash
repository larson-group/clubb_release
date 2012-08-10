#!/bin/bash

#created by janhft 07/17/2012

myhelp()
{ 
echo "------------------------------------------------------------------------------------"
echo "This script combines all eps-files contained in the input directory to one pdf file." 
echo "The pdf file will be stored in the input directory."
echo ""
echo "Usage: ./mergeEPS.bash [-fp] <indir> <outfile>"
echo ""
echo "Options:"
echo "-h -- show this page"
echo "-p -- use png graphics"
echo "-f -- force overwrite if <outfile> already exists"
echo "------------------------------------------------------------------------------------"
}

# parse options
force=false
eps=true

while getopts pfh opt
do
  case "$opt" in
    f) force=true; shift;;
    p) eps=false; shift;;
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
echo '\documentclass{beamer}'>>tmp.tex
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
