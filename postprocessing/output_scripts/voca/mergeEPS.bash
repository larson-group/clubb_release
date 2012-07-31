#!/bin/bash

#created by janhft 07/17/2012

myhelp()
{ 
echo "------------------------------------------------------------------------------------"
echo "This script combines all eps-files contained in the input directory to one pdf file." 
echo "The pdf file will be stored in the input directory."
echo ""
echo "Usage: ./mergeEPS.bash [-f|h] <indir> <outfile>"
echo ""
echo "Options:"
echo "-h -- show this page"
echo "-f -- force overwrite if <outfile> already exists"
echo "------------------------------------------------------------------------------------"
}

# parse options
force=false

while getopts fh opt
do
  case "$opt" in
    f) force=true; shift;;
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
for curFile in $indir*; do
	if [ "${curFile/*./}" = "eps" ]; then
		files[$k]=$curFile;
		let k=k+1;
	fi	
done

#write tex-file
if [ -d mergeEPStmp ]; then
	rm -r mergeEPStmp
fi
mkdir mergeEPStmp
touch mergeEPStmp/tmp.tex

#header
echo '\documentclass{beamer}'>>mergeEPStmp/tmp.tex
echo '\usepackage{graphicx}'>>mergeEPStmp/tmp.tex
echo '\begin{document}'>>mergeEPStmp/tmp.tex
echo '\thispagestyle{empty}'>>mergeEPStmp/tmp.tex

#insert figures from eps-graphics 
#figure 0
echo '\begin{figure}'>>mergeEPStmp/tmp.tex
echo '\centering'>>mergeEPStmp/tmp.tex
echo '\fbox{'>>mergeEPStmp/tmp.tex
echo '	\includegraphics[width=\textwidth]{'${files[0]%.*}'}'>>mergeEPStmp/tmp.tex
echo '}'>>mergeEPStmp/tmp.tex
echo '\end{figure}'>>mergeEPStmp/tmp.tex
echo '\clearpage'>>mergeEPStmp/tmp.tex

#remaining figures
let n=k-1
for i in `seq 1 $n`; do
	echo '\begin{figure}'>>mergeEPStmp/tmp.tex
	echo '\centering'>>mergeEPStmp/tmp.tex
	echo '\fbox{'>>mergeEPStmp/tmp.tex
	echo '	\includegraphics[width=\textwidth]{'${files[$i]%.*}'}'>>mergeEPStmp/tmp.tex
	echo '}'>>mergeEPStmp/tmp.tex
	echo '\end{figure}'>>mergeEPStmp/tmp.tex
	echo '\clearpage'>>mergeEPStmp/tmp.tex
done
echo '\end{document}'>>mergeEPStmp/tmp.tex

#convert tex->pdf
cd mergeEPStmp
pdflatex tmp.tex >> /dev/null
mv tmp.pdf $indir''$outfile
cd ..

#clean up
rm -r mergeEPStmp

if [ -f $indir''$outfile ]; then
	echo "Operation successful. Output has been written to $outfile."
else
	echo "ERROR: Something went wrong while compiling the tex-document."
fi
