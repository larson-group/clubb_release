#!/bin/bash
###############################################################################
# make_latex_docs.bash
# Author: Ryan Senkbeil
#
# Converts latex documents to html and pdf.
###############################################################################

# Get the current location so it can be restored later
RestoreLoc=`pwd`

# Get the location of this script
# readlink -f follows any sym links, dirname gets the directory.
ScriptLoc=`readlink -f $0`
ScriptLoc=`dirname $ScriptLoc`

# Go to the script's location
cd $ScriptLoc

# For every file that ends with .tex
for file in *.tex ; do
	# Get the filename without the extension
	BaseFileName=`basename $file .tex`

	# Compile the file
	latex $file

	# Create a PostScript file
	dvips -Pcmz -o $BaseFileName.ps $BaseFileName

	# Create PDF file
	ps2pdf $BaseFileName.ps

	# Create HTML file
	# -dschanen changed this to use TtH rather than latex2html
	# 18 Aug 2008
#	latex2html $file
 	tth $file
	# Remove header from index.html
#	vim -E -s $BaseFileName/index.html <<-EOF
#  		:/<!--Navigation Panel-->/,/<!--End of Navigation Panel-->/d
#		:/<!--Navigation Panel-->/,/<!--End of Navigation Panel-->/d
#		:/<!--Table of Child-Links-->/,/<!--End of Table of Child-Links-->/d
#   		:update
#   		:quit
#	EOF

	# Move the index file to ../
#	mv $BaseFileName/index.html $BaseFileName.html
	
	# Clean up
	rm -rf $BaseFileName.aux $BaseFileName.dvi $BaseFileName.log $BaseFileName
done

cd $RestoreLoc
