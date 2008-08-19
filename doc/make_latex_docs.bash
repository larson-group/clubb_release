#!/bin/bash
###############################################################################
# make_latex_docs.bash
# Author: Ryan Senkbeil Aug 2008
#
# Converts latex documents to pdf.
###############################################################################

# Get the current location so it can be restored later
RestoreLoc=`pwd`

# Get the location of this script
# readlink -f follows any sym links, dirname gets the directory.
ScriptLoc=`readlink -f $0`
ScriptLoc=`dirname $ScriptLoc`

# Go to the script's location
cd $ScriptLoc

# Clean up anything left from last time
make clean

# Make the documentation
make

# For every file that ends with .tex, create html page
for file in *.tex ; do
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
done

# Clean up
rm *.aux *.auxbbl.make *.auxdvi.make *.aux.make *.dvi *.d *.fls *.log *.ps

cd $RestoreLoc
