#!/bin/bash

# Script to compile Latex documentation file and generate Postscript and pdf outputs

for file in hoceqns
do

# Compile Latex and bibtex files

  latex ${file}
  bibtex ${file}
  latex ${file}
  latex ${file}

# Create PostScript file

  dvips -Pcmz -o ${file}.ps ${file}

# Create pdf file

  ps2pdf ${file}.ps

# Clean up

  rm ${file}.aux ${file}.bbl ${file}.blg ${file}.dvi ${file}.log

done

