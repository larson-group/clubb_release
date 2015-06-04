#!/bin/bash
# $Id$
#-------------------------------------------------------------------------------
# Description:
# Script to run the MATLAB scatter/contour plots from the command line, rather
# than having to start up MATLAB.  This is the lazy person's easy plotter.

# Figure out the directory where the script is located.
run_dir=`dirname $0`

# Change directories to the one in which the script is located.
cd $run_dir

# Change the permission on the "output" directory so that it is writable by
# all users. 
chmod 777 output

# Check for files (and paths) passed in through the command line.  The first one
# is the SAM LES 3D file (in netCDF format).  Every subsequent file is a CLUBB
# file (also in netCDF format).  When at least two files aren't present in the
# command line, use the value of 'default' for both SAM and CLUBB (the CLUBB
# input file array will have a size of 1).  The default setting is found in
# PDF_scatter_contour_plotter.m.
if [ -z $2 ]; then

   echo ""
   echo "The 'default' files for SAM LES 3D output and CLUBB output" \
        "are being used."
   echo ""

   SAM_LES_3D_file='default'
   CLUBB_file='default'

else

   # The first file listed is always the SAM LES 3D output file.
   SAM_LES_3D_file=$1

   # All the remaining files listed are CLUBB files, which will be placed into
   # the CLUBB_file array.
   shift
   CLUBB_file_index=-1
   while (( "$#" )); do
      CLUBB_file_index=$((CLUBB_file_index+1))
      CLUBB_file[$((CLUBB_file_index))]=$1
      shift
   done

   echo ""
   echo "SAM LES 3D file:  "$SAM_LES_3D_file

   if [ ${#CLUBB_file[@]} -gt 1 ]
   then
      for (( indx=1; indx<=${#CLUBB_file[@]}; indx++ ))
      do
         echo ""
         echo "CLUBB file "$indx":  "${CLUBB_file[$((indx-1))]}
      done
   else
      echo ""
      echo "CLUBB file:  "${CLUBB_file[0]}
   fi

   echo ""

fi

# Number of CLUBB output files.
num_CLUBB_files=${#CLUBB_file[@]}

# Get rid of annoying MATLAB warning message about "No protocol specified."
export DISPLAY=$HOSTNAME:0

# Run the scatter/contour plotter. 
sudo -u matlabuser /usr/local/bin/matlab -nodisplay -nodesktop -nosplash -r "PDF_scatter_contour_plotter( '$SAM_LES_3D_file', '${CLUBB_file[*]}', $num_CLUBB_files ), exit"
