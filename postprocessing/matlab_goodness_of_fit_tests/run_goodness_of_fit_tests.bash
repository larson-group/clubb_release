#!/bin/bash
# $Id$
#-------------------------------------------------------------------------------
# Description:
# Script to run the MATLAB goodness-of-fit tests from the command line.

# Figure out the directory where the script is located.
run_dir=`dirname $0`

# Change directories to the one in which the script is located.
cd $run_dir

# Check for files (and paths) passed in through the command line.  The first
# entry is the SAM LES 3D DIRECTORY (containing SAM LES 3D files in netCDF
# format).  Every subsequent entry is a CLUBB FILE (also in netCDF format).
# When at least two files aren't present in the command line, stop the program.
if [ -z $2 ]; then

   echo ""
   echo "At least two entries are necessary."
   echo ""
   exit

else

   # The first entry listed is always the SAM LES 3D output DIRECTORY.
   SAM_LES_3D_dir=$1

   # Loop through all files in the SAM 3D LES directory and enter them all into
   # the array of filenames, SAM_LES_3D_file.
   SAM_LES_3D_file_index=-1
   for file in $SAM_LES_3D_dir*
   do
      SAM_LES_3D_file_index=$((SAM_LES_3D_file_index+1))
      SAM_LES_3D_file[$((SAM_LES_3D_file_index))]=$file
   done

   # All the remaining entries listed are CLUBB FILES, which will be placed into
   # the CLUBB_file array.
   shift
   CLUBB_file_index=-1
   while (( "$#" )); do
      CLUBB_file_index=$((CLUBB_file_index+1))
      CLUBB_file[$((CLUBB_file_index))]=$1
      shift
   done

   if [ ${#SAM_LES_3D_file[@]} -gt 1 ]
   then
      for (( indx=1; indx<=${#SAM_LES_3D_file[@]}; indx++ ))
      do
         echo ""
         echo "SAM LES 3D file "$indx":  "${SAM_LES_3D_file[$((indx-1))]}
      done
   else
      echo ""
      echo "SAM LES 3D file:  "${SAM_LES_3D_file[0]}
   fi

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

# Number of SAM LES 3D output files.
num_SAM_LES_3D_files=${#SAM_LES_3D_file[@]}

# Number of CLUBB output files.
num_CLUBB_files=${#CLUBB_file[@]}

# Get rid of annoying MATLAB warning message about "No protocol specified."
export DISPLAY=$HOSTNAME:0

# Run the goodness-of-fit_tests. 
sudo -u matlabuser /usr/local/bin/matlab -nodisplay -nodesktop -nosplash -r "goodness_of_fit_tests( '${SAM_LES_3D_file[*]}', $num_SAM_LES_3D_files, '${CLUBB_file[*]}', $num_CLUBB_files ), exit"
