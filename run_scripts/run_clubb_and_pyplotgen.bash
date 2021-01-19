#!/bin/bash
#####################################################################
#
# Script that runs one or more CLUBB cases (BOMEX, RICO, etc.), 
# outputs the results to directory output/$NEW_OUTPUT_DIR,
# plots the results using pyplotgen, and stores the plots in  
# directory output/pyplots_$NEW_OUTPUT_DIR.
#
#####################################################################

# These are all the cloud cases that will be run.
CASES=( arm_97 arm bomex rico cgils_s6 cgils_s12 dycoms2_rf01 dycoms2_rf02_ds lba gabls3_night )

# This is a subdirectory of directory output 
#    that will be created and will store new CLUBB output.
NEW_OUTPUT_DIR="new"

# These are two subdirectories of directory output that contain 
#    existing CLUBB output that we'd like to include in our plots for comparison.
#    Enter OLD_OUTPUT_DIR_1="" and OLD_OUTPUT_DIR_2="" if you don't have 
#    either OLD_OUTPUT_DIR_1 or OLD_OUTPUT_DIR_2.  
OLD_OUTPUT_DIR_1="default"
OLD_OUTPUT_DIR_2="" # (if you have no 2nd output directory)

# Plotting options for pyplotgen (see postprocessing/pyplotgen/README.md)
PYPLOTGEN_OPTIONS="--plot-budgets -l --cases ${CASES[@]} -c"

# Code begins ------------------------------------------------------

# Figure out the directory where the script is located
scriptPath=`dirname $0`

# Store the current directory location so it can be restored
restoreDir=`pwd`

# Change directories to the one that this script is located in,
#   assumed to be directory run_scripts
cd $scriptPath

# Make sure that we won't overwrite existing output.
if [ -d "../output/$NEW_OUTPUT_DIR" ]; then
   printf "Error: directory output/$NEW_OUTPUT_DIR already exists.  Aborting.\n"
   exit
fi
if [ -d "../output/pyplots_$NEW_OUTPUT_DIR" ]; then
   printf "Error: directory output/pyplots_$NEW_OUTPUT_DIR already exists.  Aborting.\n"
   exit
fi
mkdir ../output/$NEW_OUTPUT_DIR

# Compile CLUBB
../compile/compile.bash

# Run all the cases specified in variable CASES.
for (( j=0; j < ${#CASES[@]}; j++ )) do
   printf "Running case ${CASES[$j]} . . . \n"
   ./run_scm.bash ${CASES[$j]} &>/dev/null 
done

printf "\nAll cases have finished running.\n"

# Move output to subdirectory $NEW_OUTPUT_DIR of directory output
for (( j=0; j < ${#CASES[@]}; j++ )) do
   mv ../output/${CASES[$j]}*.nc ../output/$NEW_OUTPUT_DIR
   mv ../output/${CASES[$j]}*.txt ../output/$NEW_OUTPUT_DIR
done

# Create plots using pyplotgen.
printf "\nCreating plots using pyplotgen . . .\n\n"
if [[ -z "$OLD_OUTPUT_DIR_1" && -z "$OLD_OUTPUT_DIR_2" ]]
then
python3 ../postprocessing/pyplotgen/pyplotgen.py $PYPLOTGEN_OPTIONS ../output/$NEW_OUTPUT_DIR -o ../output/pyplots_$NEW_OUTPUT_DIR
elif [ -z "$OLD_OUTPUT_DIR_1" ]
then
python3 ../postprocessing/pyplotgen/pyplotgen.py $PYPLOTGEN_OPTIONS ../output/$OLD_OUTPUT_DIR_2 ../output/$NEW_OUTPUT_DIR -o ../output/pyplots_$NEW_OUTPUT_DIR
elif [ -z "$OLD_OUTPUT_DIR_2" ]
then
python3 ../postprocessing/pyplotgen/pyplotgen.py $PYPLOTGEN_OPTIONS ../output/$OLD_OUTPUT_DIR_1 ../output/$NEW_OUTPUT_DIR -o ../output/pyplots_$NEW_OUTPUT_DIR
else
python3 ../postprocessing/pyplotgen/pyplotgen.py $PYPLOTGEN_OPTIONS ../output/$OLD_OUTPUT_DIR_1 ../output/$OLD_OUTPUT_DIR_2 ../output/$NEW_OUTPUT_DIR -o ../output/pyplots_$NEW_OUTPUT_DIR
fi

exit
