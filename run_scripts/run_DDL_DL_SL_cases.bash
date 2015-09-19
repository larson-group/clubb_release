#!/bin/bash
# $Id$
#-------------------------------------------------------------------------------
# Description:
# Script to run RICO, DYCOMS-II RF02 DS, and LBA using the DDL, DL, and SL
# hydrometeor PDFs.

# Figure out the directory where the script is located.
run_dir=`dirname $0`

# Change directories to the one in which the script is located.
cd $run_dir

# Clean the source code.
../compile/clean_all.bash

# Compile the source code.
../compile/compile.bash

# Current date in YYYYMMDD_hh:mm:ss.
date_string=`date +%Y%m%d_%H:%M:%S`

# Make a new directory named after the date string (this can easily be renamed
# later).
mkdir ../output/$date_string

# Make a new directory to store the DDL results.
mkdir ../output/$date_string/DDL

# Make a new directory to store the DL results.
mkdir ../output/$date_string/DL

# Make a new directory to store the SL results.
mkdir ../output/$date_string/SL

# Run the DDL cases.

# Run the RICO case.
./run_scm.bash -s ../input/stats/all_stats.in rico

# Run the DYCOMS-II RF02 DS case.
./run_scm.bash -s ../input/stats/all_stats.in dycoms2_rf02_ds

# Run the LBA case.
./run_scm.bash -s ../input/stats/all_stats.in lba

# Move the results to the DDL directory.
mv ../output/*.nc ../output/$date_string/DDL/.

# Setup the DL cases.
cp ../input/tunable_parameters/tunable_parameters.in \
   ../input/tunable_parameters/tunable_parameters.in.prev
sed "s/omicron        = 0.5/omicron        = 1.0/g" \
     ../input/tunable_parameters/tunable_parameters.in.prev \
     > ../input/tunable_parameters/tunable_parameters.in
rm ../input/tunable_parameters/tunable_parameters.in.prev

# Run the DL cases.

# Run the RICO case.
./run_scm.bash -s ../input/stats/all_stats.in rico

# Run the DYCOMS-II RF02 DS case.
./run_scm.bash -s ../input/stats/all_stats.in dycoms2_rf02_ds

# Run the LBA case.
./run_scm.bash -s ../input/stats/all_stats.in lba

# Move the results to the DL directory.
mv ../output/*.nc ../output/$date_string/DL/.

# Setup the SL cases.
cp ../src/CLUBB_core/model_flags.F90 ../src/CLUBB_core/model_flags.F90.prev
sed "s/l_use_precip_frac = .true./l_use_precip_frac = .false./g" \
     ../src/CLUBB_core/model_flags.F90.prev > ../src/CLUBB_core/model_flags.F90
rm ../src/CLUBB_core/model_flags.F90.prev
../compile/compile.bash

# Run the SL cases.

# Run the RICO case.
./run_scm.bash -s ../input/stats/all_stats.in rico

# Run the DYCOMS-II RF02 DS case.
./run_scm.bash -s ../input/stats/all_stats.in dycoms2_rf02_ds

# Run the LBA case.
./run_scm.bash -s ../input/stats/all_stats.in lba

# Move the results to the SL directory.
mv ../output/*.nc ../output/$date_string/SL/.

# Revert the changes made for the DL and SL cases in order to restore the souce
# code to its original state.
cp ../input/tunable_parameters/tunable_parameters.in \
   ../input/tunable_parameters/tunable_parameters.in.prev
sed "s/omicron        = 1.0/omicron        = 0.5/g" \
     ../input/tunable_parameters/tunable_parameters.in.prev \
     > ../input/tunable_parameters/tunable_parameters.in
rm ../input/tunable_parameters/tunable_parameters.in.prev
cp ../src/CLUBB_core/model_flags.F90 ../src/CLUBB_core/model_flags.F90.prev
sed "s/l_use_precip_frac = .false./l_use_precip_frac = .true./g" \
     ../src/CLUBB_core/model_flags.F90.prev > ../src/CLUBB_core/model_flags.F90
rm ../src/CLUBB_core/model_flags.F90.prev
