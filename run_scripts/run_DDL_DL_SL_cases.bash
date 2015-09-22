#!/bin/bash
# $Id$
#-------------------------------------------------------------------------------
# Description:
# Script to run RICO, DYCOMS-II RF02 DS, and LBA using the DDL, DL, and SL
# hydrometeor PDFs.
# When the relevant options are turned on (passed in from the command line), the
# script will additionally produce the plots of the PDFs, run the
# goodness-of-fit tests, and/or produce profiles (using plotgen) of warm
# microphysics process rates.

# Initialize options to false.
PDF_PLOTS=false
FIT_TESTS=false
MICRO_PROFILES=false

# Set options to true based on what is read in from the command line.
while (( "$#" )); do
   case "$1" in
      -p|--pdf_plots) # Plot PDFs
         PDF_PLOTS=true
         shift;;
      -f|--fit_tests) # Run the goodness-of-fit tests
         FIT_TESTS=true
         shift;;
      -m|--micro_profiles) # Produce profiles of warm microphysics
         MICRO_PROFILES=true
         shift;;
      -h|--help) # Print the help message
         echo -e "Usage: run_DDL_DL_SL_cases.bash [OPTIONS]"
         echo -e "\t-p or --pdf_plots to produce the plots of the PDFs"
         echo -e "\t-f or --fit_tests to run the goodness-of-fit tests"
         echo -e "\t-m or --micro_profiles to produce profiles of warm microphysics process rates"
         exit 1 ;;
   esac
done

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

# Options

# Produce plots of the PDFs.

# Run the goodness-of-fit tests.

# Produce profiles of warm microphysics process rates.
if [ $MICRO_PROFILES == true ];
then

   # Change directory to the one that contains the CLUBB .case files.
   cd ../postprocessing/plotgen/cases/clubb

   # Make a copy of the original files.
   cp rico.case rico.case.default
   cp dycoms2_rf02_ds.case dycoms2_rf02_ds.case.default
   cp lba.case lba.case.default

   # Copy the special files for warm microphysics profiles into the case/clubb
   # directory for use in plotgen.
   cp warm_micro_fields/rico.case rico.case
   cp warm_micro_fields/dycoms2_rf02_ds.case dycoms2_rf02_ds.case
   cp warm_micro_fields/lba.case lba.case

   # Change directory back to the plotgen main directory.
   cd ../../

   # Run plotgen
   ./plotgen.pl -c -l -e -t ../../output/$date_string/SL \
                ../../output/$date_string/DL \
                ../../output/$date_string/DDL \
                ../../output/$date_string/plot_SL_DL_DDL_mc_profiles

   # Change directory to the one that contains the CLUBB .case files.
   cd cases/clubb

   # Move the original .case files back.
   mv rico.case.default rico.case
   mv dycoms2_rf02_ds.case.default dycoms2_rf02_ds.case
   mv lba.case.default lba.case

   # Change directories to the one where the output is located.
   cd ../../../../output/$date_string

   # Produce a .maff file from the plotgen .eps output.
   zip -r plot_SL_DL_DDL_mc_profiles.maff plot_SL_DL_DDL_mc_profiles/

   # Change directories back to the one where the run script is located.
   cd ../../run_scripts

fi # Produce profiles of warm microphysics process rates.
