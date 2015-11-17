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
if [ $PDF_PLOTS == true ];
then

   # Change directory to the matlab_scatter_contour_plots directory.
   cd ../postprocessing/matlab_scatter_contour_plots

   # Copy PDF_scatter_contour_plotter.m back one directory.
   cp PDF_scatter_contour_plotter.m ../PDF_scatter_contour_plotter.m

   # Setup the PDF plots for the RICO case (rr and Nr).

   # Copy PDF_scatter_contour_plotter.m to PDF_scatter_contour_plotter.m.prev.
   cp PDF_scatter_contour_plotter.m PDF_scatter_contour_plotter.m.prev

   # Edit PDF_scatter_contour_plotter.m.
   sed "s/clubb_height_idx =/%clubb_height_idx =/g; \
        s/print_note =/%print_note =/g; s/%%/%/g; \
        s/%clubb_height_idx = 11/clubb_height_idx = 11/g; \
        s/%print_note = ''/print_note = ''/g; \
        s/plot_rt_thl  = true/plot_rt_thl  = false/g; \
        s/plot_chi_eta = true/plot_chi_eta = false/g; \
        s/plot_w_rr    = true/plot_w_rr    = false/g; \
        s/plot_w_Nr    = true/plot_w_Nr    = false/g; \
        s/plot_chi_rr  = true/plot_chi_rr  = false/g; \
        s/plot_chi_Nr  = true/plot_chi_Nr  = false/g; \
        s/plot_eta_rr  = true/plot_eta_rr  = false/g; \
        s/plot_eta_Nr  = true/plot_eta_Nr  = false/g; \
        s/plot_rr_Nr   = true/plot_rr_Nr   = false/g; \
        s/plot_w       = true/plot_w       = false/g; \
        s/plot_rt      = true/plot_rt      = false/g; \
        s/plot_thl     = true/plot_thl     = false/g; \
        s/plot_chi     = true/plot_chi     = false/g; \
        s/plot_eta     = true/plot_eta     = false/g; \
        s/plot_ri      = true/plot_ri      = false/g; \
        s/plot_Ni      = true/plot_Ni      = false/g; \
        s/plot_ln_ri   = true/plot_ln_ri   = false/g; \
        s/plot_ln_Ni   = true/plot_ln_Ni   = false/g; \
        s/plot_rs      = true/plot_rs      = false/g; \
        s/plot_Ns      = true/plot_Ns      = false/g; \
        s/plot_ln_rs   = true/plot_ln_rs   = false/g; \
        s/plot_ln_Ns   = true/plot_ln_Ns   = false/g; \
        s/plot_rg      = true/plot_rg      = false/g; \
        s/plot_Ng      = true/plot_Ng      = false/g; \
        s/plot_ln_rg   = true/plot_ln_rg   = false/g; \
        s/plot_ln_Ng   = true/plot_ln_Ng   = false/g; \
        s/auto_legend_text = true/auto_legend_text = false/g; \
        s/1:9/1:3/g; s/1:8/1:2/g; s/CLUBB DDL/DDL/g; \
        s/CLUBB DL/DL/g; s/CLUBB SL/SL/g;" \
        PDF_scatter_contour_plotter.m.prev > PDF_scatter_contour_plotter.m

   # Remove PDF_scatter_contour_plotter.m.prev.
   rm PDF_scatter_contour_plotter.m.prev

   # Run the plotter for the RICO case.
   ./run_scatter_contour_plots.bash \
      ../../sam_benchmark_runs/RICO_256x256x100_drizzle/3D_output/RICO_256x256x100_drizzle_128_0000252000_micro.nc \
      ../../output/$date_string/DDL/rico_zt.nc \
      ../../output/$date_string/DL/rico_zt.nc \
      ../../output/$date_string/SL/rico_zt.nc

   # Setup the PDF plots for the DYCOMS-II RF02 DS case (rr and Nr).

   # Copy PDF_scatter_contour_plotter.m to PDF_scatter_contour_plotter.m.prev.
   cp PDF_scatter_contour_plotter.m PDF_scatter_contour_plotter.m.prev

   # Edit PDF_scatter_contour_plotter.m.
   sed "s/clubb_height_idx =/%clubb_height_idx =/g;
        s/casename =/%casename =/g; s/%%/%/g; \
        s/%clubb_height_idx = 28/clubb_height_idx = 28/g; \
        s/%casename = 'RF02'/casename = 'RF02'/g;" \
        PDF_scatter_contour_plotter.m.prev > PDF_scatter_contour_plotter.m

   # Remove PDF_scatter_contour_plotter.m.prev.
   rm PDF_scatter_contour_plotter.m.prev

   # Run the plotter for the DYCOMS-II RF02 DS case.
   ./run_scatter_contour_plots.bash \
      ../../sam_benchmark_runs/DYCOMS_RF02_128x128x96_dr_sed/3D_output/DYCOMS_RF02_128x128x96_dr_sed_128_0000039600_micro.nc \
      ../../output/$date_string/DDL/dycoms2_rf02_ds_zt.nc \
      ../../output/$date_string/DL/dycoms2_rf02_ds_zt.nc \
      ../../output/$date_string/SL/dycoms2_rf02_ds_zt.nc

   # Setup the PDF plots for the LBA case (rr and Nr).

   # Copy PDF_scatter_contour_plotter.m to PDF_scatter_contour_plotter.m.prev.
   cp PDF_scatter_contour_plotter.m PDF_scatter_contour_plotter.m.prev

   # Edit PDF_scatter_contour_plotter.m.
   sed "s/clubb_height_idx =/%clubb_height_idx =/g;
        s/casename =/%casename =/g; s/%%/%/g; \
        s/%clubb_height_idx = 25/clubb_height_idx = 25/g; \
        s/%casename = 'LBA'/casename = 'LBA'/g;" \
        PDF_scatter_contour_plotter.m.prev > PDF_scatter_contour_plotter.m

   # Remove PDF_scatter_contour_plotter.m.prev.
   rm PDF_scatter_contour_plotter.m.prev

   # Run the plotter for the LBA case.
   ./run_scatter_contour_plots.bash \
      ../../sam_benchmark_runs/LBA_r1663_128x128x128_1km_Morrison/3D_output/LBA_128kmx128kmx128_1km_Morrison_64_0000003300_micro.nc \
      ../../output/$date_string/DDL/lba_zt.nc \
      ../../output/$date_string/DL/lba_zt.nc \
      ../../output/$date_string/SL/lba_zt.nc

   # Setup the PDF plots for the LBA case (rs, Ns, rg, and Ng).

   # Copy PDF_scatter_contour_plotter.m to PDF_scatter_contour_plotter.m.prev.
   cp PDF_scatter_contour_plotter.m PDF_scatter_contour_plotter.m.prev

   # Edit PDF_scatter_contour_plotter.m.
   sed "s/clubb_height_idx =/%clubb_height_idx =/g; s/%%/%/g; \
        s/%clubb_height_idx = 44/clubb_height_idx = 44/g; \
        s/plot_rr      = true/plot_rr      = false/g; \
        s/plot_Nr      = true/plot_Nr      = false/g; \
        s/plot_ln_rr   = true/plot_ln_rr   = false/g; \
        s/plot_ln_Nr   = true/plot_ln_Nr   = false/g; \
        s/plot_rs      = false/plot_rs      = true/g; \
        s/plot_Ns      = false/plot_Ns      = true/g; \
        s/plot_ln_rs   = false/plot_ln_rs   = true/g; \
        s/plot_ln_Ns   = false/plot_ln_Ns   = true/g; \
        s/plot_rg      = false/plot_rg      = true/g; \
        s/plot_Ng      = false/plot_Ng      = true/g; \
        s/plot_ln_rg   = false/plot_ln_rg   = true/g; \
        s/plot_ln_Ng   = false/plot_ln_Ng   = true/g;" \
        PDF_scatter_contour_plotter.m.prev > PDF_scatter_contour_plotter.m

   # Remove PDF_scatter_contour_plotter.m.prev.
   rm PDF_scatter_contour_plotter.m.prev

   # Run the plotter for the LBA case.
   ./run_scatter_contour_plots.bash \
      ../../sam_benchmark_runs/LBA_r1663_128x128x128_1km_Morrison/3D_output/LBA_128kmx128kmx128_1km_Morrison_64_0000003600_micro.nc \
      ../../output/$date_string/DDL/lba_zt.nc \
      ../../output/$date_string/DL/lba_zt.nc \
      ../../output/$date_string/SL/lba_zt.nc

   # Setup the PDF plots for the LBA case (ri and Ni).

   # Copy PDF_scatter_contour_plotter.m to PDF_scatter_contour_plotter.m.prev.
   cp PDF_scatter_contour_plotter.m PDF_scatter_contour_plotter.m.prev

   # Edit PDF_scatter_contour_plotter.m.
   sed "s/clubb_height_idx =/%clubb_height_idx =/g; s/%%/%/g; \
        s/%clubb_height_idx = 60/clubb_height_idx = 60/g; \
        s/plot_rs      = true/plot_rs      = false/g; \
        s/plot_Ns      = true/plot_Ns      = false/g; \
        s/plot_ln_rs   = true/plot_ln_rs   = false/g; \
        s/plot_ln_Ns   = true/plot_ln_Ns   = false/g; \
        s/plot_rg      = true/plot_rg      = false/g; \
        s/plot_Ng      = true/plot_Ng      = false/g; \
        s/plot_ln_rg   = true/plot_ln_rg   = false/g; \
        s/plot_ln_Ng   = true/plot_ln_Ng   = false/g; \
        s/plot_ri      = false/plot_ri      = true/g; \
        s/plot_Ni      = false/plot_Ni      = true/g; \
        s/plot_ln_ri   = false/plot_ln_ri   = true/g; \
        s/plot_ln_Ni   = false/plot_ln_Ni   = true/g;" \
        PDF_scatter_contour_plotter.m.prev > PDF_scatter_contour_plotter.m

   # Remove PDF_scatter_contour_plotter.m.prev.
   rm PDF_scatter_contour_plotter.m.prev

   # Run the plotter for the LBA case.
   ./run_scatter_contour_plots.bash \
      ../../sam_benchmark_runs/LBA_r1663_128x128x128_1km_Morrison/3D_output/LBA_128kmx128kmx128_1km_Morrison_64_0000003600_micro.nc \
      ../../output/$date_string/DDL/lba_zt.nc \
      ../../output/$date_string/DL/lba_zt.nc \
      ../../output/$date_string/SL/lba_zt.nc

   # Restore original PDF_scatter_contour_plotter.m.
   mv ../PDF_scatter_contour_plotter.m PDF_scatter_contour_plotter.m

   # Change directories back to the one where the run script is located.
   cd ../../run_scripts

fi # Produce plots of the PDFs.

# Run the goodness-of-fit tests.
if [ $FIT_TESTS == true ];
then

   # Change directory to the matlab_goodness_of_fit_tests directory.
   cd ../postprocessing/matlab_goodness_of_fit_tests

   # Copy goodness_of_fit_tests.m back one directory.
   cp goodness_of_fit_tests.m ../goodness_of_fit_tests.m

   # Setup the goodness-of-fit tests for the RICO case (rr and Nr).

   # Copy goodness_of_fit_tests.m to goodness_of_fit_tests.m.prev.
   cp goodness_of_fit_tests.m goodness_of_fit_tests.m.prev

   # Edit goodness_of_fit_tests.m.
   sed "1,100 s/clubb_start_height_idx =/%clubb_start_height_idx =/g; \
        1,100 s/clubb_stop_height_idx =/%clubb_stop_height_idx =/g; \
        1,100 s/sam_start_file_idx =/%sam_start_file_idx =/g; \
        1,100 s/sam_stop_file_idx =/%sam_stop_file_idx =/g; s/%%/%/g; \
        s/%clubb_start_height_idx = 2;/clubb_start_height_idx = 2;/g; \
        s/%clubb_stop_height_idx = 71;/clubb_stop_height_idx = 71;/g; \
        s/%sam_start_file_idx = 1;/sam_start_file_idx = 1;/g; \
        s/%sam_stop_file_idx = 3;/sam_stop_file_idx = 3;/g; \
        1,100 s/test_w     = true/test_w     = false/g; \
        1,100 s/test_rt    = true/test_rt    = false/g; \
        1,100 s/test_thl   = true/test_thl   = false/g; \
        1,100 s/test_chi   = true/test_chi   = false/g; \
        1,100 s/test_eta   = true/test_eta   = false/g;" \
        goodness_of_fit_tests.m.prev > goodness_of_fit_tests.m

   # Remove goodness_of_fit_tests.m.prev.
   rm goodness_of_fit_tests.m.prev

   # Run the goodness-of-fit tests for the RICO case.
   ./run_goodness_of_fit_tests.bash \
      ../../sam_benchmark_runs/RICO_256x256x100_drizzle/3D_output/ \
      ../../output/$date_string/DDL/rico_zt.nc \
      ../../output/$date_string/DL/rico_zt.nc \
      ../../output/$date_string/SL/rico_zt.nc \
      > ../../output/$date_string/goodness_of_fit_tests_rico_rr_Nr.txt

   # Setup the goodness-of-fit tests for the DYCOMS-II RF02 DS case (rr and Nr).

   # Copy goodness_of_fit_tests.m to goodness_of_fit_tests.m.prev.
   cp goodness_of_fit_tests.m goodness_of_fit_tests.m.prev

   # Edit goodness_of_fit_tests.m.
   sed "1,100 s/clubb_start_height_idx =/%clubb_start_height_idx =/g; \
        1,100 s/clubb_stop_height_idx =/%clubb_stop_height_idx =/g; \
        1,100 s/sam_start_file_idx =/%sam_start_file_idx =/g; \
        1,100 s/sam_stop_file_idx =/%sam_stop_file_idx =/g; \
        1,100 s/casename =/%casename =/g; s/%%/%/g; \
        s/%clubb_start_height_idx = 23;/clubb_start_height_idx = 23;/g; \
        s/%clubb_stop_height_idx = 56;/clubb_stop_height_idx = 56;/g; \
        s/%sam_start_file_idx = 1;/sam_start_file_idx = 1;/g; \
        s/%sam_stop_file_idx = 3;/sam_stop_file_idx = 3;/g; \
        s/%casename = 'RF02'/casename = 'RF02'/g;" \
        goodness_of_fit_tests.m.prev > goodness_of_fit_tests.m

   # Remove goodness_of_fit_tests.m.prev.
   rm goodness_of_fit_tests.m.prev

   # Run the goodness-of-fit tests for the DYCOMS-II RF02 DS case.
   ./run_goodness_of_fit_tests.bash \
      ../../sam_benchmark_runs/DYCOMS_RF02_128x128x96_dr_sed/3D_output/ \
      ../../output/$date_string/DDL/dycoms2_rf02_ds_zt.nc \
      ../../output/$date_string/DL/dycoms2_rf02_ds_zt.nc \
      ../../output/$date_string/SL/dycoms2_rf02_ds_zt.nc \
      > ../../output/$date_string/goodness_of_fit_tests_rf02_rr_Nr.txt

   # Setup the goodness-of-fit tests for the LBA case (rr and Nr).

   # Copy goodness_of_fit_tests.m to goodness_of_fit_tests.m.prev.
   cp goodness_of_fit_tests.m goodness_of_fit_tests.m.prev

   # Edit goodness_of_fit_tests.m.
   sed "1,100 s/clubb_start_height_idx =/%clubb_start_height_idx =/g; \
        1,100 s/clubb_stop_height_idx =/%clubb_stop_height_idx =/g; \
        1,100 s/sam_start_file_idx =/%sam_start_file_idx =/g; \
        1,100 s/sam_stop_file_idx =/%sam_stop_file_idx =/g; \
        1,100 s/casename =/%casename =/g; s/%%/%/g; \
        s/%clubb_start_height_idx = 2;/clubb_start_height_idx = 2;/g; \
        s/%clubb_stop_height_idx = 42;/clubb_stop_height_idx = 42;/g; \
        s/%sam_start_file_idx = 1;/sam_start_file_idx = 1;/g; \
        s/%sam_stop_file_idx = 6;/sam_stop_file_idx = 6;/g; \
        s/%casename = 'LBA'/casename = 'LBA'/g; \
        1,100 s/test_ri    = true/test_ri    = false/g; \
        1,100 s/test_Ni    = true/test_Ni    = false/g; \
        1,100 s/test_ri_ip = true/test_ri_ip = false/g; \
        1,100 s/test_Ni_ip = true/test_Ni_ip = false/g; \
        1,100 s/test_rs    = true/test_rs    = false/g; \
        1,100 s/test_Ns    = true/test_Ns    = false/g; \
        1,100 s/test_rs_ip = true/test_rs_ip = false/g; \
        1,100 s/test_Ns_ip = true/test_Ns_ip = false/g; \
        1,100 s/test_rg    = true/test_rg    = false/g; \
        1,100 s/test_Ng    = true/test_Ng    = false/g; \
        1,100 s/test_rg_ip = true/test_rg_ip = false/g; \
        1,100 s/test_Ng_ip = true/test_Ng_ip = false/g;" \
        goodness_of_fit_tests.m.prev > goodness_of_fit_tests.m

   # Remove goodness_of_fit_tests.m.prev.
   rm goodness_of_fit_tests.m.prev

   # Run the goodness-of-fit tests for the LBA case.
   ./run_goodness_of_fit_tests.bash \
      ../../sam_benchmark_runs/LBA_r1663_128x128x128_1km_Morrison/3D_output/ \
      ../../output/$date_string/DDL/lba_zt.nc \
      ../../output/$date_string/DL/lba_zt.nc \
      ../../output/$date_string/SL/lba_zt.nc \
      > ../../output/$date_string/goodness_of_fit_tests_lba_rr_Nr.txt

   # Setup the goodness-of-fit tests for the LBA case (rg and Ng).

   # Copy goodness_of_fit_tests.m to goodness_of_fit_tests.m.prev.
   cp goodness_of_fit_tests.m goodness_of_fit_tests.m.prev

   # Edit goodness_of_fit_tests.m.
   sed "1,100 s/clubb_start_height_idx =/%clubb_start_height_idx =/g; \
        1,100 s/clubb_stop_height_idx =/%clubb_stop_height_idx =/g; \
        1,100 s/sam_start_file_idx =/%sam_start_file_idx =/g; \
        1,100 s/sam_stop_file_idx =/%sam_stop_file_idx =/g; s/%%/%/g; \
        s/%clubb_start_height_idx = 34;/clubb_start_height_idx = 34;/g; \
        s/%clubb_stop_height_idx = 57;/clubb_stop_height_idx = 57;/g; \
        s/%sam_start_file_idx = 3;/sam_start_file_idx = 3;/g; \
        s/%sam_stop_file_idx = 6;/sam_stop_file_idx = 6;/g; \
        s/test_rr    = true/test_rr    = false/g; \
        s/test_Nr    = true/test_Nr    = false/g; \
        s/test_rr_ip = true/test_rr_ip = false/g; \
        s/test_Nr_ip = true/test_Nr_ip = false/g;
        1,100 s/test_ri    = true/test_ri    = false/g; \
        1,100 s/test_Ni    = true/test_Ni    = false/g; \
        1,100 s/test_ri_ip = true/test_ri_ip = false/g; \
        1,100 s/test_Ni_ip = true/test_Ni_ip = false/g; \
        1,100 s/test_rs    = true/test_rs    = false/g; \
        1,100 s/test_Ns    = true/test_Ns    = false/g; \
        1,100 s/test_rs_ip = true/test_rs_ip = false/g; \
        1,100 s/test_Ns_ip = true/test_Ns_ip = false/g; \
        1,100 s/test_rg    = false/test_rg    = true/g; \
        1,100 s/test_Ng    = false/test_Ng    = true/g; \
        1,100 s/test_rg_ip = false/test_rg_ip = true/g; \
        1,100 s/test_Ng_ip = false/test_Ng_ip = true/g;" \
        goodness_of_fit_tests.m.prev > goodness_of_fit_tests.m

   # Remove goodness_of_fit_tests.m.prev.
   rm goodness_of_fit_tests.m.prev

   # Run the goodness-of-fit tests for the LBA case.
   ./run_goodness_of_fit_tests.bash \
      ../../sam_benchmark_runs/LBA_r1663_128x128x128_1km_Morrison/3D_output/ \
      ../../output/$date_string/DDL/lba_zt.nc \
      ../../output/$date_string/DL/lba_zt.nc \
      ../../output/$date_string/SL/lba_zt.nc \
      > ../../output/$date_string/goodness_of_fit_tests_lba_rg_Ng.txt

   # Setup the goodness-of-fit tests for the LBA case (rs and Ns).

   # Copy goodness_of_fit_tests.m to goodness_of_fit_tests.m.prev.
   cp goodness_of_fit_tests.m goodness_of_fit_tests.m.prev

   # Edit goodness_of_fit_tests.m.
   sed "1,100 s/clubb_start_height_idx =/%clubb_start_height_idx =/g; \
        1,100 s/clubb_stop_height_idx =/%clubb_stop_height_idx =/g; \
        1,100 s/sam_start_file_idx =/%sam_start_file_idx =/g; \
        1,100 s/sam_stop_file_idx =/%sam_stop_file_idx =/g; s/%%/%/g; \
        s/%clubb_start_height_idx = 38;/clubb_start_height_idx = 38;/g; \
        s/%clubb_stop_height_idx = 54;/clubb_stop_height_idx = 54;/g; \
        s/%sam_start_file_idx = 5;/sam_start_file_idx = 5;/g; \
        s/%sam_stop_file_idx = 6;/sam_stop_file_idx = 6;/g; \
        s/test_rr    = true/test_rr    = false/g; \
        s/test_Nr    = true/test_Nr    = false/g; \
        s/test_rr_ip = true/test_rr_ip = false/g; \
        s/test_Nr_ip = true/test_Nr_ip = false/g; \
        1,100 s/test_ri    = true/test_ri    = false/g; \
        1,100 s/test_Ni    = true/test_Ni    = false/g; \
        1,100 s/test_ri_ip = true/test_ri_ip = false/g; \
        1,100 s/test_Ni_ip = true/test_Ni_ip = false/g; \
        1,100 s/test_rs    = false/test_rs    = true/g; \
        1,100 s/test_Ns    = false/test_Ns    = true/g; \
        1,100 s/test_rs_ip = false/test_rs_ip = true/g; \
        1,100 s/test_Ns_ip = false/test_Ns_ip = true/g; \
        1,100 s/test_rg    = true/test_rg    = false/g; \
        1,100 s/test_Ng    = true/test_Ng    = false/g; \
        1,100 s/test_rg_ip = true/test_rg_ip = false/g; \
        1,100 s/test_Ng_ip = true/test_Ng_ip = false/g;" \
        goodness_of_fit_tests.m.prev > goodness_of_fit_tests.m

   # Remove goodness_of_fit_tests.m.prev.
   rm goodness_of_fit_tests.m.prev

   # Run the goodness-of-fit tests for the LBA case.
   ./run_goodness_of_fit_tests.bash \
      ../../sam_benchmark_runs/LBA_r1663_128x128x128_1km_Morrison/3D_output/ \
      ../../output/$date_string/DDL/lba_zt.nc \
      ../../output/$date_string/DL/lba_zt.nc \
      ../../output/$date_string/SL/lba_zt.nc \
      > ../../output/$date_string/goodness_of_fit_tests_lba_rs_Ns.txt

   # Setup the goodness-of-fit tests for the LBA case (ri and Ni).

   # Copy goodness_of_fit_tests.m to goodness_of_fit_tests.m.prev.
   cp goodness_of_fit_tests.m goodness_of_fit_tests.m.prev

   # Edit goodness_of_fit_tests.m.
   sed "1,100 s/clubb_start_height_idx =/%clubb_start_height_idx =/g; \
        1,100 s/clubb_stop_height_idx =/%clubb_stop_height_idx =/g; \
        1,100 s/sam_start_file_idx =/%sam_start_file_idx =/g; \
        1,100 s/sam_stop_file_idx =/%sam_stop_file_idx =/g; s/%%/%/g; \
        s/%clubb_start_height_idx = 59;/clubb_start_height_idx = 59;/g; \
        s/%clubb_stop_height_idx = 65;/clubb_stop_height_idx = 65;/g; \
        s/%sam_start_file_idx = 6;/sam_start_file_idx = 6;/g; \
        s/%sam_stop_file_idx = 6;/sam_stop_file_idx = 6;/g; \
        s/test_rr    = true/test_rr    = false/g; \
        s/test_Nr    = true/test_Nr    = false/g; \
        s/test_rr_ip = true/test_rr_ip = false/g; \
        s/test_Nr_ip = true/test_Nr_ip = false/g; \
        1,100 s/test_ri    = false/test_ri    = true/g; \
        1,100 s/test_Ni    = false/test_Ni    = true/g; \
        1,100 s/test_ri_ip = false/test_ri_ip = true/g; \
        1,100 s/test_Ni_ip = false/test_Ni_ip = true/g; \
        1,100 s/test_rs    = true/test_rs    = false/g; \
        1,100 s/test_Ns    = true/test_Ns    = false/g; \
        1,100 s/test_rs_ip = true/test_rs_ip = false/g; \
        1,100 s/test_Ns_ip = true/test_Ns_ip = false/g; \
        1,100 s/test_rg    = true/test_rg    = false/g; \
        1,100 s/test_Ng    = true/test_Ng    = false/g; \
        1,100 s/test_rg_ip = true/test_rg_ip = false/g; \
        1,100 s/test_Ng_ip = true/test_Ng_ip = false/g;" \
        goodness_of_fit_tests.m.prev > goodness_of_fit_tests.m

   # Remove goodness_of_fit_tests.m.prev.
   rm goodness_of_fit_tests.m.prev

   # Run the goodness-of-fit tests for the LBA case.
   ./run_goodness_of_fit_tests.bash \
      ../../sam_benchmark_runs/LBA_r1663_128x128x128_1km_Morrison/3D_output/ \
      ../../output/$date_string/DDL/lba_zt.nc \
      ../../output/$date_string/DL/lba_zt.nc \
      ../../output/$date_string/SL/lba_zt.nc \
      > ../../output/$date_string/goodness_of_fit_tests_lba_ri_Ni.txt

   # Restore original goodness_of_fit_tests.m.
   mv ../goodness_of_fit_tests.m goodness_of_fit_tests.m

   # Change directories back to the one where the run script is located.
   cd ../../run_scripts

fi # Run the goodness-of-fit tests.

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
   ./plotgen.pl -c -l -e -t ../../output/$date_string/DDL \
                ../../output/$date_string/DL \
                ../../output/$date_string/SL \
                ../../output/$date_string/plot_DDL_DL_SL_mc_profiles

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
