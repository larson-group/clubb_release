#!/bin/csh

# Note:  The NCO package must be loaded before running this script using:
#        module load nco

# Specify the start year and end year of the E3SM run that you want the
# climatology produced for.
setenv START_YEAR 1
setenv END_YEAR 7

# Specify the map used for regridded climatology.
setenv REGRID_MAP /home/ac.zender/data/maps/map_ne30pg2_to_cmip6_180x360_aave.20200201.nc

# Specify the case name.
setenv CASE_NAME default.alpha22.F2010.clubb_c_K10h_0p35.ne30pg2_r05_oECv3.anvil

# The INPUT_DIR is the directory where the input to ncclimo is stored.
# This is the directory where the raw output from the E3SM model run is stored,
# since the raw model output is used as the input to ncclimo.
setenv INPUT_DIR /lcrc/group/acme/ac.griffin/E3SM_simulations/${CASE_NAME}/run/

# The OUTPUT_DIR is the directory where the climatology files that are output
# from ncclimo are stored.
setenv OUTPUT_DIR ${INPUT_DIR}climatology_yrs_${START_YEAR}_${END_YEAR}/

# The REGRID_OUTPUT_DIR is the directory where the regridded climatology files
# that are output from ncclimo using the map specificied in REGRID_MAP are
# stored.
setenv REGRID_OUTPUT_DIR ${OUTPUT_DIR}regridded_climo/

# Call ncclimo to generate climatology and regridded climatology files.
# The "-a sdd" option is for a seasonally discontinuous December.
ncclimo -a sdd -s ${START_YEAR} -e ${END_YEAR} -r ${REGRID_MAP} -c ${CASE_NAME} -i ${INPUT_DIR} -o ${OUTPUT_DIR} -O ${REGRID_OUTPUT_DIR}
