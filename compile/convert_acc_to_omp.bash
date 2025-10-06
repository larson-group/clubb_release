#!/bin/bash

# This script takes a branch of CLUBB with openacc directives and adds openmp directives.
# To run this script, go to CLUBB's compile directory of the branch you wish to convert and type
# ./convert_acc_to_omp.bash

# For further information about this script, see
#    https://github.com/larson-group/clubb/issues/1138

# Move to script directory
cd "$(dirname "$0")"

if [ ! -d $PWD/intel-application-migration-tool-for-openacc-to-openmp/ ]; then
	git clone https://github.com/huebleruwm/intel-application-migration-tool-for-openacc-to-openmp
fi
ACC_TO_OMP_SCRIPT=$PWD/intel-application-migration-tool-for-openacc-to-openmp/src/intel-application-migration-tool-for-openacc-to-openmp

cd ../src

# Remove all instances of "default(present)" from !$acc directives.
find . -type f -exec sed -i "s: default(present)::g" {} \;

# Convert the openacc code to openmp directives using the converter script.
# Using -async=ignore prevents the script from converting openacc async/wait
# statements into openmp nowait/taskwait statements. This only affects SILHS
# but it causes code crashes in SILHS cases without this option.
$ACC_TO_OMP_SCRIPT -overwrite-input -async=ignore -no-declare-mapper *.F90
$ACC_TO_OMP_SCRIPT -overwrite-input -async=ignore -no-declare-mapper CLUBB_core/*.F90
$ACC_TO_OMP_SCRIPT -overwrite-input -async=ignore -no-declare-mapper SILHS/*.F90
$ACC_TO_OMP_SCRIPT -overwrite-input -async=ignore -no-declare-mapper Benchmark_cases/*.F90
$ACC_TO_OMP_SCRIPT -overwrite-input -async=ignore -no-declare-mapper Microphys/*.F90


# Remove the extra generated files
find . -type f \( -name "*.original" -o -name "*.report" -o -name "*.translated" \) -delete

# Replace the scripts usage of "OPENACC2OPENMP_ORIGINAL_OPENMP" for 
# ifdefs around existing openmp code, with "OPENMP_CPU".
find . -type f -exec sed -i "s:OPENACC2OPENMP_ORIGINAL_OPENMP:OPENMP_CPU:g"  {} \;

# The script mistranslates "acc host_data use_device(rand_pool)" to 
# "omp target update from(rand_pool)" but we want this to be "omp 
# target data use_device_ptr(rand_pool)"
sed -i "s:omp target update from(rand_pool):omp target data use_device_ptr(rand_pool):g" SILHS/latin_hypercube_driver_module.F90

# To match "omp target data use_device_ptr(rand_pool)", we need a 
# corresponding "omp end target data" to follow. We want only to 
# change the first instance of "omp target data use_device_ptr(rand_pool)" 
# though, so we have this nonsense bit ":a;$!{N;ba};" that somehow does that 
sed -i ':a;$!{N;ba}; s:omp target update to(rand_pool):omp end target data:' SILHS/latin_hypercube_driver_module.F90

# Return the "default(present)" clauses to all openacc statements
find . -type f -exec sed -i "s:acc parallel loop gang vector collapse(4):acc parallel loop gang vector collapse(4) default(present):g" {} \;
find . -type f -exec sed -i "s:acc parallel loop gang vector collapse(3):acc parallel loop gang vector collapse(3) default(present):g" {} \;
find . -type f -exec sed -i "s:acc parallel loop gang vector collapse(2):acc parallel loop gang vector collapse(2) default(present):g" {} \;
find . -type f -exec sed -i "s:acc parallel loop gang vector\s*$:acc parallel loop gang vector default(present):g"  {} \;

# Remove the translated code messages from the bottom of files
find . -type f -exec sed -i "s:^.*Code was translated using.*$::g"  {} \;
