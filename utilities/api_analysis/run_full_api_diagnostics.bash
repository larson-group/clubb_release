#!/bin/bash

####################################################################################
# run_full_api_diagnostics v1.3
#
# Run every api/* test.
#
# File History:
#  v1.3: Updated nightly tests variables
#  v1.2: Converted svn commands to git
#  v1.1: Merged the bitten and nightly version of this script
#  v1.0: Initial Release
#
# Required Arguments
#  -None
# Optional Arguments
#  $1 = "-nightly"
#         With this argument, the tables are produced using the nightly paths.
#         Without this argument, no tables are produced (bitten requires no tables).
####################################################################################

#This variable allows the script to return to where it was run from.
restoreLoc=`pwd`

# Figure out the directory where the script is located.
run_dir=`dirname $0`

# Change directories to the one in which the script is located.
cd $run_dir

## Checkout CLUBB and the Host Models
#if [[ "$1" == "-nightly" ]]
#then
#    # Setup the paths to the host models
#    clubbDir=$clubbSource
#    samDir=$samSource
#    wrfDir=$wrfSource
#    camDir=$camSource
#
#    # Setup the checkout scripts
#    echo "Checking out CLUBB"
#    $buildDir/checkout_fresh.bash "CLUBB"
#    echo "Checking out SAM"
#    $buildDir/checkout_fresh.bash "SAM"
#    echo "Checking out WRF"
#    $buildDir/checkout_fresh.bash "WRF"
#    echo "Checking out CAM"
#    $buildDir/checkout_fresh.bash "CAM"
#else
    # Setup the paths to the host models
    clubbDir="CLUBB"
    samDir="SAM"
    wrfDir="WRF"
    camDir="CAM"

    echo "Checking out CLUBB"
    git clone https://github.com/larson-group/clubb.git $clubbDir
    echo "Checking out SAM"
    git clone https://github.com/larson-group/sam_clubb.git $samDir
    echo "Checking out WRF"
    git clone https://github.com/larson-group/wrf.git $wrfDir
    echo "Checking out CAM"
    git clone https://github.com/larson-group/cam.git $camDir
#fi

echo "Moving CLUBB_core"
mv $clubbDir/src/CLUBB_core CLUBB_core

echo "Removing .git Folders"
find CLUBB_core -type d -name .git -exec rm -rf {} \;
find $clubbDir -type d -name .git -exec rm -rf {} \;
find $samDir -type d -name .git -exec rm -rf {} \;
find $wrfDir -type d -name .git -exec rm -rf {} \;
find $camDir -type d -name .git -exec rm -rf {} \;

find CLUBB_core -type d -name .gitignore -exec rm -rf {} \;
find $clubbDir -type d -name .gitignore -exec rm -rf {} \;
find $samDir -type d -name .gitignore -exec rm -rf {} \;
find $wrfDir -type d -name .gitignore -exec rm -rf {} \;
find $camDir -type d -name .gitignore -exec rm -rf {} \;

echo "Copying the SILHS API into the CLUBB API"
cat $clubbDir/src/SILHS/silhs_api_module.F90 CLUBB_core/clubb_api_module.F90 > tempApiFile
mv tempApiFile CLUBB_core/clubb_api_module.F90

echo "Running the Usage Analyzer"
python usage_analyzer.py CLUBB_core/clubb_api_module.F90 $samDir $wrfDir/WRF $camDir > log/usageAnalyzerTable.txt

echo "Removing API from CLUBB_core"
rm  CLUBB_core/clubb_api_module.F90

echo "Removing API from SILHS"
rm  $clubbDir/src/SILHS/silhs_api_module.F90

echo "Removing G_Unit_Tests from CLUBB"
rm -rf $clubbDir/src/G_unit_test_types
rm $clubbDir/src/G_unit_tests.F90

echo "Moving SILHS to CLUBB_core"
mv $clubbDir/src/SILHS/* CLUBB_core/

echo "Testing CLUBB_Standalone's API Commitment"
python api_commitment_test.py -cpu CLUBB_core $clubbDir > clubb_standalone_modules.txt

echo "Testing CLUBB_core's API Commitment"
python api_commitment_test.py -cpu CLUBB_core CLUBB_core > clubb_core_modules.txt

echo "Testing SAM's API Commitment"
python api_commitment_test.py -cpu CLUBB_core $samDir --exclude-dir CLUBB SILHS > sam_modules.txt

echo "Testing CAM's API Commitment"
python api_commitment_test.py -cpu CLUBB_core $camDir --exclude-dir spcam cime clubb silhs > cam_modules.txt

echo "Testing WRF's API Commitment"
python api_commitment_test.py -cpu CLUBB_core $wrfDir/WRF --exclude-dir clubb silhs > wrf_modules.txt

python create_module_table.py CLUBB_core > log/apiCommitmentTable.txt
echo "Removing Checkouts"
rm -rf $clubbDir
rm -rf $samDir
rm -rf $wrfDir
rm -rf $camDir
rm -rf CLUBB_core
    
echo "Testing API Commitment"
sam_modules="sam_modules.txt"
wrf_modules="wrf_modules.txt"
cam_modules="cam_modules.txt"
if [[ -s "$sam_modules" ]] || [[ -s "$wrf_modules" ]] || [[ -s "$cam_modules" ]] ; then
   echo "ERROR: A host model is circumventing the API."
   echo "Please inspect the nightly test API Commitment Table."
   echo "MODELS AT FAULT INCLUDE:"

   if [[ -s "$sam_modules" ]] ; then
        echo " - SAM Model"
        cat $sam_modules
   fi
   if [[ -s "$wrf_modules" ]] ; then
        echo " - WRF Model"
        cat $wrf_modules
   fi
   if [[ -s "$cam_modules" ]] ; then
        echo " - CAM Model"
        cat $cam_modules
   fi
   exitCode=1
else 
   echo "All host models passed."
   exitCode=0
fi


echo "Removing Dependencies"
rm -rf clubb_standalone_modules.txt
rm -rf clubb_core_modules.txt
rm -rf sam_modules.txt
rm -rf wrf_modules.txt
rm -rf cam_modules.txt

cd $restoreLoc

exit $exitCode

