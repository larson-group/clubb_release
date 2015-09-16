#!/bin/bash

####################################################################################
# run_full_api_diagnostics v1.0
#
# Run every api/* test.
#
# File History:
#  v1.1: Merged the bitten and nightly version of this script
#  v1.0: Initial Release
#
# Required Arguments
#  -None
# Optional Arguments
#  -nightly
#       With this argument, the tables are produced using the nightly paths.
#       Without this argument, no tables are produced (bitten requires no tables).
####################################################################################

#This variable allows the script to reurn to where it was run from.
restoreLoc=`pwd`

# Figure out the directory where the script is located.
run_dir=`dirname $0`

# Change directories to the one in which the script is located.
cd $run_dir

# Checkout CLUBB and the Host Models
if [ "$1" == "-nightly" ]
then
    # Setup the checkout scripts
    . ../nightly_config.sh
    echo "Checking out CLUBB"
    ../checkout_clubb_fresh.bash
    echo "Checking out SAM"
    ../checkout_sam_fresh.bash
    echo "Checking out WRF"
    ../checkout_wrf_fresh.bash
    echo "Checking out CAM"
    ../checkout_cam_fresh.bash
else
    echo "Checking out CLUBB"
    svn co http://carson.math.uwm.edu/repos/clubb_repos/trunk CLUBB
    echo "Checking out SAM"
    svn co http://carson.math.uwm.edu/repos/sam_repos/trunk SAM
    echo "Checking out WRF"
    svn co http://carson.math.uwm.edu/repos/wrf/trunk WRF
    echo "Checking out CAM"
    svn co --username vlarson@uwm.edu https://svn-ccsm-models.cgd.ucar.edu/cam1/branches/subcol_SILHS_UWM CAM
fi

# Setup the paths to the host models
if [ "$1" == "-nightly" ]
then
    clubbDir="../../../clubb"
    samDir="../../../SAM_CLUBB"
    wrfDir="../../../WRF_CLUBB"
    camDir="../../../CAM"
else
    clubbDir="CLUBB"
    samDir="SAM"
    wrfDir="WRF"
    camDir="CAM"
fi

echo "Moving CLUBB_core"
mv $clubbDir/src/CLUBB_core CLUBB_core

echo "Removing .svn Folders"
find CLUBB_core -type d -name .svn -exec rm -rf {} \;
find $clubbDir -type d -name .svn -exec rm -rf {} \;
find $samDir -type d -name .svn -exec rm -rf {} \;
find $wrfDir -type d -name .svn -exec rm -rf {} \;
find $camDir -type d -name .svn -exec rm -rf {} \;

if [ "$1" == "-nightly" ]
then
    echo "Copying the SILHS API into the CLUBB API"
    cat $clubbDir/src/SILHS/silhs_api_module.F90 CLUBB_core/clubb_api_module.F90 > tempApiFile
    mv tempApiFile CLUBB_core/clubb_api_module.F90

    echo "Running the Usage Analyzer"
    python usage_analyzer.py CLUBB_core/clubb_api_module.F90 $samDir $wrfDir $camDir > ../text_output/usageAnalyzerTable.html
fi

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
python api_commitment_test.py -cpu CLUBB_core $samDir --exclude-dir="CLUBB","SILHS" > sam_modules.txt

echo "Testing CAM's API Commitment"
python api_commitment_test.py -cpu CLUBB_core $camDir --exclude-dir="clubb","silhs" > cam_modules.txt

echo "Testing WRF's API Commitment"
python api_commitment_test.py -cpu CLUBB_core $wrfDir --exclude-dir="clubb","silhs" > wrf_modules.txt

if [ "$1" == "-nightly" ]
then
    python create_module_table.py CLUBB_core> ../text_output/apiCommitmentTable.html
    rm -rf CLUBB_core
else
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
    if [ -s "$sam_modules" ] || [ -s "$wrf_modules" ] || [ -s "$cam_modules" ] ; then
        echo "ERROR: A host model is circumventing the API." 
        echo "Please inspect the nightly test API Commitment Table."
        exitCode=1
    else 
        echo "All host models passed."
        exitCode=0
    fi
fi

echo "Removing Dependencies"
rm -rf clubb_standalone_modules.txt
rm -rf clubb_core_modules.txt
rm -rf sam_modules.txt
rm -rf wrf_modules.txt
rm -rf cam_modules.txt

cd $restoreLoc

exit $exitCode
