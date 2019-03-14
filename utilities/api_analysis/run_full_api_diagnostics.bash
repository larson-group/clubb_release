#!/bin/bash

####################################################################################
# run_full_api_diagnostics v1.2
#
# Run every api/* test.
#
# File History:
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
#  $1 = <log file>
#         $outputLog will be set to the log file specified (used by bitten tests).
####################################################################################

#This variable allows the script to return to where it was run from.
restoreLoc=`pwd`

# Figure out the directory where the script is located.
run_dir=`dirname $0`

# Change directories to the one in which the script is located.
cd $run_dir

# Checkout CLUBB and the Host Models
if [ "$1" == "-nightly" ]
then
    # The output log file
    outputLog=$logLocation/current/nightly_tests.log

    # Setup the checkout scripts
    echo "Checking out CLUBB" >> $outputLog
    ../checkout_clubb_fresh.bash
    echo "Checking out SAM" >> $outputLog
    ../checkout_sam_fresh.bash
    echo "Checking out WRF" >> $outputLog
    ../checkout_wrf_fresh.bash
    echo "Checking out CAM" >> $outputLog
    ../checkout_cam_fresh.bash
else
    # The output log file
    outputLog=$1

    echo "Checking out CLUBB" >> $outputLog
    git clone https://github.com/larson-group/clubb.git CLUBB
    echo "Checking out SAM" >> $outputLog
    git clone https://github.com/larson-group/sam_clubb.git SAM
    echo "Checking out WRF" >> $outputLog
    git clone https://github.com/larson-group/wrf.git WRF
    echo "Checking out CAM" >> $outputLog
    git clone https://github.com/larson-group/cam.git CAM
fi

# Setup the paths to the host models
if [ "$1" == "-nightly" ]
then
    clubbDir=$clubbSource
    samDir=$samSource
    wrfDir=$wrfSource
    camDir=$camSource
else
    clubbDir="CLUBB"
    samDir="SAM"
    wrfDir="WRF/WRF"
    camDir="CAM"
fi

echo "Moving CLUBB_core" >> $outputLog
mv $clubbDir/src/CLUBB_core CLUBB_core

echo "Removing .git Folders" >> $outputLog
find CLUBB_core -type d -name .git -exec rm -rf {} \; >> $outputLog
find $clubbDir -type d -name .git -exec rm -rf {} \; >> $outputLog
find $samDir -type d -name .git -exec rm -rf {} \; >> $outputLog
find $wrfDir -type d -name .git -exec rm -rf {} \; >> $outputLog
find $camDir -type d -name .git -exec rm -rf {} \; >> $outputLog

find CLUBB_core -type d -name .gitignore -exec rm -rf {} \; >> $outputLog
find $clubbDir -type d -name .gitignore -exec rm -rf {} \; >> $outputLog
find $samDir -type d -name .gitignore -exec rm -rf {} \; >> $outputLog
find $wrfDir -type d -name .gitignore -exec rm -rf {} \; >> $outputLog
find $camDir -type d -name .gitignore -exec rm -rf {} \; >> $outputLog

if [ "$1" == "-nightly" ]
then
    echo "Copying the SILHS API into the CLUBB API" >> $outputLog
    cat $clubbDir/src/SILHS/silhs_api_module.F90 CLUBB_core/clubb_api_module.F90 > tempApiFile
    mv tempApiFile CLUBB_core/clubb_api_module.F90

    echo "Running the Usage Analyzer" >> $outputLog
    python usage_analyzer.py CLUBB_core/clubb_api_module.F90 $samDir $wrfDir $camDir > ../text_output/usageAnalyzerTable.html
fi

echo "Removing API from CLUBB_core" >> $outputLog
rm  CLUBB_core/clubb_api_module.F90 >> $outputLog

echo "Removing API from SILHS" >> $outputLog
rm  $clubbDir/src/SILHS/silhs_api_module.F90

echo "Removing G_Unit_Tests from CLUBB" >> $outputLog
rm -rf $clubbDir/src/G_unit_test_types >> $outputLog
rm $clubbDir/src/G_unit_tests.F90 >> $outputLog

echo "Moving SILHS to CLUBB_core" >> $outputLog
mv $clubbDir/src/SILHS/* CLUBB_core/ >> $outputLog

echo "Testing CLUBB_Standalone's API Commitment" >> $outputLog
python api_commitment_test.py -cpu CLUBB_core $clubbDir > clubb_standalone_modules.txt

echo "Testing CLUBB_core's API Commitment" >> $outputLog
python api_commitment_test.py -cpu CLUBB_core CLUBB_core > clubb_core_modules.txt

echo "Testing SAM's API Commitment" >> $outputLog
python api_commitment_test.py -cpu CLUBB_core $samDir --exclude-dir="CLUBB","SILHS" > sam_modules.txt

echo "Testing CAM's API Commitment" >> $outputLog
python api_commitment_test.py -cpu CLUBB_core $camDir --exclude-dir="spcam","cime","clubb","silhs" > cam_modules.txt

echo "Testing WRF's API Commitment" >> $outputLog
python api_commitment_test.py -cpu CLUBB_core $wrfDir --exclude-dir="clubb","silhs" > wrf_modules.txt

if [ "$1" == "-nightly" ]
then
    python create_module_table.py CLUBB_core> ../text_output/apiCommitmentTable.html
    rm -rf CLUBB_core     
else
    echo "Removing Checkouts" >> $outputLog
    rm -rf $clubbDir
    rm -rf $samDir
    rm -rf $wrfDir
    rm -rf $camDir
    rm -rf CLUBB_core
    
    echo "Testing API Commitment" >> $outputLog
    sam_modules="sam_modules.txt"
    wrf_modules="wrf_modules.txt"
    cam_modules="cam_modules.txt"
    if [ -s "$sam_modules" ] || [ -s "$wrf_modules" ] || [ -s "$cam_modules" ] ; then
        echo "ERROR: A host model is circumventing the API." >> $outputLog
        echo "Please inspect the nightly test API Commitment Table." >> $outputLog
        echo "MODELS AT FAULT INCLUDE:"

                if [ -s "$sam_modules" ] ; then
                        echo " - SAM Model"
                        cat sam_modules.txt             
                fi
                if [ -s "$wrf_modules" ] ; then 
                        echo " - WRF Model"
                        cat wrf_modules.txt
                fi
                if [ -s "$cam_modules" ] ; then
                        echo " - CAM Model"
                        cat cam_modules.txt
                fi


        exitCode=1
    else 
        echo "All host models passed." >> $outputLog
        exitCode=0
    fi
fi

echo "Removing Dependencies" >> $outputLog
rm -rf clubb_standalone_modules.txt
rm -rf clubb_core_modules.txt
rm -rf sam_modules.txt
rm -rf wrf_modules.txt
rm -rf cam_modules.txt

cd $restoreLoc

exit $exitCode

