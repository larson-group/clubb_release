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
analyzeUsageTxt="log/usageAnalyzerTable.txt"
clubbStandaloneTxt="clubb_standalone_modules.txt"
clubbCoreTxt="clubb_core_modules.txt"
samTxt="sam_modules.txt"
camTxt="cam_modules.txt"
wrfTxt="wrf_modules.txt"
tableOut="log/apiCommitmentTable.txt"
clubbStandalone="$clubbDir/src/clubb_standalone.F90"

echo "###Checking out CLUBB###"
git clone git@github.com:larson-group/clubb.git $clubbDir
echo ""
echo "###Checking out SAM###"
git clone git@github.com:larson-group/sam_clubb.git $samDir
echo ""
echo "###Checking out WRF###"
git clone git@github.com:larson-group/wrf.git $wrfDir
echo ""
echo "###Checking out CAM###"
git clone git@github.com:larson-group/cam.git $camDir
echo ""
#fi

echo "###Moving CLUBB_core###"
mv $clubbDir/src/CLUBB_core CLUBB_core
echo ""

echo "###Removing .git Folders###"
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
echo ""

# All four calls to api_commitment_test.py and the one to usage_analyzer.py create output files
# that serve as input for the create_module_table.py script
# In addition, the exit codes of the 3 host model calls to api_commitment indicate success or failure,
# where success means only clubb_api_module is imported,
# while failure means another imported CLUBB module was found.
# The script either return 0 (success) or 1 (failure)
echo "###Copying the SILHS API into the CLUBB API###"
cat $clubbDir/src/SILHS/silhs_api_module.F90 CLUBB_core/clubb_api_module.F90 > tempApiFile
mv tempApiFile CLUBB_core/clubb_api_module.F90
echo ""

#echo "###Running the Usage Analyzer###"
## This script will create a table that shows which API element of CLUBB is used in which module in
## each of the 3 host models SAM, WRF, and CAM
## For Jenkins we don't need this. Skip!
# Create log folder if necessary
#if [ ! -d "log" ]; then mkdir log; fi
#python usage_analyzer.py CLUBB_core/clubb_api_module.F90 $samDir $wrfDir/WRF $camDir > $analyzeUsageTxt
#echo ""

echo "###Removing API from CLUBB_core###"
rm  CLUBB_core/clubb_api_module.F90
echo ""

echo "###Removing API from SILHS###"
rm  $clubbDir/src/SILHS/silhs_api_module.F90
echo ""

echo "###Removing G_Unit_Tests from CLUBB###"
rm -rf $clubbDir/src/G_unit_test_types
rm $clubbDir/src/G_unit_tests.F90
echo ""

echo "###Moving SILHS to CLUBB_core###"
mv $clubbDir/src/SILHS/* CLUBB_core/
echo ""

#echo "###Testing CLUBB_Standalone's API Commitment###"
# We figured out that clubb_standalone is not a good emulation of a host model
# because it does not directly call advance_clubb_core and does too much extra stuff.
# So this test does not work :(
#python api_commitment_test.py -cpu CLUBB_core $clubbStandalone > $clubbStandaloneTxt
#standaloneExit=$?
#cat $clubbStandaloneTxt
#echo ""

#echo "###Testing CLUBB_core's API Commitment###"
#python api_commitment_test.py -cpu CLUBB_core CLUBB_core > $clubbCoreTxt
#cat $clubbCoreTxt
#echo ""

echo "###Testing SAM's API Commitment###"
python api_commitment_test.py -cpu CLUBB_core $samDir --exclude-dir CLUBB SILHS > $samTxt
# $? is a builtin bash variable storing the exit code of the last command
# Here, this is the same value that is passed to sys.exit() within the api_commitment_test.py script
samExit=$?
# Print output to console (and therefore the Jenkins site)
cat $samTxt
echo ""

echo "###Testing CAM's API Commitment###"
python api_commitment_test.py -cpu CLUBB_core $camDir --exclude-dir spcam cime clubb silhs > $camTxt
# $? is a builtin bash variable storing the exit code of the last command
# Here, this is the same value that is passed to sys.exit() within the api_commitment_test.py script
camExit=$?
# Print output to console (and therefore the Jenkins site)
cat $camTxt
echo ""

echo "###Testing WRF's API Commitment###"
python api_commitment_test.py -cpu CLUBB_core $wrfDir/WRF --exclude-dir clubb silhs > $wrfTxt
# $? is a builtin bash variable storing the exit code of the last command
# Here, this is the same value that is passed to sys.exit() within the api_commitment_test.py script
wrfExit=$?
# Print output to console (and therefore the Jenkins site)
cat $wrfTxt
echo ""

#echo "###Creating Table Output###"
## For Jenkins we don't need this. Skip!
#python create_module_table.py CLUBB_core $clubbStandaloneTxt $clubbCoreTxt $samTxt $camTxt $wrfTxt > $tableOut
#echo ""

echo "###Removing Checkouts###"
rm -rf $clubbDir
rm -rf $samDir
rm -rf $wrfDir
rm -rf $camDir
rm -rf CLUBB_core
echo ""

echo "###Testing API Commitment###"
# Check exit states for all 3 api_commitment calls
# Since 0 means success, any sum value greater than 0 means one of the three scripts found an issue
# The details can be found in the output of the scripts above.
echo ""
if [[ $standaloneExit+$samExit+$camExit+$wrfExit > 0 ]] ; then
   echo "#########################################################"
   echo "######ERROR: A host model is circumventing the API.######"
   echo "#########################################################"
#   echo "Please inspect the nightly test API Commitment Table."
#   echo "MODELS AT FAULT INCLUDE:"

#   if [[ -s "$sam_modules" ]] ; then
#        echo " - SAM Model"
#        cat $sam_modules
#   fi
#   if [[ -s "$wrf_modules" ]] ; then
#        echo " - WRF Model"
#        cat $wrf_modules
#   fi
#   if [[ -s "$cam_modules" ]] ; then
#        echo " - CAM Model"
#        cat $cam_modules
#   fi
   exitCode=1
else
   echo "###################################"
   echo "######All host models passed.######"
   echo "###################################"
   exitCode=0
fi
echo ""
echo ""

echo "###Removing Dependencies###"
rm -rf $clubbStandaloneTxt
rm -rf $clubbCoreTxt
rm -rf $samTxt
rm -rf $camTxt
rm -rf $wrfTxt

cd $restoreLoc

exit $exitCode

