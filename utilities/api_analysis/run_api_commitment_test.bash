# Runs the API Commitment Tests and Returns 0 if all host models only use the CLUBB API
# Author: Orion Schemenauer

# Figure out the directory where the script is located.
run_dir=`dirname $0`

# Change directories to the one in which the script is located.
cd $run_dir

echo "Checking out CLUBB"
svn co http://carson.math.uwm.edu/repos/clubb_repos/trunk CLUBB

echo "Checking out SAM"
svn co http://carson.math.uwm.edu/repos/sam_repos/trunk SAM

echo "Checking out WRF"
svn co http://carson.math.uwm.edu/repos/wrf/trunk WRF

echo "Checking out CAM"
svn co https://svn-ccsm-models.cgd.ucar.edu/cam1/branches/subcol_SILHS_UWM CAM

echo "Moving CLUBB_core"
mv CLUBB/src/CLUBB_core CLUBB_core

echo "Removing .svn Folders"
find . -type d -name .svn -exec rm -rf {} \;

echo "Removing API from CLUBB_core"
rm  CLUBB_core/clubb_api_module.F90

echo "Removing API from SILHS"
rm  CLUBB/src/SILHS/silhs_api_module.F90

echo "Removing G_Unit_Tests from CLUBB"
rm -rf CLUBB/src/G_unit_test_types
rm CLUBB/src/G_unit_tests.F90

echo "Moving SILHS to CLUBB_core"
mv CLUBB/src/SILHS/* CLUBB_core/

echo "Testing CLUBB_Standalone's API Commitment"
python api_commitment_test.py -cpu CLUBB_core CLUBB > clubb_standalone_modules.txt

echo "Testing CLUBB_core's API Commitment"
python api_commitment_test.py -cpu CLUBB_core CLUBB_core > clubb_core_modules.txt

echo "Testing SAM's API Commitment"
python api_commitment_test.py -cpu CLUBB_core SAM/SRC --exclude-dir="CLUBB","SILHS" > sam_modules.txt

echo "Testing CAM's API Commitment"
python api_commitment_test.py -cpu CLUBB_core CAM --exclude-dir="clubb","silhs" > cam_modules.txt

echo "Testing WRF's API Commitment"
python api_commitment_test.py -cpu CLUBB_core WRF --exclude-dir="clubb","silhs" > wrf_modules.txt

python create_module_table.py CLUBB_core> apiCommitment_.html

echo "Removing Checkouts"
rm -rf CLUBB
rm -rf SAM
rm -rf WRF
rm -rf CAM
rm -rf CLUBB_core

# Are any host model files > 0 size? If so, they failed
sam_modules="sam_modules.txt"
wrf_modules="wrf_modules.txt"
cam_modules="cam_modules.txt"
if [ -s "$sam_modules" ] || [ -s "$wrf_modules" ] || [ -s "$cam_modules" ] ; then
    echo "ERROR: A host model is circumavigating the API." 
    echo "Please inspect the nightly api_commitment_chart."
    exit 1
fi

#If it made it this far - all is well
echo "All host models passed."
exit 0
