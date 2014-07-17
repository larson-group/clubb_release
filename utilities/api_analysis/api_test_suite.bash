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

echo "removing .svn folder from CLUBB_core"
rm -rf CLUBB_core/.svn

echo "removing api from CLUBB_core"
rm  CLUBB_core/clubb_api_module.F90

echo "removing api from SILHS"
rm  CLUBB/src/SILHS/silhs_api_module.F90

echo "removing G_Unit_Tests from CLUBB"
rm -rf CLUBB/src/G_unit_test_types
rm CLUBB/src/G_unit_tests.F90

echo "moving SILHS to CLUBB_core"
mv CLUBB/src/SILHS/* CLUBB_core/

echo "testing CLUBB_standalone"
python api_commitment_test.py -cpu CLUBB_core CLUBB > clubb_standalone_modules.txt

echo "testing CLUBB_core"
python api_commitment_test.py -cpu CLUBB_core CLUBB_core > clubb_core_modules.txt

echo "testing SAM"
python api_commitment_test.py -cpu CLUBB_core SAM/SRC --exclude-dir="CLUBB","SILHS" > sam_modules.txt

echo "testing CAM"
python api_commitment_test.py -cpu CLUBB_core CAM --exclude-dir="clubb","silhs" > cam_modules.txt

echo "testing WRF"
python api_commitment_test.py -cpu CLUBB_core WRF --exclude-dir="clubb","silhs" > wrf_modules.txt

python create_module_table.py CLUBB_core> output.html

echo "Removing Dependencies"
rm -rf CLUBB
rm -rf SAM
rm -rf WRF
rm -rf CAM
rm -rf CLUBB_core

rm -rf clubb_standalone_modules.txt
rm -rf clubb_core_modules.txt
rm -rf sam_modules.txt
rm -rf wrf_modules.txt
rm -rf cam_modules.txt
