if [ ! -f /CLUBB/ ] || [ ! -f /SAM/ ] || [ ! -f /CAM/ ] || [ ! -f /WRF/ ] ; then
    echo "This script requires a freshly checked out version of CLUBB, SAM, CAM, and WRF, named respectively."
fi

echo "Moving CLUBB_core"
mv CLUBB/src/CLUBB_core CLUBB_core

echo "removing .svn folder from CLUBB_core"
rm -rf CLUBB_core/.svn

echo "testing CLUBB_standalone"
python api_commitment_test.py -cpu CLUBB_core CLUBB > clubb_standalone_modules.txt

echo "testing CLUBB_core"
python api_commitment_test.py -cpu CLUBB_core CLUBB_core > clubb_core_modules.txt

echo "testing SAM"
python api_commitment_test.py -cpu CLUBB_core SAM --exclude-dir="CLUBB","SILHS" > sam_modules.txt

echo "testing CAM"
python api_commitment_test.py -cpu CLUBB_core CAM --exclude-dir="clubb","silhs" > cam_modules.txt

echo "testing WRF"
python api_commitment_test.py -cpu CLUBB_core WRF --exclude-dir="clubb","silhs" > wrf_modules.txt

python create_module_table.py CLUBB_core> output.html

