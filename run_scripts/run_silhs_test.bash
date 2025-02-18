#!/bin/bash
#################################################################################
# $Id$
#
# Desciption:
#   Script to run the test whether SILHS is still converging. It runs RUN_CASE
#   twice once with N_SMALL number of sample points and once with N_LARGE. Then
#   it calculates the root mean squared errors (RMSE) to see if lh_AKm converges
#   to AKm.
#################################################################################

# Set case to run
RUN_CASE="rico"
RUN_CASE_INPUT="../input/case_setups/${RUN_CASE}_model.in"

# Set the number of sample points for the two runs
N_SMALL=8
N_LARGE=1000

# Define output directories for both runs
OUTPUT_SMALL="../output/small"
OUTPUT_LARGE="../output/large"

# Figure out the directory where the script is located
scriptPath=`dirname $0`

# Store the current directory location so it can be restored
restoreDir=`pwd`

# Change directories to the one the script is located in
cd $scriptPath

# Copy the original input file so it can be restored
cp $RUN_CASE_INPUT ../input/case_setups/tmp_model.in

# Get line number after which to input the SILHS relevant settings
LINNUM="$(grep -n "&microphysics_setting" $RUN_CASE_INPUT | head -n 1 | cut -d: -f1)"

sed -i "$LINNUM a lh_microphys_type = \"non-interactive\"\nlh_num_samples = 0\nlh_sequence_length = 1\nl_silhs_KK_convergence_adj_mean = .true.\nl_local_kk = .false." $RUN_CASE_INPUT

echo -n "Running $RUN_CASE case with small number of sample points... "
sed -i "/lh_num_samples = .*/c\lh_num_samples = $N_SMALL" $RUN_CASE_INPUT
./run_scm.bash --netcdf $RUN_CASE &> /dev/null
echo "Done!"

mkdir $OUTPUT_SMALL
mv ../output/*.* $OUTPUT_SMALL

echo -n "Running $RUN_CASE case with large number of sample points... "
sed -i "/lh_num_samples = .*/c\lh_num_samples = $N_LARGE" $RUN_CASE_INPUT
./run_scm.bash --netcdf $RUN_CASE &> /dev/null
echo "Done!"

mkdir $OUTPUT_LARGE
mv ../output/*.* $OUTPUT_LARGE

export RUN_CASE N_SMALL N_LARGE OUTPUT_SMALL OUTPUT_LARGE

# Calculate RMSEs
python3 -c '
from netCDF4 import Dataset;
import numpy as np;
import os;
import sys;
try:
    f = Dataset(os.environ["OUTPUT_SMALL"]+"/"+os.environ["RUN_CASE"]+"_lh_zt.nc", mode="r");
    AKm = np.asarray(np.asmatrix(f.variables["AKm"][:]));
    lh_AKm = np.asarray(np.asmatrix(f.variables["lh_AKm"][:]));
    f.close();
except:
    sys.exc_info()[0];
    sys.exit(1);
if AKm.shape!=lh_AKm.shape or len(AKm.shape)!=2 or AKm.shape[0]==0 or AKm.shape[1]==0:
    sys.exit(1);
shape_small = AKm.shape;
RMSE_small = np.sqrt(1/AKm.shape[1]*sum([(np.mean(lh_AKm[:, k])-np.mean(AKm[:, k]))**2 for k in range(AKm.shape[1])]));
try:
    f = Dataset(os.environ["OUTPUT_LARGE"]+"/"+os.environ["RUN_CASE"]+"_lh_zt.nc", mode="r");
    AKm = np.asarray(np.asmatrix(f.variables["AKm"][:]));
    lh_AKm = np.asarray(np.asmatrix(f.variables["lh_AKm"][:]));
    f.close();
except:
    sys.exc_info()[0];
    sys.exit(1);
if AKm.shape!=lh_AKm.shape or len(AKm.shape)!=2 or AKm.shape[0]==0 or AKm.shape[1]==0:
    sys.exit(1);
shape_large = AKm.shape;
RMSE_large = np.sqrt(1/AKm.shape[1]*sum([(np.mean(lh_AKm[:, k])-np.mean(AKm[:, k]))**2 for k in range(AKm.shape[1])]));
if shape_small!=shape_large:
    sys.exit(1);
if RMSE_small<0 or RMSE_large<0 or RMSE_small>1e-5 or RMSE_large>1e-5:
    sys.exit(1);
print("N_small =", os.environ["N_SMALL"]);
print("N_large =", os.environ["N_LARGE"]);
print("RMSE_small =", RMSE_small);
print("RMSE_large =", RMSE_large);
N_SMALL = int(os.environ["N_SMALL"])
N_LARGE = int(os.environ["N_LARGE"])
if RMSE_large/RMSE_small < 4*np.sqrt(N_SMALL)/np.sqrt(N_LARGE):
    sys.exit(0);
else:
    print(f"FAIL: RMSE_large/RMSE_small = {RMSE_large/RMSE_small} is greater than 4*sqrt(N_SMALL)/sqrt(N_LARGE) = {4*np.sqrt(N_SMALL)/np.sqrt(N_LARGE)}")
    sys.exit(1)
'

if [ $? -ne 0 ]; then
  echo "SILHS failed to converge!"
  RESULT=1
else
  echo "SILHS converged!"
  RESULT=0
fi

mv ../input/case_setups/tmp_model.in $RUN_CASE_INPUT
rm -rf $OUTPUT_SMALL
rm -rf $OUTPUT_LARGE

cd $restoreDir

exit $RESULT
