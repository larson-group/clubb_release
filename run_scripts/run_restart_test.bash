#!/bin/bash
#################################################################################
# $Id$
#
# Desciption:
#   Script to test whether CLUBB produces bit-for-bit identical results at
#   the last time step when restarted from halfway through a test case.
#   First runs the full test case, then restarts from halfway through and 
#   compares results from a specified variable.
#################################################################################

# Set case to run and choose variable to test
RUN_CASE="bomex"
VAR_TO_TEST="thlm"

# Specify input file
RUN_CASE_INPUT="../input/case_setups/${RUN_CASE}_model.in"

# Figure out the directory where the script is located
scriptPath=`dirname $0`

# Store the current directory location so it can be restored
restoreDir=`pwd`

# Change directories to the one the script is located in
cd $scriptPath

# Copy the original input file so it can be restored
cp $RUN_CASE_INPUT ../input/case_setups/tmp_model.in

# Figure out the halfway point (also works for cases that don't start at t=0)
sed -i 's/!time_initial\s*=\s*.*//g' $RUN_CASE_INPUT
sed -i 's/!time_final\s*=\s*.*//g' $RUN_CASE_INPUT
time_init=$(grep -Po 'time_initial\s*= \K\d+\.*\d*' $RUN_CASE_INPUT)
time_final=$(grep -Po 'time_final\s*= \K\d+\.*\d*' $RUN_CASE_INPUT)
sum="$( bc <<<"$time_init + $time_final" )"
halfway_time="$( bc <<< "scale=1; $sum / 2")"

# run standard case
echo -n "Running standard $RUN_CASE case... "
./run_scm.bash $RUN_CASE &> /dev/null
echo "Done!"

# Create restart folder and move 90-min run there
mkdir ../restart
mv ../output/*.* ../restart/

# Set case file for restart runs
sed -i -e 's/time_initial\s*=\s*.*/time_initial = '$halfway_time'/g' \
       -e 's/l_restart\s*=\s*.*/l_restart = .true./g' \
       -e 's/restart_path_case\s*=\s*.*/restart_path_case = "restart\/'$RUN_CASE'"/g' \
       -e 's/time_restart\s*=\s*.*/time_restart = '$halfway_time'/g' $RUN_CASE_INPUT

# run restart case
echo -n "Running restart $RUN_CASE case from halfway point... "
./run_scm.bash $RUN_CASE &> /dev/null
echo "Done!"

#export variables to system environment for python use
export RUN_CASE VAR_TO_TEST

python3 -c '
from netCDF4 import Dataset;
import numpy as np;
import os;
import sys;
dataset1=Dataset("../restart/"+os.environ["RUN_CASE"]+"_zt.nc")
dataset2=Dataset("../output/"+os.environ["RUN_CASE"]+"_zt.nc")
var_to_test1=np.mean(np.array(dataset1.variables[os.environ["VAR_TO_TEST"]]),(2,3))
var_to_test2=np.mean(np.array(dataset2.variables[os.environ["VAR_TO_TEST"]]),(2,3))
diff = var_to_test1[-1,:] - var_to_test2[-1,:]
if all(diff) == 0.0:
    sys.exit(0);
else:
    sys.exit(1)
'

if [ $? -ne 0 ]; then
  echo "Results not bit-for-bit!"
  RESULT=1
else
  echo "Results bit-for-bit."
  RESULT=0
fi

# replace original model.in file
mv ../input/case_setups/tmp_model.in $RUN_CASE_INPUT

# cleanup the output
rm -rf ../restart
rm -rf ../output/*

cd $restoreDir

exit $RESULT



