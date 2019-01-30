#!/bin/bash
#######################################################################
# $Id$
#
# Script to test CLUBB's convergence by running a specified case at various time step
# lengths.  An example call:
#
# ./timestep_convergence_test.bash bomex
#
#######################################################################

# Figure out the directory where the script is located
scriptPath=`dirname $0`

# Store the current directory location so it can be restored
restoreDir=`pwd`

# Change directories to the one the script is located in
cd $scriptPath

model_file='../input/case_setups/'${!#}'_model.in'

TEST_TIMESTEP[0]=300      # Time step of output to disk 
                                 # and longest computational time step tested
TEST_TIMESTEP[1]=150       # This and all subsequent time steps
TEST_TIMESTEP[2]=100       # should divide evenly into TEST_TIMESTEP[0].
TEST_TIMESTEP[3]=60  # 20 minute time step.
TEST_TIMESTEP[4]=30  # 30 minute time step.
TEST_TIMESTEP[5]=20  # 45 minute time step.
TEST_TIMESTEP[6]=10  # 45 minute time step.
TEST_TIMESTEP[7]=5  # 60 minute (one hour) time step.
TEST_TIMESTEP[8]=3  # 60 minute (one hour) time step.
TEST_TIMESTEP[9]=2  # 60 minute (one hour) time step.
TEST_TIMESTEP[10]=1  # 60 minute (one hour) time step.
TEST_TIMESTEP[11]=0.5  # 60 minute (one hour) time step.
TEST_TIMESTEP[12]=0.3  # 60 minute (one hour) time step.
TEST_TIMESTEP[13]=0.2  # 60 minute (one hour) time step.
TEST_TIMESTEP[14]=0.1  # 60 minute (one hour) time step.

# Set the time interval of disk output
sed -i -e 's/stats_tsamp\s*=\s*.*/stats_tsamp = '${TEST_TIMESTEP[0]}'/g' \
       -e 's/stats_tout\s*=\s*.*/stats_tout = '${TEST_TIMESTEP[0]}'/g' $model_file
echo -e "\nOutputting to disk at a "${TEST_TIMESTEP[0]}"-second time step."

for (( x=0; x < "${#TEST_TIMESTEP[@]}"; x++ )); do
   echo -e "\nRunning a case at a "${TEST_TIMESTEP[$x]}"-second time step."

   # Set computational time step
   sed -i -e 's/dt_main\s*=\s*.*/dt_main = '${TEST_TIMESTEP[$x]}' /g' \
          -e 's/dt_rad\s*=\s*.*/dt_rad =   '${TEST_TIMESTEP[$x]}' /g' \
          $model_file
   # Convert decimal portion so that we can use time step in the output file names
   #integer_timestep=${TEST_TIMESTEP[$x]%.*}
   integer_timestep=${TEST_TIMESTEP[$x]/./p}
   sed -i -r 's/fname_prefix\s*=\s*"(.*?)"/fname_prefix = "'${!#}'_'$integer_timestep'" /g' \
          $model_file

   ../run_scripts/run_scm.bash --netcdf ${!#} 2>&1 >/dev/null
done

cd $restoreDir
