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

case_name=${!#}

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

echo -e "\nOutputting to disk at a "${TEST_TIMESTEP[0]}"-second time step."

for (( x=0; x < "${#TEST_TIMESTEP[@]}"; x++ )); do
   echo -e "\nRunning a case at a "${TEST_TIMESTEP[$x]}"-second time step."

   integer_timestep=${TEST_TIMESTEP[$x]/./p}
   python3 ../run_scripts/run_scm.py "$case_name" -override \
       "dt_main=${TEST_TIMESTEP[$x]},dt_rad=${TEST_TIMESTEP[$x]},stats_tsamp=${TEST_TIMESTEP[0]},stats_tout=${TEST_TIMESTEP[0]},fname_prefix='${case_name}_${integer_timestep}',stats_output_filename='${case_name}_${integer_timestep}_stats.nc'"
done

cd $restoreDir
