#!/bin/bash
#######################################################################
# $Id$
#
# Script to test CLUBB by running all cases at various time step
# lengths.  This script calls:
#
# ./run_standalone_all.bash --timestep_test {time step length}
#
# multiple times, using a different time step length each time.
#
#######################################################################

TEST_TIMESTEP[0]=300.0   # 5 minute time step.
TEST_TIMESTEP[1]=600.0   # 10 minute time step.
TEST_TIMESTEP[2]=900.0   # 15 minute time step.
TEST_TIMESTEP[3]=1200.0  # 20 minute time step.
TEST_TIMESTEP[4]=1800.0  # 30 minute time step.
TEST_TIMESTEP[5]=2700.0  # 45 minute time step.
TEST_TIMESTEP[6]=3600.0  # 60 minute (one hour) time step.

for (( x=0; x < "${#TEST_TIMESTEP[@]}"; x++ )); do

   ./run_standalone_all.bash --timestep_test ${TEST_TIMESTEP[$x]}

done
