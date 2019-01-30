#!/bin/bash
# $Id$
#-------------------------------------------------------------------------------
# Description:
# Script to run G_unit_tests, which is the unit testing framework.

# Figure out the directory where the script is located.
run_dir=`dirname $0`

# Change directories to the one in which the script is located.
cd $run_dir

# Copy the namelist to the local directory.
cp "../input_misc/G_unit_tests.in" "G_unit_tests.in"

# Run G_unit_tests
./../bin/G_unit_tests

if [ $? != 0 ]; then
  # Remove the namelist from the local directory.
  rm "G_unit_tests.in"
  printf "Error: A G unit test has failed.\n"
  exit 1
else
  # Remove the namelist from the local directory.
  rm "G_unit_tests.in"
  printf "All tests have succeeded.\n"
  exit 0
fi
