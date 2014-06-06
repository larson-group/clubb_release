#!/bin/bash
# $Id$
#-------------------------------------------------------------------------------
# Description:
# Script to run G_unit_tests, which is the unit testing framework.

# Copy the namelist to the local directory.
cp "../input_misc/G_unit_tests.in" "G_unit_tests.in"

# Run G_unit_tests
../bin/G_unit_tests

# Remove the namelist from the local directory.
rm "G_unit_tests.in"
