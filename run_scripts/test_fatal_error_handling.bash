#!/bin/bash

###############################################################################
# $Id$
#
#   Description:
#       This script injects errors into clubb then runs a case and captures the 
#       error output. The injected errors need to be labeled in the fortran code
#       with a name, then added to the testName array in this script along with
#       the name of the variable to be tested to the testVariable array.
#
#   Example:
#
#       In fortran code:    
#                           ! advance_wp2_wp3_bad_wp2
#
#       In this script:
#                           testName[1]="advance_wp2_wp3_bad_wp2"
#                           testVariable[1]="wp2"
#
#       This results in "! advance_wp2_wp3_bad_wp2" in the code being replaced
#       with a line that sets "wp2" to NaN, then the error messages that clubb
#       produces are saved in a file named "advance_wp2_wp3_bad_wp2" in run_scripts
#
#   Notes:
#           THIS SCRIPT NEEDS TO BE RUN FROM WITHIN 'run_scripts' OR MODIFIED TO
#           WORK IN A DIFFERENT DIRECTORY
#
#       To add error tests, simply add the lines that are explained in the example
#       above, changing the test name to something unique andthe variable to the one
#       that will be tested. Remeber to change the index as well. 
#
###############################################################################


declare -a testName
declare -a testVariable

########### TEST DEFINITIONS  ##################
testName[1]="advance_wp2_wp3_bad_wp2"
testVariable[1]="wp2"

testName[2]="advance_xm_wpxp_bad_wp2"
testVariable[2]="wp2"

# Save compiler flags so they can be restored after the test
initial_debug_flags=$(grep "^DEBUG=.*" ../compile/config/linux_x86_64_gfortran.bash)

# Change compiler flags so the compliation ignores the NaN, and any floating point errors
sed -i "s:^DEBUG=\":DEBUG=\"-fno-range-check :g " ../compile/config/linux_x86_64_gfortran.bash
sed -i "s:-ffpe-trap=invalid,zero,overflow::g " ../compile/config/linux_x86_64_gfortran.bash

for(( j=1; j<=${#testName[@]}; j++)); do

    # Find test name and replace with line that sets test variable to NaN
    sed -i "s:! ${testName[j]}:${testVariable[j]} = transfer( 2143289344, 1.0 ):g" ../src/CLUBB_core/advance_clubb_core_module.F90

    #compile
    ./../compile/compile.bash

    # Run a case (any non error causing case will work) and write error output to a file
    ./run_scm.bash args = "arm_97" > /dev/null 2> ${testName[j]}

    # Undo previous sed 
    sed -i "s:${testVariable[j]} = transfer( 2143289344, 1.0 ):! ${testName[j]}:g" ../src/CLUBB_core/advance_clubb_core_module.F90

done

# Restore compiler flags
sed -i "s:^DEBUG=.*:$initial_debug_flags:g" ../compile/config/linux_x86_64_gfortran.bash
