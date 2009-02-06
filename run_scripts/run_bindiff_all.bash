#!/bin/bash
#-------------------------------------------------------------------------------
#
# run_bindiff-all.bash
# --------------------
#
# Description:
# This bash script is written to check for any binary differences between two
# sets of GrADS binary data (*.dat) files.
#
# This bash script also checks for any differences between the GrADS 
# control (*.ctl) files; for if the control files contain different 
# amounts/orders of variables, different amounts of altitudes, or different 
# amounts of time outputs, the binary data files will differ, regardless if the 
# case results are the same.  For example, if the CLUBB code is compiled; and 
# then the user runs BOMEX twice -- once with all_stats.in and once with 
# nobudgets_stats.in -- the control and binary data files will be different even
# though it is the exact same code.
#
# This script loops over all cases found in run_standalone-all.bash.  This 
# script is useful for testing if CLUBB code changes made any binary difference 
# in output results for any case.  For example, simply changing a variable name 
# in the code or a changing subroutine name shouldn't result in any binary 
# difference for any of the cases.  This script can also be used to check if a 
# CLUBB code change made tiny differences in the results.  This is useful if the
# differences are too small to be seen on plotgen.
#
# In order to use this script, two directories that contain GrADS control files
# and binary data files must be entered on the command line.
#
# Brian Griffin; August 23, 2008.
#-------------------------------------------------------------------------------


# Loop over all cases that are found in run_standalone-all.bash.
# Note:  The cases listed in RUN_CASE must always match the cases found in 
#        RUN_CASE in run_standalone-all.bash.
RUN_CASE=( \
	arm arm_97 atex bomex clex9_nov02 clex9_oct14 cobra dycoms2_rf01
        dycoms2_rf02_do dycoms2_rf02_ds	dycoms2_rf02_nd dycoms2_rf02_so \
        fire gabls2 jun25_altocu lba mpace_a mpace_b nov11_altocu rico wangara )


# The user needs to enter the paths/names for two directories on the command line.
if [ -z $1 ]; then
   echo 'Two directories need to be entered on the command line.'
   exit
elif [ -z $2 ]; then
   echo 'Directory 1 is '$1'.'
   echo 'Two directories need to be entered on the command line.'
   exit
else
   echo 'Directory 1 is '$1'.'
   echo 'Directory 2 is '$2'.'
fi

# The first command-line entry is 'dir1'.
dir1=$1

# The second command-line entry is 'dir2'.
dir2=$2


# diff the GrADS control (*.ctl) files and the GrADS binary data (*.dat) files
# for each statistical output type (zt, zm, and sfc) for each case in 
# run_standalone-all.bash.
# This will loop over all runs in sequence.
for (( x=0; x < "${#RUN_CASE[@]}"; x++ )); do


   # State which case is being diffed.
   echo 'Diffing '"${RUN_CASE[$x]}"' GrADS control (*.ctl) and binary data (*.dat) files'


   # Compare the zt GrADS control (*_zt.ctl) files
   diff $dir1/"${RUN_CASE[$x]}"'_zt.ctl' $dir2/"${RUN_CASE[$x]}"'_zt.ctl'

   # Compare the zt GrADS binary data (*_zt.dat) files
   diff $dir1/"${RUN_CASE[$x]}"'_zt.dat' $dir2/"${RUN_CASE[$x]}"'_zt.dat'


   # Compare the zm GrADS control (*_zm.ctl) files
   diff $dir1/"${RUN_CASE[$x]}"'_zm.ctl' $dir2/"${RUN_CASE[$x]}"'_zm.ctl'

   # Compare the zm GrADS binary data (*_zm.dat) files
   diff $dir1/"${RUN_CASE[$x]}"'_zm.dat' $dir2/"${RUN_CASE[$x]}"'_zm.dat'


   # Compare the sfc GrADS control (*_sfc.ctl) files
   diff $dir1/"${RUN_CASE[$x]}"'_sfc.ctl' $dir2/"${RUN_CASE[$x]}"'_sfc.ctl'

   # Compare the sfc GrADS binary data (*_sfc.dat) files
   diff $dir1/"${RUN_CASE[$x]}"'_sfc.dat' $dir2/"${RUN_CASE[$x]}"'_sfc.dat'


done
