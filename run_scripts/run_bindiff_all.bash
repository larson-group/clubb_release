#!/bin/bash
# $Id$
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
# then the user runs BOMEX twice -- once with standard_stats.in and once with 
# nobudgets_stats.in -- the control and binary data files will be different even
# though it is the exact same code.
#
# This script loops over all cases found in run_scm_all.bash.  This 
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
#
# Ryan Senkbeil added if statements to check if files exist before diffing them
# July 6, 2009.
#-------------------------------------------------------------------------------

# Figure out the directory where the script is located
scriptPath=`dirname $0`

# Store the current directory location so it can be restored
restoreDir=`pwd`

# Change directories to the one the script is located in
cd $scriptPath

# If true, differences were detected. This is used to set the exit status
differences=false

declare -a RUN_CASE

a=0
while read line
do
    # If the line is not commented out (does not start with '!')
    if [[ $line != !* ]] && [[ ! -z $line ]] ; then
        RUN_CASE[$a]=$line
        a=$(($a+1));
    fi
done < "RUN_CASES"

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
dir1=$restoreDir/$1

# The second command-line entry is 'dir2'.
dir2=$restoreDir/$2

# diff the GrADS control (*.ctl) files and the GrADS binary data (*.dat) files
# for each statistical output type (zt, zm, and sfc) for each case in 
# run_standalone-all.bash.
# This will loop over all runs in sequence.
for (( x=0; x < "${#RUN_CASE[@]}"; x++ )); do

    # State which case is being diffed.
    echo 'Diffing '"${RUN_CASE[$x]}"' GrADS control (*.ctl) and binary data (*.dat) files'

    if [ -e $dir1/"${RUN_CASE[$x]}"'_zt.ctl' -a -e $dir2/"${RUN_CASE[$x]}"'_zt.ctl' ] ; then
        # Compare the zt GrADS control (*_zt.ctl) files
        diffZtCtl=$(diff $dir1/"${RUN_CASE[$x]}"'_zt.ctl' $dir2/"${RUN_CASE[$x]}"'_zt.ctl')

        if [ -n "$diffZtCtl" ] ; then
            differences=true
            echo "*** Differences detected in ${RUN_CASE[$x]}_zt.ctl! ***" >&2
        fi
    fi

    if [ -e $dir1/"${RUN_CASE[$x]}"'_zt.dat' -a -e $dir2/"${RUN_CASE[$x]}"'_zt.dat' ] ; then
        # Compare the zt GrADS binary data (*_zt.dat) files
        diffZtDat=$(diff $dir1/"${RUN_CASE[$x]}"'_zt.dat' $dir2/"${RUN_CASE[$x]}"'_zt.dat')

        if [ -n "$diffZtDat" ] ; then
            differences=true
            echo "*** Differences detected in ${RUN_CASE[$x]}_zt.dat! ***" >&2
        fi
    fi

    if [ -e $dir1/"${RUN_CASE[$x]}"'_zm.ctl' -a -e $dir2/"${RUN_CASE[$x]}"'_zm.ctl' ] ; then
        # Compare the zm GrADS control (*_zm.ctl) files
        diffZmCtl=$(diff $dir1/"${RUN_CASE[$x]}"'_zm.ctl' $dir2/"${RUN_CASE[$x]}"'_zm.ctl')

        if [ -n "$diffZmCtl" ] ; then
            differences=true
            echo "*** Differences detected in ${RUN_CASE[$x]}_zm.ctl! ***" >&2
        fi
    fi

    if [ -e $dir1/"${RUN_CASE[$x]}"'_zm.dat' -a -e $dir2/"${RUN_CASE[$x]}"'_zm.dat' ] ; then
        # Compare the zm GrADS binary data (*_zm.dat) files
        diffZmDat=$(diff $dir1/"${RUN_CASE[$x]}"'_zm.dat' $dir2/"${RUN_CASE[$x]}"'_zm.dat')

        if [ -n "$diffZmDat" ] ; then
            differences=true
            echo "*** Differences detected in ${RUN_CASE[$x]}_zm.dat! ***" >&2
        fi
    fi

    if [ -e $dir1/"${RUN_CASE[$x]}"'_sfc.ctl' -a -e $dir2/"${RUN_CASE[$x]}"'_sfc.ctl' ] ; then
        # Compare the sfc GrADS control (*_sfc.ctl) files
        diffSfcCtl=$(diff $dir1/"${RUN_CASE[$x]}"'_sfc.ctl' $dir2/"${RUN_CASE[$x]}"'_sfc.ctl')

        if [ -n "$diffSfcCtl" ] ; then
            differences=true
            echo "*** Differences detected in ${RUN_CASE[$x]}_sfc.ctl! ***" >&2
        fi
    fi

    if [ -e $dir1/"${RUN_CASE[$x]}"'_sfc.dat' -a -e $dir2/"${RUN_CASE[$x]}"'_sfc.dat' ] ; then
        # Compare the sfc GrADS binary data (*_sfc.dat) files
        diffSfcDat=$(diff $dir1/"${RUN_CASE[$x]}"'_sfc.dat' $dir2/"${RUN_CASE[$x]}"'_sfc.dat')

        if [ -n "$diffSfcDat" ] ; then
            differences=true
            echo "*** Differences detected in ${RUN_CASE[$x]}_sfc.dat! ***" >&2
        fi
    fi
done

# Determine exit status and exit
if [ $differences == "true" ] ; then
    echo -e "\nThere were some differences detected!"
    exit 1
else
    echo -e "\nThere were no differences detected!"
    exit 0
fi

cd $restoreDir
