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
#
# Kenneth Connor modified the script so the output directories don't need to be
# subdirectories of the location you called the script from. Also added checks
# to confirm that the directories exist and if they have any data. The script
# now only outputs that it is checking a case if the data for that case exists.
# June 17, 2011.
#
# Brian Griffin upgraded this script to compare netCDF output files.
# June 5, 2014.
#-------------------------------------------------------------------------------

# Figure out the directory where the script is located
scriptPath=`dirname $0`

# Store the current directory location so it can be restored
restoreDir=`pwd`

# If true, differences were detected. This is used to set the exit status
differences=false
differencesCtl=false
differencesDat=false
differencesNc=false

# If true, there was no data found. This is used to set the exit status
noData=true

# The user needs to enter the paths/names for two directories on the command line.
if [ -z $1 ]; then
    echo 'Two directories need to be entered on the command line.'
    exit
elif [ -z $2 ]; then
    echo 'Directory 1 is '$(readlink -f "$1")'.'
    echo 'Two directories need to be entered on the command line.'
    exit
else
    # The first command-line entry is 'dir1'.
    dir1=$(readlink -f "$1")

    # The second command-line entry is 'dir2'.
    dir2=$(readlink -f "$2")
    
    echo 'Directory 1 is '$dir1'.'
    echo 'Directory 2 is '$dir2'.'

    if [ ! -d "$dir1" ]; then
      echo $dir1' does not exist or is not a directory.'
      exit
    elif [ ! -d "$dir2" ]; then
      echo $dir2' does not exist or is not a directory.'
      exit
    fi
fi

# Change directories to the one the script is located in
cd $scriptPath

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

# diff the GrADS control (*.ctl) files and the GrADS binary data (*.dat) files,
# or netCDF (*.nc) files, for each statistical output type (zt, zm, and sfc) for
# each case in found in both directories that are part of RUN_CASES.
# This will loop over all runs in sequence.
for (( x=0; x < "${#RUN_CASE[@]}"; x++ )); do

    declare -a DIFF_LIST

    dataFound=false
    dataFoundGrADS=false
    dataFoundnetCDF=false

    if [ -e $dir1/"${RUN_CASE[$x]}"'_zt.ctl' -a -e $dir2/"${RUN_CASE[$x]}"'_zt.ctl' ] ; then
        # Compare the zt GrADS control (*_zt.ctl) files
        diffZtCtl=$(diff $dir1/"${RUN_CASE[$x]}"'_zt.ctl' $dir2/"${RUN_CASE[$x]}"'_zt.ctl')

        if [ -n "$diffZtCtl" ] ; then
            differences=true
            differencesCtl=true
            DIFF_LIST[${#DIFF_LIST[@]}]="${RUN_CASE[$x]}_zt.ctl"
        fi
        
        dataFound=true
        dataFoundGrADS=true
    fi

    if [ -e $dir1/"${RUN_CASE[$x]}"'_zt.dat' -a -e $dir2/"${RUN_CASE[$x]}"'_zt.dat' ] ; then
        # Compare the zt GrADS binary data (*_zt.dat) files
        diffZtDat=$(diff $dir1/"${RUN_CASE[$x]}"'_zt.dat' $dir2/"${RUN_CASE[$x]}"'_zt.dat')

        if [ -n "$diffZtDat" ] ; then
            differences=true
            differencesDat=true
            DIFF_LIST[${#DIFF_LIST[@]}]="${RUN_CASE[$x]}_zt.dat"
        fi
        
        dataFound=true
        dataFoundGrADS=true
    fi

    if [ -e $dir1/"${RUN_CASE[$x]}"'_zt.nc' -a -e $dir2/"${RUN_CASE[$x]}"'_zt.nc' ] ; then
        # Compare the zt netCDF (*_zt.nc) files
        diffZtNc=$(python diff_netcdf_outputs.py $dir1/"${RUN_CASE[$x]}"'_zt.nc' $dir2/"${RUN_CASE[$x]}"'_zt.nc')

        if [ -n "$diffZtNc" ] ; then
            differences=true
            differencesNc=true
            DIFF_LIST[${#DIFF_LIST[@]}]="${RUN_CASE[$x]}_zt.nc"
        fi
        
        dataFound=true
        dataFoundnetCDF=true
    fi

    if [ -e $dir1/"${RUN_CASE[$x]}"'_zm.ctl' -a -e $dir2/"${RUN_CASE[$x]}"'_zm.ctl' ] ; then
        # Compare the zm GrADS control (*_zm.ctl) files
        diffZmCtl=$(diff $dir1/"${RUN_CASE[$x]}"'_zm.ctl' $dir2/"${RUN_CASE[$x]}"'_zm.ctl')

        if [ -n "$diffZmCtl" ] ; then
            differences=true
            differencesCtl=true
            DIFF_LIST[${#DIFF_LIST[@]}]="${RUN_CASE[$x]}_zm.ctl"
        fi
        
        dataFound=true
        dataFoundGrADS=true
    fi

    if [ -e $dir1/"${RUN_CASE[$x]}"'_zm.dat' -a -e $dir2/"${RUN_CASE[$x]}"'_zm.dat' ] ; then
        # Compare the zm GrADS binary data (*_zm.dat) files
        diffZmDat=$(diff $dir1/"${RUN_CASE[$x]}"'_zm.dat' $dir2/"${RUN_CASE[$x]}"'_zm.dat')

        if [ -n "$diffZmDat" ] ; then
            differences=true
            differencesDat=true
            DIFF_LIST[${#DIFF_LIST[@]}]="${RUN_CASE[$x]}_zm.dat"
        fi
        
        dataFound=true
        dataFoundGrADS=true
    fi

    if [ -e $dir1/"${RUN_CASE[$x]}"'_zm.nc' -a -e $dir2/"${RUN_CASE[$x]}"'_zm.nc' ] ; then
        # Compare the zm netCDF (*_zm.nc) files
        diffZmNc=$(python diff_netcdf_outputs.py $dir1/"${RUN_CASE[$x]}"'_zm.nc' $dir2/"${RUN_CASE[$x]}"'_zm.nc')

        if [ -n "$diffZmNc" ] ; then
            differences=true
            differencesNc=true
            DIFF_LIST[${#DIFF_LIST[@]}]="${RUN_CASE[$x]}_zm.nc"
        fi
        
        dataFound=true
        dataFoundnetCDF=true
    fi

    if [ -e $dir1/"${RUN_CASE[$x]}"'_sfc.ctl' -a -e $dir2/"${RUN_CASE[$x]}"'_sfc.ctl' ] ; then
        # Compare the sfc GrADS control (*_sfc.ctl) files
        diffSfcCtl=$(diff $dir1/"${RUN_CASE[$x]}"'_sfc.ctl' $dir2/"${RUN_CASE[$x]}"'_sfc.ctl')

        if [ -n "$diffSfcCtl" ] ; then
            differences=true
            differencesCtl=true
            DIFF_LIST[${#DIFF_LIST[@]}]="${RUN_CASE[$x]}_sfc.ctl"
        fi
        
        dataFound=true
        dataFoundGrADS=true
    fi

    if [ -e $dir1/"${RUN_CASE[$x]}"'_sfc.dat' -a -e $dir2/"${RUN_CASE[$x]}"'_sfc.dat' ] ; then
        # Compare the sfc GrADS binary data (*_sfc.dat) files
        diffSfcDat=$(diff $dir1/"${RUN_CASE[$x]}"'_sfc.dat' $dir2/"${RUN_CASE[$x]}"'_sfc.dat')

        if [ -n "$diffSfcDat" ] ; then
            differences=true
            differencesDat=true
            DIFF_LIST[${#DIFF_LIST[@]}]="${RUN_CASE[$x]}_sfc.dat"
        fi
        
        dataFound=true
        dataFoundGrADS=true
    fi

    if [ -e $dir1/"${RUN_CASE[$x]}"'_sfc.nc' -a -e $dir2/"${RUN_CASE[$x]}"'_sfc.nc' ] ; then
        # Compare the sfc netCDF (*_sfc.nc) files
        diffSfcNc=$(python diff_netcdf_outputs.py $dir1/"${RUN_CASE[$x]}"'_sfc.nc' $dir2/"${RUN_CASE[$x]}"'_sfc.nc')

        if [ -n "$diffSfcNc" ] ; then
            differences=true
            differencesNc=true
            DIFF_LIST[${#DIFF_LIST[@]}]="${RUN_CASE[$x]}_sfc.nc"
        fi
        
        dataFound=true
        dataFoundnetCDF=true
    fi

    if [ $dataFound == "true" ] ; then
        if [ $dataFoundGrADS == "true" ] ; then
            echo 'Diffing '"${RUN_CASE[$x]}"' GrADS control (*.ctl) and binary data (*.dat) files'
        fi
        if [ $dataFoundnetCDF == "true" ] ; then
            echo 'Diffing '"${RUN_CASE[$x]}"' netCDF (*.nc) files'
        fi

      for (( y=0; y < "${#DIFF_LIST[@]}"; y++ )); do
        echo "*** Differences detected in ${DIFF_LIST[$y]}! ***" >&2
      done

      noData=false
    fi
    
    # Clear the array
    DIFF_LIST=( )

done

# Determine exit status and exit
if [ $differences == "true" ] ; then
    echo -e "\nThere were some differences detected!"
    if [ $differencesCtl == "true" ] ; then
       echo -e "There were differences detected in GrADS control (*.ctl) files."
    fi
    if [ $differencesDat == "true" ] ; then
       echo -e "There were differences detected in GrADS binary data (*.dat) files."
    fi
    if [ $differencesNc == "true" ] ; then
       echo -e "There were differences detected in netCDF (*.nc) files."
    fi
    exit 1
elif [ $noData == "true" ] ; then
    echo -e "\nNo data was found!"
    exit 1
else
    echo -e "\nThere were no differences detected!"
    exit 0
fi

cd $restoreDir
