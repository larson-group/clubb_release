#!/bin/bash

############################################################
# This script is meant to make recreatign the storer
# simulations as easy as possible. Simply run
# this script without any parameters.
#
# Author: Nicolas Strike - Jun 2019
############################################################

SAMfiles=( ARM9707.nc TWPICE.nc LBA.nc RICO.nc DYCOMS_RF02_DO.nc )


checkForFile(){
    FILE=$1
    if [ -f "SAMoutput/$FILE" ]; then
        echo "Sam output file $FILE found."
    else 
        echo "Sam output file $FILE not found. Downloading now."
        cd SAMoutput
        wget https://carson.math.uwm.edu/larson-group/storer/Storeretal_SAM_output/$FILE
        cd ..
    fi
}

for file in "${SAMfiles[@]}"
do
    checkForFile $file
done

wait
bash compile/compile.bash
wait
bash run_scripts/run_scm_all.bash
