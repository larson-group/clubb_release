#!/bin/csh

#$ -pe ic.postp 2
#$ -l h_cpu=10:00:00
#$ -N sort
#$ -o /archive/cjg/calibration/revised/exp/post.log

# Start bash scripts

cd $HOME/archive

echo Sorting results
./sortresults.bash

echo Averaging best simulations
./average.bash

# Done

exit 0
