#!/bin/bash
# $Id: generate_seed.bash 5197 2011-05-24 20:09:18Z dschanen@uwm.edu $
# Usage:
#  Create a new rand_seed.dat 
#  Used by clubb_tuner to seed the Merssene twister pseudo-random number generator.
#----------------------------------------------------------------------
# Output:
# rand_seed.dat: text file that may be edited by hand.
#----------------------------------------------------------------------
SEED="../input_misc/tuner/rand_seed.dat"
echo "Generating a new seed at: " $SEED
dd if=/dev/urandom of=rand_seed.bin bs=4 count=34 &> /dev/null
../bin/int2txt rand_seed.bin > $SEED
rm -f rand_seed.bin
