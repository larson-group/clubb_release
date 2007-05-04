#!/bin/bash
# $Id: sortresults.bash,v 1.2 2007-05-04 18:08:05 dschanen Exp $
################################################################################
# Description:
# Search the log files in each subdirectory of an ensemble tuning run and then
# sort the results from lowest J (cost function) to the highest.
# Author: dschanen 4 May 2007
################################################################################

DIRS=`ls -d ens_*`

for dir in $DIRS; do
	# The position in the array is the iteration of the directory
	let index=`echo $dir | tr -dc '[:digit:]'`

	# Prepare the string by removing the $$ and blank spaces
	LIST[$index]=`grep -i '\$' $dir/tune.log | tr -d '$''[:blank:]'`
done

#echo "${LIST[@]}"

# Put the results in a file in the seqence they were run
printf "Cost function \t Iteration \n" > results.txt
for (( i=1; i <= "${#LIST[@]}"; i++ )); do
	printf "%11.4f \t %d \n" "${LIST[$i]}" "$i" >> results.txt
done

# Sort from lowest to highest cost function J
sort -n results.txt

# Or top ten results
#sort -n results.txt | tail
