#!/bin/bash

# Script to sort results of a series of tuning run by the value of the
# optimal error of the solution found.

for exp in bomex_404 combined_404 dycoms2_rf01_404 \
           bomex_409 combined_409 dycoms2_rf01_409
do

echo Processing results for ${exp}
cd ${exp}

file=${exp}_summary.dat
tmp=tmp.dat
rm -f $tmp

n=0
nmax=400
for (( x = 1; x <= 6; x++))
do

# Check that there is only one each of the following files: *.log, error*.in

  if [ `ls -1 ${exp}_${x}/error*.in | wc -l` -ne 1 ]; then
    echo WARNING: incorrect number of error files for ${exp}_${x}
    exit 1
  fi

  if [ `ls -1 ${exp}_${x}/*.log | wc -l` -ne 1 ]; then
    echo WARNING: incorrect number of log files for ${exp}_${x}
    exit 1
  fi

# Extract the optimized constants if they resulted in a valid
# standalone run. Occasionally, some do not.

  dmget ${exp}_${x}/*.log
  if [ `grep wasn ${exp}_${x}/*.log | wc -l` -eq 0 ]; then
    ( printf "%11.4f," `grep 'optimal error' ${exp}_${x}/*.log | cut -c24-35`
      printf "%11d," ${x}
      printf "%11.4f," `grep 'C1 ' ${exp}_${x}/error*.in | cut -c18-31`
      printf "%11.4f," `grep 'C1b ' ${exp}_${x}/error*.in | cut -c18-31`
      printf "%11.4f," `grep 'C1c ' ${exp}_${x}/error*.in | cut -c18-31`
      printf "%11.4f," `grep 'C2 ' ${exp}_${x}/error*.in | cut -c18-31`
      printf "%11.4f," `grep 'C2b ' ${exp}_${x}/error*.in | cut -c18-31`
      printf "%11.4f," `grep 'C2c ' ${exp}_${x}/error*.in | cut -c18-31`
      printf "%11.4f," `grep 'C5 ' ${exp}_${x}/error*.in | cut -c18-31`
      printf "%11.4f," `grep 'C6rt ' ${exp}_${x}/error*.in | cut -c18-31`
      printf "%11.4f," `grep 'C6rtb ' ${exp}_${x}/error*.in | cut -c18-31`
      printf "%11.4f," `grep 'C6rtc ' ${exp}_${x}/error*.in | cut -c18-31`
      printf "%11.4f," `grep 'C6thl ' ${exp}_${x}/error*.in | cut -c18-31`
      printf "%11.4f," `grep 'C6thlb ' ${exp}_${x}/error*.in | cut -c18-31`
      printf "%11.4f," `grep 'C6thlc ' ${exp}_${x}/error*.in | cut -c18-31`
      printf "%11.4f," `grep 'C7 ' ${exp}_${x}/error*.in | cut -c18-31`
      printf "%11.4f," `grep 'C7b ' ${exp}_${x}/error*.in | cut -c18-31`
      printf "%11.4f," `grep 'C7c ' ${exp}_${x}/error*.in | cut -c18-31`
      printf "%11.4f," `grep 'C8 ' ${exp}_${x}/error*.in | cut -c18-31`
      printf "%11.4f," `grep 'C11 ' ${exp}_${x}/error*.in | cut -c18-31`
      printf "%11.4f," `grep 'C11b ' ${exp}_${x}/error*.in | cut -c18-31`
      printf "%11.4f," `grep 'C11c ' ${exp}_${x}/error*.in | cut -c18-31`
      printf "%11.4f," `grep 'beta ' ${exp}_${x}/error*.in | cut -c18-31`
      printf "%11.4f," `grep 'gamma_coef ' ${exp}_${x}/error*.in | cut -c18-31`
      printf "%11.4f," `grep 'gamma_coefb ' ${exp}_${x}/error*.in | cut -c18-31`
      printf "%11.4f," `grep 'gamma_coefc ' ${exp}_${x}/error*.in | cut -c18-31`
      printf "%11.4f," `grep 'c_K ' ${exp}_${x}/error*.in | cut -c18-31`
      printf "%11.4e," `grep 'mu ' ${exp}_${x}/error*.in | cut -c18-31`
      printf "\n" ) >> $tmp
    ((n++))
  fi

  if [ $n -ge $nmax ]; then
    echo Extracted $n valid simulations for ${exp}
    break
  fi

done

echo Sorting results for ${exp}
sort -b -g $tmp > $file
rm -f $tmp

cd ..

done
