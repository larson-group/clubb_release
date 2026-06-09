#!/bin/bash

caseName="twp_ice"

for ENSEMBLE in ensembles/*
do
	ensNum=${ENSEMBLE##ensembles/ensemble}
	echo "Running $ensNum..."

	#Create a directory to store the output
	folderName=ensemble$ensNum
	mkdir "output/$folderName"

	#Setup CLUBB using the files we created
	cp "$ENSEMBLE/$caseName"_sounding.in ../../input/case_setups
	cp "$ENSEMBLE/$caseName"_forcings.in ../../input/case_setups
	cp "$ENSEMBLE/$caseName"_sfc.in ../../input/case_setups

	#Clean up any old output files before running CLUBB
	rm -f ../../output/"$caseName"_zt.nc
	rm -f ../../output/"$caseName"_zm.nc
	rm -f ../../output/"$caseName"_sfc.nc

	#Run CLUBB
	../../run_scripts/run_scm.py $caseName

	#Harvest the output files
	cp ../../output/twp_ice_zt.nc output/"$folderName"
	cp ../../output/twp_ice_zm.nc output/"$folderName"
	cp ../../output/twp_ice_sfc.nc output/"$folderName"
done
