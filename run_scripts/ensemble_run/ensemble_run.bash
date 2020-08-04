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
	cp "$ENSEMBLE/$caseName"_surface.in ../../input/case_setups

	#Clean up any old output files before running CLUBB
	rm -f ../../output/"$caseName"_zt.ctl
	rm -f ../../output/"$caseName"_zt.dat
	rm -f ../../output/"$caseName"_zm.ctl
	rm -f ../../output/"$caseName"_zm.dat
	rm -f ../../output/"$caseName"_sfc.ctl
	rm -f ../../output/"$caseName"_sfc.dat

	#Run CLUBB
	../../run_scripts/run_scm.bash $caseName --grads

	#Harvest the output files
	cp ../../output/twp_ice_zt.ctl output/"$folderName"
	cp ../../output/twp_ice_zt.dat output/"$folderName"
	cp ../../output/twp_ice_zm.ctl output/"$folderName"
	cp ../../output/twp_ice_zm.dat output/"$folderName"
	cp ../../output/twp_ice_sfc.ctl output/"$folderName"
	cp ../../output/twp_ice_sfc.dat output/"$folderName"
done
