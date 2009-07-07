#!/bin/bash

clubbDir="/home/nielsenb/clubb"

for ENSEMBLE in output/*
do
	ensNum=${ENSEMBLE##output/ensemble}

	folderName=ensemble$ensNum

	returnDir=`pwd`

	cd $clubbDir/postprocessing/output_scripts/twp_ice

	#Convert some output files
	echo "quit" | (sudo -u matlabuser matlab -nodisplay -nodesktop -r twp_ice_profiles_creator"( $ensNum )")
	echo "quit" | (sudo -u matlabuser matlab -nodisplay -nodesktop -r twp_ice_timeseries_creator"( $ensNum )")

	cd $returnDir

	mv "$clubbDir"/postprocessing/output_scripts/twp_ice/*.nc output/"$folderName"
	cp output/"$folderName"/*.nc submission
done
