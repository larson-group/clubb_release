#!/bin/bash



for ENSEMBLE in output/*
do
	ensNum=${ENSEMBLE##output/ensemble}

	folderName=ensemble$ensNum

	#Convert some output files
	echo "quit" | (sudo -u matlabuser matlab -nodisplay -nodesktop -r twp_ice_profiles_creator"( $ensNum )")
	echo "quit" | (sudo -u matlabuser matlab -nodisplay -nodesktop -r twp_ice_timeseries_creator"( $ensNum )")

	mv ../../postprocessing/output_scripts/twp_ice/*.nc submission
done
