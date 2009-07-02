#!/bin/bash

clubbDir="~/clubb"

for ENSEMBLE in ensembles/*
do
	echo "Running $ENSEMBLE..."
	ensNum=${ENSEMBLE##ensembles/ensmem}

	#Create a directory to store the output
	folderName=ensemble$ensNum
	mkdir "output/$folderName"

	#Create the sounding and forcing data
	dataName="twpice_p$ensNum"_layer.dat
	cp "$ENSEMBLE/$dataName" .
	sudo -u matlabuser matlab -nodisplay -nodesktop -r $clubbDir/utilities/case_setup/twp_ice/convert_arm_input_data"('twp_ice', '$dataName', 0, 1004.4 )"
	mv -f "twp_ice_sounding.in" "output/$folderName"
	mv -f "twp_ice_forcings.in" "output/$folderName"
	mv -f "twp_ice_surface.in" "output/$folderName"
	rm $dataName
done
