#!/bin/bash

for ENSEMBLE in input/*
do
	echo "Running $ENSEMBLE..."
	ensNum=${ENSEMBLE##input/ensmem}

	#Create a directory to store the output
	folderName=ensemble$ensNum
	mkdir "ensembles/$folderName"

	#Create the sounding and forcing data
	dataName="twpice_p$ensNum"_layer.dat
	cp "$ENSEMBLE/$dataName" .
	sudo -u matlabuser matlab -nodisplay -nodesktop -r convert_arm_input_data"('twp_ice', '$dataName', 0, 1004.4 )"
	mv -f "twp_ice_sounding.in" "ensembles/$folderName"
	mv -f "twp_ice_forcings.in" "ensembles/$folderName"
	mv -f "twp_ice_sfc.in" "ensembles/$folderName"
	rm $dataName
done
