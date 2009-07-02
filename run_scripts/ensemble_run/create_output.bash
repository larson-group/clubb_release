#!/bin/bash

echo "quit" | (sudo -u matlabuser matlab -nodisplay -nodesktop -r twp_ice_profiles_creator"( $1 )")

echo "quit" | (sudo -u matlabuser matlab -nodisplay -nodesktop -r twp_ice_timeseries_creator"( $1 )")
