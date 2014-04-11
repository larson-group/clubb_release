#!/bin/bash

echo "quit" | (sudo -u matlabuser matlab -nodisplay -nodesktop -r twp_ice_profiles_creator"( 0 )")

echo "quit" | (sudo -u matlabuser matlab -nodisplay -nodesktop -r twp_ice_timeseries_creator"( 0 )")
