#!/bin/bash

echo "quit" | (sudo -u matlabuser matlab -nodisplay -nodesktop -r cloud_feedback_output_plot)
