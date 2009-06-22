#!/bin/bash
###############################################################################
# This script edits the matlab file used to generate the output required
# for the Cloud Feedback cases. The script uses vim editing to switch between
# the different cases.
###############################################################################

# Function to run the output generator
gen_output()
{
	sudo -u matlabuser matlab -nodisplay -nodesktop -r cloud_feedback_output_creator
}

# Assume the file is configured to run the S6 case and run it.
gen_output

# Configure the file to run the S6_p2k case
vim -E -s create_output.bash <<-EOF
	:%s/cloud_feedback_s6/cloud_feedback_s6_p2k/	
	:update
	:quit
EOF

gen_output

# Configure the file to run the S11 case
vim -E -s create_output.bash <<-EOF
	:%s/cloud_feedback_s6_p2k/cloud_feedback_s11/	
	:update
	:quit
EOF

gen_output

# Configure the file to run the S11_p2k case
vim -E -s create_output.bash <<-EOF
	:%s/cloud_feedback_s11/cloud_feedback_s11_p2k/	
	:update
	:quit
EOF

gen_output

# Configure the file to run the S12 case
vim -E -s create_output.bash <<-EOF
	:%s/cloud_feedback_s11_p2k/cloud_feedback_s12/	
	:update
	:quit
EOF

gen_output

# Configure the file to run the S12_p2k case
vim -E -s create_output.bash <<-EOF
	:%s/cloud_feedback_s12/cloud_feedback_s12_p2k/	
	:update
	:quit
EOF

gen_output

# Revert the file to S6
vim -E -s create_output.bash <<-EOF
	:%s/cloud_feedback_s12_p2k/cloud_feedback_s6/	
	:update
	:quit
EOF
