#!/bin/bash
###############################################################################
# This script edits the matlab files used to verify output required for the
# Cloud Feedback cases. It uses vim editing to switch between the different
# cases.
###############################################################################

# Function to create pdf files
gen_output()
{
	echo "quit" | (sudo -u matlabuser matlab -nodisplay -nodesktop -r cloud_feedback_profiles_plot)
	echo "quit" | (sudo -u matlabuser matlab -nodisplay -nodesktop -r cloud_feedback_timeseries_plot)
}


# Configure the file to run the S6 case
vim -E -s cloud_feedback_profiles_plot <<-EOF
	:%s/curr_case\s*=\s*\'[a-zA-Z0-9]*\'/curr_case = \'s6\'/
	:update
	:quit
EOF

vim -E -s cloud_feedback_timeseries_plot <<-EOF
	:%s/curr_case\s*=\s*\'[a-zA-Z0-9]*\'/curr_case = \'s6\'/
	:update
	:quit
EOF

gen_output

# Configure the file to run the S6_p2k case
vim -E -s cloud_feedback_profiles_plot <<-EOF
	:%s/curr_case\s*=\s*\'[a-zA-Z0-9]*\'/curr_case = \'s6_p2k\'/
	:update
	:quit
EOF

vim -E -s cloud_feedback_timeseries_plot <<-EOF
	:%s/curr_case\s*=\s*\'[a-zA-Z0-9]*\'/curr_case = \'s6_p2k\'/
	:update
	:quit
EOF

gen_output

# Configure the file to run the S11 case
vim -E -s cloud_feedback_profiles_plot <<-EOF
	:%s/curr_case\s*=\s*\'[a-zA-Z0-9]*\'/curr_case = \'s11\'/
	:update
	:quit
EOF

vim -E -s cloud_feedback_timeseries_plot <<-EOF
	:%s/curr_case\s*=\s*\'[a-zA-Z0-9]*\'/curr_case = \'s11\'/
	:update
	:quit
EOF

gen_output

# Configure the file to run the S11_p2k case
vim -E -s cloud_feedback_profiles_plot <<-EOF
	:%s/curr_case\s*=\s*\'[a-zA-Z0-9]*\'/curr_case = \'s11_p2k\'/
	:update
	:quit
EOF

vim -E -s cloud_feedback_timeseries_plot <<-EOF
	:%s/curr_case\s*=\s*\'[a-zA-Z0-9]*\'/curr_case = \'s11_p2k\'/
	:update
	:quit
EOF

gen_output

# Configure the file to run the S12 case
vim -E -s cloud_feedback_profiles_plot <<-EOF
	:%s/curr_case\s*=\s*\'[a-zA-Z0-9]*\'/curr_case = \'s12\'/
	:update
	:quit
EOF

vim -E -s cloud_feedback_timeseries_plot <<-EOF
	:%s/curr_case\s*=\s*\'[a-zA-Z0-9]*\'/curr_case = \'s12\'/
	:update
	:quit
EOF

gen_output

# Configure the file to run the S12_p2k case
vim -E -s cloud_feedback_profiles_plot <<-EOF
	:%s/curr_case\s*=\s*\'[a-zA-Z0-9]*\'/curr_case = \'s12_p2k\'/
	:update
	:quit
EOF

vim -E -s cloud_feedback_timeseries_plot <<-EOF
	:%s/curr_case\s*=\s*\'[a-zA-Z0-9]*\'/curr_case = \'s12_p2k\'/
	:update
	:quit
EOF
