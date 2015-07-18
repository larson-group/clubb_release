#!/usr/bin/env python
import sys
import os

model_in_file = sys.argv[1]
model_out_file = sys.argv[2]

state = 0

with open(model_in_file, 'r') as f_in:
    with open(model_out_file, 'w') as f_out:
        for line in f_in:
            if state == 0 and line.rstrip() == '&microphysics_setting':
                state = 1
            elif state == 1 and line.rstrip() == '/':
                f_out.write('! Addition for using SILHS\n')
                f_out.write('lh_microphys_type = "non-interactive",\n')
                f_out.write('l_local_kk = .false.,\n')
                f_out.write('l_lh_importance_sampling = .true.,\n')
                f_out.write('l_lh_clustered_sampling = .true.,\n')
                f_out.write('l_rcm_in_cloud_k_lh_start = .true.,\n')
                f_out.write('importance_prob_thresh = 1.0e-8,\n')
                f_out.write('l_lh_limit_weights = .true.,\n')
                f_out.write('cluster_allocation_strategy = 3,\n')
                f_out.write('lh_num_samples = 8,\n')
                f_out.write('lh_seed = 1,\n')
                f_out.write('lh_sequence_length = 1\n\n')
                state = 2
            f_out.write(line)
