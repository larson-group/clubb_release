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
                f_out.write('eight_cluster_presc_probs%cloud_precip_comp1      = 0.125,\n')
                f_out.write('eight_cluster_presc_probs%cloud_precip_comp2      = 0.125,\n')
                f_out.write('eight_cluster_presc_probs%nocloud_precip_comp1    = 0.125,\n')
                f_out.write('eight_cluster_presc_probs%nocloud_precip_comp2    = 0.125,\n')
                f_out.write('eight_cluster_presc_probs%cloud_noprecip_comp1    = 0.125,\n')
                f_out.write('eight_cluster_presc_probs%cloud_noprecip_comp2    = 0.125,\n')
                f_out.write('eight_cluster_presc_probs%nocloud_noprecip_comp1  = 0.125,\n')
                f_out.write('eight_cluster_presc_probs%nocloud_noprecip_comp2  = 0.125\n')

                state = 2
            f_out.write(line)
