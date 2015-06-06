#!/usr/bin/env python
from __future__ import print_function
import netCDF4
import numpy as np
import pylab as pl
import sys
import os

#######################################################################
# This script will take several directories output from the scipt
# silhs_varying_sp_output.sh, given as command line parameters to
# the script, and compute the average RMSE for each number of sample
# points for each of the directories, and plot this information.
#
# This script is to be used only with non-interactive SILHS runs.
#
# Usage:
#  ./silhs_rms_plot_mult_sim.py [options] path1 [path2 [...]]
#######################################################################

case_name = 'rico_lh'
time1 = 0
time2 = 4320

# Average over all height levels. If false, RMSE will be computed only at
# k_lh_start
l_all_height_avg = False

clubb_var_str = 'rrm_mc'
silhs_var_str  = 'lh_rrm_mc'

plot_title = clubb_var_str

output_file = ''

#-------------------------------------------------------------------------

silhs_dirs = []

# Read command line arguments
i = 1
while i < len(sys.argv):
    if sys.argv[i] == '--plot_title':
        i = i + 1
        plot_title = sys.argv[i]
    elif sys.argv[i] == '--time1':
        i = i + 1
        time1 = int(sys.argv[i])
    elif sys.argv[i] == '--time2':
        i = i + 1
        time2 = int(sys.argv[i])
    elif sys.argv[i] == '--case_name':
        i = i + 1
        case_name = sys.argv[i]
    elif sys.argv[i] == '--clubb_var_str':
        i = i + 1
        clubb_var_str = sys.argv[i]
    elif sys.argv[i] == '--silhs_var_str':
        i = i + 1
        silhs_var_str = sys.argv[i]
    elif sys.argv[i] == '--output_file':
        i = i + 1
        output_file = sys.argv[i]
    else:
        silhs_dirs.append(sys.argv[i])

    i = i + 1

if len(silhs_dirs) == 0:
    print("Usage: ./silhs_rms_plot_mult_sim.py [options] "
          "path1 [path2 [...]]", file=sys.stderr)
    sys.exit(1)

sim_points_all = list()
for entry in os.listdir(silhs_dirs[0]):
    if (entry[:6] == 'silhs_'):
        sim_points_all.append(int(entry[6:]))
sim_points_all.sort()

clubb_var = netCDF4.Dataset(silhs_dirs[0]+'/silhs_'+str(sim_points_all[0])+ \
    '/'+case_name+'_zt.nc').variables[clubb_var_str]

rms_all = list()
for i in range(len(silhs_dirs)):
    rms_all.append(np.empty(len(sim_points_all)))

n_timesteps = time2-time1

if not l_all_height_avg:
    k_lh_start = netCDF4.Dataset(silhs_dirs[0]+'/silhs_'+str(sim_points_all[0])+ \
      '/'+case_name+'_lh_sfc.nc').variables['k_lh_start']
    k_lh_start = k_lh_start[:,0,0,0]
else:
    n_heights = clubb_var.shape[1]

clubb_var = clubb_var[:,:,0,0]

for n_i in range(0,len(sim_points_all)):
    num_samples = sim_points_all[n_i]
    print("Trying " + str(num_samples) + " sample points.")
    for d_i in range(0,len(silhs_dirs)):
        silhs_dir = silhs_dirs[d_i]
        rms_val = 0.0
        silhs_var = netCDF4.Dataset(silhs_dir+'/silhs_'+str(num_samples)+ \
            '/'+case_name+'_lh_zt.nc').variables[silhs_var_str]
        # Copy to memory for better performance
        silhs_var = silhs_var[:,:,0,0]

        for t in range(time1,time2):
            if l_all_height_avg:
                # Exclude the model lower level
                for k in range(1,n_heights):
                    rms_val = rms_val + (clubb_var[t,k]-silhs_var[t,k])**2
            else:
                k = int(round(k_lh_start[t]))-1
                rms_val = rms_val + (clubb_var[t,k]-silhs_var[t,k])**2

        rms_val = rms_val / n_timesteps
        if l_all_height_avg:
            rms_val = rms_val / n_heights
        rms_val = np.sqrt(rms_val)
        rms_all[d_i][n_i] = rms_val

for d_i in range(len(silhs_dirs)):
    pl.plot(sim_points_all, rms_all[d_i], label=os.path.basename(silhs_dirs[d_i]))

# Create a 1/N line
## rms_1_over_n = np.empty(len(sim_points_all))
## rms_1_over_n[0] = rms_all[2][0]
## for n_i in range(1,len(sim_points_all)):
##     rms_1_over_n[n_i] = rms_1_over_n[n_i-1] * sim_points_all[n_i-1] / sim_points_all[n_i]
## pl.plot(sim_points_all, rms_1_over_n, label='1 over N')

# Create a 1/sqrt(N) line
## rms_1_over_sqrt_n = np.empty(len(sim_points_all))
## rms_1_over_sqrt_n[0] = rms_all[2][0]
## for n_i in range(1,len(sim_points_all)):
##     rms_1_over_sqrt_n[n_i] = rms_1_over_sqrt_n[n_i-1] * \
##          np.sqrt(float(sim_points_all[n_i-1]) / sim_points_all[n_i])
## pl.plot(sim_points_all, rms_1_over_sqrt_n, label='1 over sqrt N')

pl.xlabel('Number of Sample Points')
pl.ylabel('RMSE of SILHS Estimate')
pl.xscale('log')
pl.yscale('log')
pl.title(plot_title)
# This will change the axes so that the plot 'tight'ly hugs the data
#pl.axis('tight')

pl.legend()

if output_file == '':
    pl.show()
else:
    # Output to disk
    pl.savefig(output_file)
