#!/usr/bin/env python
from __future__ import print_function
import netCDF4
import numpy as np
import pylab as pl
import sys
import os

#######################################################################
# TODO: NOT FINISHED!!
#
# This script will take several directories output from the scipt
# silhs_varying_sp_output.sh, given as command line parameters to
# the script, and compute the average RMSE for each number of sample
# points for each of the directories, and plot this information. This
# version of the script creates one plot with four subpanels.
#
# This script is to be used only with non-interactive SILHS runs.
#
# Usage:
#  ./silhs_rms_4pan_mult_sim.py [options] path1 [path2 [...]]
#######################################################################

case_name = 'rico_lh'
time1 = 0
time2 = 4320

clubb_var_strs  = [ 'rrm_auto',    'rrm_accr',    'rrm_cond',    'rrm_mc_nonadj' ]
silhs_var_strs  = [ 'lh_rrm_auto', 'lh_rrm_accr', 'lh_rrm_evap', 'lh_rrm_mc_nonadj' ]

plot_sup_title = ''

output_file = ''

#-------------------------------------------------------------------------

from pylab import rcParams
rcParams['figure.figsize'] = 11, 8

silhs_dirs = []

# Read command line arguments
i = 1
while i < len(sys.argv):
    if sys.argv[i] == '--plot_sup_title':
        i = i + 1
        plot_sup_title = sys.argv[i]
    elif sys.argv[i] == '--time1':
        i = i + 1
        time1 = int(sys.argv[i])
    elif sys.argv[i] == '--time2':
        i = i + 1
        time2 = int(sys.argv[i])
    elif sys.argv[i] == '--case_name':
        i = i + 1
        case_name = sys.argv[i]
    elif sys.argv[i] == '--output_file':
        i = i + 1
        output_file = sys.argv[i]
    else:
        silhs_dirs.append(sys.argv[i])

    i = i + 1

if len(silhs_dirs) == 0:
    print("Usage: silhs_rms_4pan_mult_sim.py [options] "
          "path1 [path2 [...]]", file=sys.stderr)
    sys.exit(1)

sim_points_all = list()
for entry in os.listdir(silhs_dirs[0]):
    if (entry[:6] == 'silhs_'):
        sim_points_all.append(int(entry[6:]))
sim_points_all.sort()

dir_names = list()
for d in silhs_dirs:
    dir_names.append(os.path.basename(d))

lines = list()

for plot_num in range(4):

    clubb_var = netCDF4.Dataset(silhs_dirs[0]+'/silhs_'+str(sim_points_all[0])+ \
      '/'+case_name+'_zt.nc').variables[clubb_var_strs[plot_num]]

    rms_all = list()
    for i in range(len(silhs_dirs)):
        rms_all.append(np.empty(len(sim_points_all)))

    n_timesteps = time2-time1

    k_lh_start = netCDF4.Dataset(silhs_dirs[0]+'/silhs_'+str(sim_points_all[0])+ \
      '/'+case_name+'_lh_sfc.nc').variables['k_lh_start']

    # Copy to memory for faster access
    k_lh_start = k_lh_start[:,0,0,0]
    clubb_var = clubb_var[:,:,0,0]

    for n_i in range(0,len(sim_points_all)):
        num_samples = sim_points_all[n_i]

        for d_i in range(0,len(silhs_dirs)):
            silhs_dir = silhs_dirs[d_i]
            rms_val = 0.0
            silhs_var = netCDF4.Dataset(silhs_dir+'/silhs_'+str(num_samples)+ \
                '/'+case_name+'_lh_zt.nc').variables[silhs_var_strs[plot_num]]
            # Copy to memory for better performance
            silhs_var = silhs_var[:,:,0,0]

            for t in range(time1,time2):
                k = int(round(k_lh_start[t]))-1
                rms_val = rms_val + (clubb_var[t,k]-silhs_var[t,k])**2

            rms_val = rms_val / n_timesteps
            rms_val = np.sqrt(rms_val)
            rms_all[d_i][n_i] = rms_val
    pl.subplot(2, 2, plot_num+1)

    for d_i in range(len(silhs_dirs)):
        format_str = ''
        if os.path.basename(silhs_dirs[d_i]) == 'cloud_weighted':
            format_str = 'b-'
        elif os.path.basename(silhs_dirs[d_i]) == 'prescribed':
            format_str = 'g--'

        line, = pl.plot(sim_points_all, rms_all[d_i], format_str, \
                label=os.path.basename(silhs_dirs[d_i]))
        if plot_num == 0:
            lines.append(line)

    if plot_num >= 2:
        pl.xlabel('Number of Sample Points')
    if plot_num == 0 or plot_num == 2:
        pl.ylabel('RMSE of SILHS estimate')
    pl.xscale('log')
    pl.yscale('log')
    if clubb_var_strs[plot_num] == 'rrm_auto':
        eq = '$\\left(\\frac{\partial r_r}{\\partial t}\\right)_\\mathrm{auto}$'
        pl.title('Rain autoconversion tendency ' + eq)
    elif clubb_var_strs[plot_num] == 'rrm_accr':
        eq = '$\\left(\\frac{\\partial r_r}{\\partial t}\\right)_\\mathrm{accr}$'
        pl.title('Rain accretion tendency ' + eq)
    elif clubb_var_strs[plot_num] == 'rrm_cond':
        eq = '$\\left(\\frac{\\partial r_r}{\\partial t}\\right)_\\mathrm{evap}$'
        pl.title('Rain evaporation tendency ' + eq)
    elif clubb_var_strs[plot_num] == 'rrm_mc_nonadj':
        eq = '$\\left(\\frac{\\partial r_r}{\\partial t}\\right)$'
        pl.title('Total rain tendency ' + eq)

pl.figlegend( lines, dir_names, 'lower center', ncol=2, fontsize=7 )

if output_file == '':
    pl.show()
else:
    # Output to disk
    pl.savefig(output_file)
