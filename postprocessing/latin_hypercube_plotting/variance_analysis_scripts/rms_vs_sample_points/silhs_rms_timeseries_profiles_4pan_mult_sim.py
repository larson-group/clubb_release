#!/usr/bin/env python
from __future__ import print_function
import netCDF4
import numpy as np
import matplotlib.pyplot as pl
import sys
import os
import glob

#######################################################################
# This script will take several directories output from the scipt
# silhs_varying_sp_output.sh, given as command line parameters to
# the script, and compute the average RMSE for each number of sample
# points for each of the directories, and plot this information. This
# version of the script creates one plot with four subpanels.
#
# This script is to be used only with non-interactive SILHS runs.
#
# Usage:
#  ./silhs_rms_timeseries_profiles_4pan_mult_sim.py [options] path1 [path2 [...]]
#######################################################################

case_name = 'rico_lh'
time1 = 0
time2 = 864

clubb_var_strs  = [ 'rrm_auto',    'rrm_accr',    'rrm_cond',    'rrm_mc_nonadj' ]
silhs_var_strs  = [ 'lh_rrm_auto', 'lh_rrm_accr', 'lh_rrm_evap', 'lh_rrm_mc_nonadj' ]

# 0: RMS (default)
# 1: timeseries
# 2: profiles
mode = 0

plot_sup_title = ''

output_file = 'out.svg'

num_pts_for_timeseries = '16'
num_pts_for_profiles = '16'

#-------------------------------------------------------------------------

pl.rcParams['figure.figsize'] = 10, 10

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
    elif sys.argv[i] == '--rms':
        mode = 0
    elif sys.argv[i] == '--timeseries':
        mode = 1
    elif sys.argv[i] == '--profiles':
        mode = 2
    else:
        silhs_dirs.append(sys.argv[i])

    i = i + 1

if len(silhs_dirs) == 0:
    print("Usage: silhs_rms_timeseries_profiles_4pan_mult_sim.py [options] "
          "path1 [path2 [...]]", file=sys.stderr)
    sys.exit(1)

sim_points_all = set()
for entry in os.listdir(silhs_dirs[0]):
    if (entry[:6] == 'silhs_'):
        sim_points_all.add(int(entry[6:].split('_')[0]))
# Convert set to sorted list
sim_points_all = sorted(sim_points_all)

dir_names = list()
for d in silhs_dirs:
    dir_names.append(os.path.basename(d))

lines = list()

k_lh_start = netCDF4.Dataset(silhs_dirs[0]+'/silhs_'+str(sim_points_all[0])+'_1'+ \
    '/'+case_name+'_lh_sfc.nc').variables['k_lh_start']
# Copy to memory
k_lh_start = k_lh_start[:,0,0,0]

for plot_num in range(4):

    clubb_nc = netCDF4.Dataset(silhs_dirs[0]+'/silhs_'+str(sim_points_all[0])+'_1'+ \
        '/'+case_name+'_zt.nc')
    clubb_var = clubb_nc.variables[clubb_var_strs[plot_num]]
    # Copy to memory for faster access
    clubb_var = clubb_var[:,:,0,0]

    if mode == 1:
        silhs_vars = list()
        silhs_vars_plt = list()
        for silhs_dir in silhs_dirs:
            silhs_vars.append(netCDF4.Dataset(silhs_dir+'/silhs_'+num_pts_for_timeseries+'_1'+\
            '/'+case_name+'_lh_zt.nc').variables[silhs_var_strs[plot_num]][:,:,0,0])
            silhs_vars_plt.append(np.empty(time2-time1))

        for t in range(time1,time2):
            k = int(round(k_lh_start[t])) - 1
            for u in range(len(silhs_vars_plt)):
                silhs_vars_plt[u][t-time1] = silhs_vars[u][t,k]

    elif mode == 0:
        rms_all = list()
        for i in range(len(silhs_dirs)):
            rms_all.append(np.empty(len(sim_points_all)))

        n_timesteps = time2-time1

        for n_i in range(0,len(sim_points_all)):
            num_samples = sim_points_all[n_i]

            for d_i in range(0,len(silhs_dirs)):
                silhs_dir = silhs_dirs[d_i]
                rms_val = 0.0
                n_seed_dirs = 0
                for seed_dir in glob.glob(silhs_dir+'/silhs_'+str(num_samples)+'_*'):
                    n_seed_dirs += 1
                    silhs_var = netCDF4.Dataset(seed_dir+'/'+case_name+'_lh_zt.nc') \
                        .variables[silhs_var_strs[plot_num]][:,:,0,0]

                    for t in range(time1,time2):
                        k = int(round(k_lh_start[t]))-1
                        rms_val = rms_val + (clubb_var[t,k]-silhs_var[t,k])**2

                rms_val = np.sqrt(rms_val / (n_timesteps*n_seed_dirs))
                rms_all[d_i][n_i] = rms_val

    elif mode == 2:
        altitude = clubb_nc.variables['altitude'][:]
        # Exclude lower boundary
        altlow = 1
        for althigh in range(2,len(altitude)):
            if altitude[althigh] > 4000.0:
                break

    pl.subplot(2, 2, plot_num+1)

    # Plot analytic line
    if mode == 2:
        line, = pl.plot(np.average(clubb_var[time1:time2,altlow:althigh], axis=0), \
                    altitude[altlow:althigh], 'k-', linewidth=2)
        if plot_num == 0:
            lines.append(line)
            dir_names.insert(0, 'analytic')

    for d_i in range(len(silhs_dirs)):
        format_str = ''
        if os.path.basename(silhs_dirs[d_i]) == 'LH-only':
            format_str = 'm-'
        elif os.path.basename(silhs_dirs[d_i]) == '2Cat-Cld':
            format_str = 'b--'
        elif os.path.basename(silhs_dirs[d_i]) == '2Cat-CldPcp':
            format_str = 'r:'
        elif os.path.basename(silhs_dirs[d_i]) == '8Cat':
            format_str = 'c*-'
        markerSize = 3
        if mode == 1:
            line, = pl.plot(range(time1+1,time2+1), silhs_vars_plt[d_i], format_str, \
                            markersize=markerSize)
        elif mode == 0:
            line, = pl.plot(sim_points_all, rms_all[d_i], format_str, markersize=markerSize)
        elif mode == 2:
            for seed_dir in glob.glob(silhs_dirs[d_i]+'/silhs_'+num_pts_for_profiles+'_*'):
                line, = pl.plot(np.average(netCDF4.Dataset(seed_dir+'/'+case_name+'_lh_zt.nc') \
                    .variables[silhs_var_strs[plot_num]][time1:time2,altlow:althigh,0,0], axis=0), \
                    altitude[altlow:althigh], format_str, markersize=markerSize)
        if plot_num == 0:
            lines.append(line)

    if plot_num >= 2:
        if mode == 1:
            pl.xlabel('Time [min]')
        elif mode == 0:
            pl.xlabel('Number of Sample Points')
        elif mode == 2:
            pl.xlabel('Tendency $[\\rm{kg}\\ \\rm{kg}^{-1}\\ \\rm{s}^{-1}]$')
    if plot_num == 0 or plot_num == 2:
        if mode == 1:
            pl.ylabel('Tendency $[\\rm{kg}\\ \\rm{kg}^{-1}\\ \\rm{s}^{-1}]$')
        elif mode == 0:
            pl.ylabel('RMSE of SILHS estimate')
        elif mode == 2:
            pl.ylabel('Height [m]')
    if mode == 0:
        pl.xscale('log')
        pl.yscale('log')
    if clubb_var_strs[plot_num] == 'rrm_auto':
        eq = '$\\left(\\frac{\partial r_r}{\\partial t}\\right)_\\mathrm{auto}$'
        pl.title('Autoconversion ' + eq)
    elif clubb_var_strs[plot_num] == 'rrm_accr':
        eq = '$\\left(\\frac{\\partial r_r}{\\partial t}\\right)_\\mathrm{accr}$'
        pl.title('Accretion ' + eq)
    elif clubb_var_strs[plot_num] == 'rrm_cond':
        eq = '$\\left(\\frac{\\partial r_r}{\\partial t}\\right)_\\mathrm{evap}$'
        pl.title('Evaporation ' + eq)
    elif clubb_var_strs[plot_num] == 'rrm_mc_nonadj':
        eq = '$\\left(\\frac{\\partial r_r}{\\partial t}\\right)$'
        pl.title('Total rain tendency ' + eq)

lgd = pl.figlegend( lines, dir_names, 'lower center', ncol=2, fontsize=9 )

# Output to disk
pl.savefig(output_file, bbox_extra_artists=(lgd,), bbox_inches='tight')
