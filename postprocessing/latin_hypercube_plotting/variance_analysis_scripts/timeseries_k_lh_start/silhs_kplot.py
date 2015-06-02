#!/usr/bin/env python
from __future__ import print_function
import netCDF4
import numpy as np
import pylab as pl
import sys
import os

#---------------------------------------------------------------------
# Plots a SILHS variable and its analytic equivalent as a time series
# at k_lh_start, for various SILHS directories.
#
# Usage:
#   silhs_kplot.py dir1 [dir2 [...]]
#---------------------------------------------------------------------

clubb_var = 'rrm_mc'
silhs_var = 'lh_rrm_mc'

if len(sys.argv) <= 1:
    print("Usage: ./silhs_kplot.py dir1 [dir2 [...]]", file=sys.stderr)
    sys.exit(1)

silhs_dirs = sys.argv[1:]

clubb_nc     = netCDF4.Dataset(silhs_dirs[0] + '/rico_lh_zt.nc')

silhs_sfc_nc = netCDF4.Dataset(silhs_dirs[0] + '/rico_lh_lh_sfc.nc')

silhs_ncs = list()
for silhs_dir in silhs_dirs:
    silhs_ncs.append(netCDF4.Dataset(silhs_dir + '/rico_lh_lh_zt.nc'))

clubb_var = clubb_nc.variables[clubb_var]
silhs_vars = list()
for silhs_nc in silhs_ncs:
    silhs_vars.append(silhs_nc.variables[silhs_var])

# In CLUBB (as of the time of writing this script), hydrometeor means such
# as "rrm" are output one timestep earlier than the corresponding SILHS
# estimate (in this case, "lh_rrm"). Thus, when (and only when) comparing a
# hydrometeor mean and the corresponding SILHS estimate, set this flag to
# True.
l_time_shift = False

k_lh_start = silhs_sfc_nc.variables['k_lh_start']

time1 = 3000
time2 = 4320

clubb_var_plt = np.empty(time2-time1)
silhs_vars_plt = list()
for silhs_var in silhs_vars:
    silhs_vars_plt.append(np.zeros(time2-time1))

for t in range(time1,time2):
    k = int(round(k_lh_start[t,0,0,0])) - 1
    if l_time_shift:
        clubb_var_plt[t-time1] = clubb_var[t-1,k,0,0]
    else:
        clubb_var_plt[t-time1] = clubb_var[t,k,0,0]

    for u in range(0,len(silhs_vars_plt)):
        silhs_vars_plt[u][t-time1] = silhs_vars[u][t,k,0,0]

pl.plot(range(time1,time2), clubb_var_plt[:], label='analytic')
for u in range(0,len(silhs_vars_plt)):
    pl.plot(range(time1,time2), silhs_vars_plt[u][:], label=os.path.basename(silhs_dirs[u]))

pl.legend()
pl.show()
