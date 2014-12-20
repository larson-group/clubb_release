import netCDF4
import numpy as np
import pylab as pl

clubb_nc     = netCDF4.Dataset('old/rico_lh_zt.nc')

silhs_files = ['old/rico_lh_lh_zt.nc', 'test/rico_lh_lh_zt.nc']
silhs_labels = [ 'old', 'new']
silhs_sfc_nc = netCDF4.Dataset('old/rico_lh_lh_sfc.nc')

silhs_ncs = list()
for silhs_file in silhs_files:
    silhs_ncs.append(netCDF4.Dataset(silhs_file))

clubb_var = clubb_nc.variables['rrm_cond']

l_time_shift = False

silhs_vars = list()
for silhs_nc in silhs_ncs:
    silhs_vars.append(silhs_nc.variables['lh_rrm_evap'])
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
    pl.plot(range(time1,time2), silhs_vars_plt[u][:], label=silhs_labels[u])

pl.show()
pl.legend()