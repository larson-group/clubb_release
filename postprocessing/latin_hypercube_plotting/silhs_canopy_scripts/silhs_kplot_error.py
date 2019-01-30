import netCDF4
import numpy as np
import pylab as pl

clubb_nc     = netCDF4.Dataset('cldwt/rico_silhs_zt.nc')

silhs_files = [ 'cldwt/rico_silhs_lh_zt.nc' ]
silhs_labels = [ 'cldwt' ]
silhs_sfc_nc = netCDF4.Dataset('cldwt/rico_silhs_lh_sfc.nc')

silhs_ncs = list()
for silhs_file in silhs_files:
    silhs_ncs.append(netCDF4.Dataset(silhs_file))

clubb_var = clubb_nc.variables['rrm']
variance_var = clubb_nc.variables['rrp2_zt']
l_time_shift = True

silhs_vars = list()
for silhs_nc in silhs_ncs:
    silhs_vars.append(silhs_nc.variables['lh_rrm'])

k_lh_start = silhs_sfc_nc.variables['k_lh_start']

time1 = 3000
time2 = 4320

clubb_var_plt = np.empty(time2-time1)
variance_var_plt = np.empty(time2-time1)

for t in range(time1,time2):
    k = int(round(k_lh_start[t,0,0,0])) - 1
    clubb_var_plt[t-time1] = np.abs(clubb_var[t-1,k,0,0] - silhs_vars[0][t,k,0,0])
    variance_var_plt[t-time1] = np.sqrt(variance_var[t,k,0,0]) / 10.0

pl.plot(range(time1,time2), clubb_var_plt[:], label='error')
pl.plot(range(time1,time2), variance_var_plt[:], label='standard deviation')

pl.show()
pl.legend()

clubb_nc.close()
for silhs_nc in silhs_ncs:
    silhs_nc.close()
silhs_sfc_nc.close()