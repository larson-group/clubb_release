import netCDF4
import numpy as np
import pylab as pl

silhs_files = [ 'rcm/rico_lh_lh_sfc.nc', 'rcm_in_cloud/rico_lh_lh_sfc.nc' ]
silhs_labels = [ 'rcm', 'rcm_in_cloud' ]

clubb_nc = netCDF4.Dataset('rcm/rico_lh_zt.nc')
altitude = clubb_nc.variables['altitude']

silhs_ncs = list()
for silhs_file in silhs_files:
    silhs_ncs.append(netCDF4.Dataset(silhs_file))

silhs_vars = list()
for silhs_nc in silhs_ncs:
    silhs_vars.append(silhs_nc.variables['k_lh_start'])

time1 = 0
time2 = 4320

silhs_vars_plt = list()
for silhs_var in silhs_vars:
    silhs_vars_plt.append(np.zeros(time2-time1))

for t in range(time1,time2):
    for u in range(0,len(silhs_vars_plt)):
        silhs_vars_plt[u][t-time1] = altitude[int(silhs_vars[u][t,0,0,0])-1]

for u in range(0,len(silhs_vars_plt)):
    pl.plot(range(time1,time2), silhs_vars_plt[u][:], label=silhs_labels[u])

pl.legend()
pl.show()
