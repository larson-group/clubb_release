import netCDF4
import numpy as np
import pylab as pl

clubb_nc       = netCDF4.Dataset('rico_lh_zt.nc')
silhs_nc       = netCDF4.Dataset('rico_lh_lh_zt.nc')
silhs_sfc_nc   = netCDF4.Dataset('rico_lh_lh_sfc.nc')
silhs_2D_nl_nc = netCDF4.Dataset('rico_lh_nl_lh_sample_points_2D.nc')

clubb_var = clubb_nc.variables['rrm_auto']
silhs_var = silhs_nc.variables['lh_rrm_auto']

k_lh_start = silhs_sfc_nc.variables['k_lh_start']

time1 = 3000
time2 = 4320

clubb_var_plt = np.empty(time2-time1)
silhs_var_plt = np.empty(time2-time1)

for t in range(time1,time2):
    k = int(round(k_lh_start[t,0,0,0])) - 1
    clubb_var_plt[t-time1] = clubb_var[t,k,0,0]
    silhs_var_plt[t-time1] = silhs_var[t,k,0,0]

pl.plot(range(time1,time2), clubb_var_plt[:], label="analytic")
pl.plot(range(time1,time2), silhs_var_plt[:], label="SILHS")
pl.show()
pl.legend()

clubb_nc.close()
silhs_nc.close()
silhs_sfc_nc.close()