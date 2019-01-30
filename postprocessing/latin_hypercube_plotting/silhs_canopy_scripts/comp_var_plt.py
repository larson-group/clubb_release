import netCDF4
import numpy as np
import pylab as pl
import random

clubb_nc     = netCDF4.Dataset('rico_silhs_zt.nc')
silhs_nc     = netCDF4.Dataset('rico_silhs_lh_zt.nc')
silhs_sfc_nc = netCDF4.Dataset('rico_silhs_lh_sfc.nc')

clubb_var = clubb_nc.variables['mixt_frac']
silhs_var = silhs_nc.variables['lh_mixt_frac']

k_lh_start = silhs_sfc_nc.variables['k_lh_start']

time1 = 4000
time2 = 4320

clubb_var_plt = np.empty(time2-time1)
silhs_var_plt = np.empty(time2-time1)
comp_var_plt  = np.empty(time2-time1)

for t in range(time1,time2):
    k = int(round(k_lh_start[t,0,0,0])) - 1
    clubb_var_plt[t-time1] = clubb_var[t,k,0,0]
    silhs_var_plt[t-time1] = silhs_var[t,k,0,0]
    # Code added for new line
    num_samples = 100
    n_lt = 0
    mf_element = clubb_var_plt[t-time1]
    for i in range(0,num_samples):
        if ( random.random() < mf_element ):
            n_lt = n_lt + 1
    comp_var_plt[t-time1] = float(n_lt) / num_samples

pl.plot(range(time1,time2), clubb_var_plt[:])
pl.plot(range(time1,time2), silhs_var_plt[:])
pl.plot(range(time1,time2), comp_var_plt[:])
pl.show()

clubb_nc.close()
silhs_nc.close()
silhs_sfc_nc.close()