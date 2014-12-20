import netCDF4
import numpy as np
import pylab as pl

sim_points_all = [50, 100, 200, 300, 400, 500, 600, 800, 1000, 1200, 1500, 2000]

clubb_var_str = 'rrm_auto'
silhs_var_str  = 'lh_rrm_auto'

clubb_var = netCDF4.Dataset('silhs_'+str(sim_points_all[0])+'/rico_lh_zt.nc').variables[clubb_var_str]

k_lh_start = netCDF4.Dataset('silhs_k_lh_start/rico_lh_lh_sfc.nc').variables['k_lh_start']

rms_all = list()

for num_samples in sim_points_all:
    rms_val = 0.0
    silhs_var = netCDF4.Dataset('silhs_'+str(num_samples)+'/rico_lh_lh_zt.nc').variables[silhs_var_str]
    n_timesteps = silhs_var.shape[1]
    for t in range(0,n_timesteps):
        k = int(round(k_lh_start[t,0,0,0]))-1
        rms_val = rms_val + (clubb_var[t,k,0,0]-silhs_var[t,k,0,0])**2
    rms_val = rms_val / n_timesteps
    rms_val = np.sqrt(rms_val)
    rms_all.append(rms_val)

pl.plot(sim_points_all, rms_all)
pl.show()