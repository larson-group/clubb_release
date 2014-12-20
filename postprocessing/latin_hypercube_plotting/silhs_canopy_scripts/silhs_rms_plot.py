import netCDF4
import numpy as np
import pylab as pl

sim_points_all = [50, 100, 200, 300, 400, 500, 600, 800, 1000, 1200, 1500, 2000]

clubb_var_str = 'rrm_cond'
silhs_var_str  = 'lh_rrm_evap'

clubb_var = netCDF4.Dataset('silhs_'+str(sim_points_all[0])+'/rico_lh_zt.nc').variables[clubb_var_str]

k_lh_start = netCDF4.Dataset('silhs_k_lh_start/rico_lh_lh_sfc.nc').variables['k_lh_start']

rms_all = list()
rms_1_over_n = list()
rms_1_over_sqrt_n = list()

for n in range(0,len(sim_points_all)):
    num_samples = sim_points_all[n]
    rms_val = 0.0
    silhs_var = netCDF4.Dataset('silhs_'+str(num_samples)+'/rico_lh_lh_zt.nc').variables[silhs_var_str]
    n_timesteps = silhs_var.shape[1]
    for t in range(0,n_timesteps):
        k = int(round(k_lh_start[t,0,0,0]))-1
        rms_val = rms_val + (clubb_var[t,k,0,0]-silhs_var[t,k,0,0])**2
    rms_val = rms_val / n_timesteps
    rms_val = np.sqrt(rms_val)
    rms_all.append(rms_val)
    if n == 0:
      rms_1_over_n.append(rms_val)
      rms_1_over_sqrt_n.append(rms_val)
    else:
      rms_1_over_n.append(rms_1_over_n[n-1] * (float(sim_points_all[n-1]) / sim_points_all[n]))
      rms_1_over_sqrt_n.append(rms_1_over_sqrt_n[n-1] * np.sqrt(float(sim_points_all[n-1]) / sim_points_all[n]))

pl.plot(sim_points_all, rms_all, label='SILHS')
pl.plot(sim_points_all, rms_1_over_n, label='1 over N')
pl.plot(sim_points_all, rms_1_over_sqrt_n, label='1 over sqrt N')

pl.xlabel('Number of Sample Points')
pl.ylabel('Root Mean Square of Absolute Error')
pl.legend()
pl.show()
