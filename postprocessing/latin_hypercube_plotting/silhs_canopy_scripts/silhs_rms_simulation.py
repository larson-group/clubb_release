import netCDF4
import numpy as np
import pylab as pl

clubb_nc = netCDF4.Dataset('')
silhs_nc = netCDF4.Dataset('')

clubb_var_str = 'rrm_auto'
silhs_var_str = 'lh_rrm_auto'

clubb_var = clubb_nc.variables[clubb_var_str]
silhs_var = silhs_nc.variables[silhs_var_str]

k_lh_start = netCDF4.Dataset('').variables['k_lh_start']

rms_val = 0.0
n_timesteps = clubb_var.shape[1]

for t in range(0,n_timesteps):
    k = int(round(k_lh_start[t,0,0,0]))-1
    rms_val = rms_val + (clubb_var[t,k,0,0]-silhs_var[t,k,0,0])**2

rms_val = rms_val / n_timesteps
rms_val = np.sqrt(rms_val)

print("Average RMS absolute error: " + str(rms_val))
