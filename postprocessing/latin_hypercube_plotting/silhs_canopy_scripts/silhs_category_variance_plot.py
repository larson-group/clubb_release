import netCDF4
import numpy as np
import pylab as pl

num_importance_categories = 8

nc     = netCDF4.Dataset('rico_lh_lh_zt.nc')
sfc_nc = netCDF4.Dataset('rico_lh_lh_sfc.nc')

cat_labels = ['c_p_1', 'c_p_2', 'nc_p_1', 'nc_p_2', 'c_np_1', \
              'c_np_2', 'nc_np_1', 'nc_np_2']

silhs_var_cat = list()
for i in range(1,num_importance_categories+1):
  silhs_var_cat.append(nc.variables['silhs_var_cat_'+str(i)])

k_lh_start = sfc_nc.variables['k_lh_start']

tstart = 0
tend   = 864

variance_fractions = np.zeros((num_importance_categories,tend-tstart))
averages = np.zeros(num_importance_categories)

n_timesteps = [0] * num_importance_categories

for t in range(tstart,tend):

  # Determine k_lh_start
  k = int(k_lh_start[t,0,0,0]) - 1

  # Find the total variance
  total_variance = 0.0
  for i in range(0,num_importance_categories):
    varnce = silhs_var_cat[i][t,k,0,0]
    if varnce >= 0.0:
      total_variance = total_variance + varnce

  # Determine variance in each category as a fraction
  if (total_variance > 0.0):
    for i in range(0,num_importance_categories):
      varnce = silhs_var_cat[i][t,k,0,0]
      if varnce >= 0.0:
        f = varnce / total_variance
        averages[i] = averages[i] + f
        n_timesteps[i] += 1
      else:
        f = 0.0
      variance_fractions[i,t-tstart] = f

averages = np.divide(averages, n_timesteps)
print(averages)

# Plots!!!
for i in range(0,num_importance_categories):
  pl.plot(range(tstart,tend), variance_fractions[i,:], label=cat_labels[i])

pl.legend()
pl.show()
