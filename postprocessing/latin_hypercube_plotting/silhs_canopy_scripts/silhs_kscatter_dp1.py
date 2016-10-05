import netCDF4
import numpy as np
import matplotlib.pyplot as plt

uni_nc = netCDF4.Dataset('rico_silhs_u_lh_sample_points_2D.nc')
silhs_sfc_nc = netCDF4.Dataset('rico_silhs_lh_sfc.nc')

dp1 = uni_nc.variables['dp1']
k_lh_start = silhs_sfc_nc.variables['k_lh_start']

time = 4000

k = int(round(k_lh_start[time,0,0,0])) - 1
print ('k = ',k)
points = dp1[time,k,:,0]

ones = np.empty(100)
ones[:] = 1
plt.scatter(points, ones)
plt.show()

uni_nc.close()
silhs_sfc_nc.close()