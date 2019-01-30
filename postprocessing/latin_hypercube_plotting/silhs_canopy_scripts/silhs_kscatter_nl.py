import netCDF4
import matplotlib.pyplot as plt

nl_nc = netCDF4.Dataset('rico_silhs_nl_lh_sample_points_2D.nc')
silhs_sfc_nc = netCDF4.Dataset('rico_silhs_lh_sfc.nc')

chi = nl_nc.variables['chi']
rr  = nl_nc.variables['rr']
k_lh_start = silhs_sfc_nc.variables['k_lh_start']

time = 3589

k = int(round(k_lh_start[time,0,0,0])) - 1
print ('k = ',k)
num_samps = len(nl_nc.dimensions['latitude']) # I call hack!
print ('num_samps = ', num_samps)
plt.scatter(chi[time,k,:,0], rr[time,k,:,0])
plt.show()

nl_nc.close()
silhs_sfc_nc.close()