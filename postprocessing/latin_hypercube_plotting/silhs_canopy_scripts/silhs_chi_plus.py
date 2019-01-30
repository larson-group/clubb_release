import netCDF4
import numpy as np
import pylab as pl

silhs_sfc_nc = netCDF4.Dataset('rico_silhs_lh_sfc.nc')
silhs_2D_nl_nc = netCDF4.Dataset('rico_silhs_nl_lh_sample_points_2D.nc')

num_samps = num_samps = len(silhs_2D_nl_nc.dimensions['latitude'])

k_lh_start = silhs_sfc_nc.variables['k_lh_start']

chi_2D = silhs_2D_nl_nc.variables['chi']

time1 = 3000
time2 = 4320

chi_plus = np.zeros(time2-time1)
num_pos = 0

for time in range(time1,time2):
    print(time)
    t = time-time1
    k = int(round(k_lh_start[time,0,0,0])) - 1
    for s in range(0,num_samps):
        if (chi_2D[time,k,s,0] > 0):
            num_pos = num_pos+1
            chi_plus[t] = chi_plus[t] + chi_2D[time,k,s,0]
    if (num_pos != 0):
        chi_plus[t] = chi_plus[t] / num_pos
    

pl.plot(range(time1,time2), chi_plus[:], label="analytic")
pl.show()
pl.legend()

silhs_2D_nl_nc.close()