#!/usr/bin/env python
import netCDF4
import sys
import numpy
import matplotlib
import matplotlib.pyplot as plt

silhs_dir = sys.argv[1]

l_print = True
l_plot = True

k_lh_start = netCDF4.Dataset(silhs_dir+'/silhs_256_1/rico_lh_lh_sfc.nc') \
             .variables['k_lh_start'][:,0,0,0]
n_timesteps = len(k_lh_start)

samp_fracs = numpy.empty((n_timesteps,8))
for i in range(8):
    samp_frac_var = netCDF4.Dataset(silhs_dir+'/silhs_256_1/rico_lh_lh_zt.nc') \
                    .variables['lh_samp_frac_'+str(i+1)][:,:,0,0]
    for t in range(0,n_timesteps):
        samp_fracs[t,i] = samp_frac_var[t,int(k_lh_start[t])-1]

if l_print:
    samp_fracs_time_avg = numpy.average(samp_fracs, axis=0)
    for i in range(8):
        print(samp_fracs_time_avg[i])
if l_plot:
    labels = ['c_p_1', 'c_p_2', 'nc_p_1', 'nc_p_2', 'c_np_1', 'c_np_2', 'nc_np_1', 'nc_np_2']
    colors = ['red', 'cyan', 'magenta', 'yellow', 'black', 'orange', 'brown', 'pink']
    for i in range(8):
        plt.plot(range(n_timesteps), samp_fracs[:,i], colors[i], label=labels[i])
    plt.legend()
    plt.show()
