import netCDF4
import numpy as np
import pylab as pl

clubb_nc     = netCDF4.Dataset('cldwt/rico_silhs_zt.nc')

silhs_files = [ 'cldwt/rico_silhs_lh_zt.nc' ]
silhs_labels = [ 'cldwt' ]
silhs_sfc_nc = netCDF4.Dataset('cldwt/rico_silhs_lh_sfc.nc')

silhs_ncs = list()
for silhs_file in silhs_files:
    silhs_ncs.append(netCDF4.Dataset(silhs_file))

#############
silhs_2D_u = netCDF4.Dataset('cldwt/rico_silhs_u_lh_sample_points_2D.nc')
silhs_2D_nl = netCDF4.Dataset('cldwt/rico_silhs_nl_lh_sample_points_2D.nc')
dp1 = silhs_2D_u.variables['dp1']
rr_nl = silhs_2D_nl.variables['rr']
mf1 = clubb_nc.variables['mixt_frac']
#############

clubb_var = clubb_nc.variables['mu_rr_1']
l_time_shift = False

silhs_vars = list()
for silhs_nc in silhs_ncs:
    silhs_vars.append(silhs_nc.variables['lh_rrm'])

k_lh_start = silhs_sfc_nc.variables['k_lh_start']

time1 = 300
time2 = 400

clubb_var_plt = np.empty(time2-time1)
silhs_vars_plt = list()
for silhs_var in silhs_vars:
    silhs_vars_plt.append(np.empty(time2-time1))

for t in range(time1,time2):
    k = int(round(k_lh_start[t,0,0,0])) - 1
    if l_time_shift:
        clubb_var_plt[t-time1] = clubb_var[t-1,k,0,0]
    else:
        clubb_var_plt[t-time1] = clubb_var[t,k,0,0]

#    for u in range(0,len(silhs_vars_plt)):
#        silhs_vars_plt[u][t-time1] = silhs_vars[u][t,k,0,0]
    samples = []
    for i in range(0,100):
        if dp1[t,k,i,0] < mf1[t,k,0,0] and rr_nl[t,k,i,0] > 0:
            samples.append(rr_nl[t,k,i,0])
    avg = np.average(samples)
    silhs_vars_plt[0][t-time1] = avg

pl.plot(range(time1,time2), clubb_var_plt[:], label='analytic')
for u in range(0,len(silhs_vars_plt)):
    pl.plot(range(time1,time2), silhs_vars_plt[u][:], label=silhs_labels[u])

pl.show()
pl.legend()

clubb_nc.close()
for silhs_nc in silhs_ncs:
    silhs_nc.close()
silhs_sfc_nc.close()