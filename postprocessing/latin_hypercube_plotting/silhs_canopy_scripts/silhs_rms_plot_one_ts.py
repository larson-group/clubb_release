import netCDF4
import numpy as np
import pylab as pl

sim_points_all = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, \
                  1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400, 2500, 2600, \
                  2700, 2800, 2900, 3000]

clubb_var_str = 'rc_1'
silhs_var_str  = 'lh_rcm'

cnc = netCDF4.Dataset('silhs_'+str(sim_points_all[0])+'_1/rico_lh_zt.nc')
clubb_var = cnc.variables[clubb_var_str]
cf1 = cnc.variables['cloud_frac_1']
pf1 = cnc.variables['precip_frac_1']
mf  = cnc.variables['mixt_frac']

k_lh_start = netCDF4.Dataset('silhs_'+str(sim_points_all[0])+'_1/rico_lh_lh_sfc.nc').variables['k_lh_start']

rms_all = list()
rms_1_over_n = list()
rms_1_over_sqrt_n = list()

n_ensembles = 100

n_interp = 0

for n in range(0,len(sim_points_all)):
    num_samples = sim_points_all[n]
    print("n = " + str(num_samples))
    rms_val = 0.0
    k = int(round(k_lh_start[0,0,0,0]))-1
    clubb_val = clubb_var[0,k,0,0]*pf1[0,k,0,0]*mf[0,k,0,0]#*cf1[0,k,0,0]
    for j in range(1,n_ensembles+1):
        silhs_val = netCDF4.Dataset('silhs_'+str(num_samples)+'_'+str(j)+'/rico_lh_lh_zt.nc').variables[silhs_var_str][0,k,0,0]
        rms_val = rms_val + (clubb_val-silhs_val)**2
    rms_val = rms_val / n_ensembles
    rms_val = np.sqrt(rms_val)
    rms_all.append(rms_val)
    if n == n_interp:
        rms_1_over_n.append(rms_val)
        rms_1_over_sqrt_n.append(rms_val)
    elif n > n_interp:
        rms_1_over_n.append(rms_1_over_n[n-1] * (float(sim_points_all[n-1]) / sim_points_all[n]))
        rms_1_over_sqrt_n.append(rms_1_over_sqrt_n[n-1] * np.sqrt(float(sim_points_all[n-1]) / sim_points_all[n]))
    else:
        rms_1_over_n.append(0.0)
        rms_1_over_sqrt_n.append(0.0)

###########################################
# Go and compute the rest
n=n_interp-1
while n>=0:
    rms_1_over_n[n] = rms_1_over_n[n+1] * (float(sim_points_all[n+1]) / sim_points_all[n])
    rms_1_over_sqrt_n[n] = rms_1_over_sqrt_n[n+1] * np.sqrt(float(sim_points_all[n+1]) / sim_points_all[n])
    n = n - 1
###########################################

pl.figure()
pl.plot(sim_points_all, rms_all, label='SILHS')
pl.plot(sim_points_all, rms_1_over_n, label='1 over N')
pl.plot(sim_points_all, rms_1_over_sqrt_n, label='1 over sqrt N')

pl.xlabel('Number of Sample Points')
pl.ylabel('Root Mean Square of Absolute Error')
pl.xscale('log')
pl.yscale('log')
pl.axis('tight')

pl.legend()
pl.show()
