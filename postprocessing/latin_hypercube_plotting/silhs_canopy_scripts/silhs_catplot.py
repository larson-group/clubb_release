import netCDF4
import numpy as np
import pylab as pl

clubb_nc     = netCDF4.Dataset('rico_silhs_zt.nc')
silhs_sfc_nc = netCDF4.Dataset('rico_silhs_lh_sfc.nc')

cf1 = clubb_nc.variables['cloud_frac_1']
cf2 = clubb_nc.variables['cloud_frac_2']
pf1 = clubb_nc.variables['precip_frac_1']
pf2 = clubb_nc.variables['precip_frac_2']
mxf = clubb_nc.variables['mixt_frac']

k_lh_start = silhs_sfc_nc.variables['k_lh_start']

time1 = 3000
time2 = 4320

vars_plt = list()
for i in range(0,8):
    vars_plt.append(np.empty(time2-time1))

category_labels = ['c_p_1', 'c_p_2', 'c_np_1', 'c_np_2', 'nc_p_1', 'nc_p_2', 'nc_np_1', 'nc_np_2']

for t in range(time1,time2):
    k = int(round(k_lh_start[t,0,0,0])) - 1
    cf1_t = cf1[t,k,0,0]
    cf2_t = cf2[t,k,0,0]
    pf1_t = pf1[t,k,0,0]
    pf2_t = pf2[t,k,0,0]
    mxf_t = mxf[t,k,0,0]
    vars_plt[0][t-time1] = cf1_t*pf1_t*mxf_t
    vars_plt[1][t-time1] = cf2_t*pf2_t*(1-mxf_t)
    vars_plt[2][t-time1] = cf1_t*(1-pf1_t)*mxf_t
    vars_plt[3][t-time1] = cf2_t*(1-pf2_t)*(1-mxf_t)
    vars_plt[4][t-time1] = (1-cf1_t)*pf1_t*mxf_t
    vars_plt[5][t-time1] = (1-cf2_t)*pf2_t*(1-mxf_t)
    vars_plt[6][t-time1] = (1-cf1_t)*(1-pf1_t)*mxf_t
    vars_plt[7][t-time1] = (1-cf2_t)*(1-pf2_t)*(1-mxf_t)
pl.figure()
for u in range(0,len(vars_plt)):
    if (u==7):
        pl.plot(range(time1,time2), vars_plt[u][:], 'orange', label=category_labels[u])
    else:
        pl.plot(range(time1,time2), vars_plt[u][:], label=category_labels[u])

pl.show()
pl.legend()

clubb_nc.close()
silhs_sfc_nc.close()