import netCDF4
import numpy as np
import pylab as pl

clubb_nc     = netCDF4.Dataset('rico_silhs_zt.nc')
silhs_nc     = netCDF4.Dataset('rico_silhs_lh_zt.nc')

#clubb_var = clubb_nc.variables['rrm']
cf1 = clubb_nc.variables['cloud_frac1']
cf2 = clubb_nc.variables['cloud_frac2']
rr1 = clubb_nc.variables['rr1']
rr2 = clubb_nc.variables['rr2']
pf1 = clubb_nc.variables['precip_frac_1']
pf2 = clubb_nc.variables['precip_frac_2']
mf = clubb_nc.variables['mixt_frac']
precip_frac = clubb_nc.variables['precip_frac']
lh_cloud_frac = silhs_nc.variables['lh_cloud_frac']
lh_mixt_frac = silhs_nc.variables['lh_mixt_frac']
lh_precip_frac = silhs_nc.variables['lh_precip_frac']
lh_rrm = silhs_nc.variables['lh_rrm']

time1 = 3600
time2 = 3601

nz = len(clubb_nc.dimensions['altitude'])

clubb_var_plt = np.zeros(nz)
silhs_var_plt = np.zeros(nz)

def clubb_var(t,k,x,y):
    return rr1[t,k,0,0]*pf1[t,k,0,0]*mf[t,k,0,0]+rr2[t,k,0,0]*pf2[t,k,0,0]*\
        (1.0-mf[t,k,0,0])
    #return cf1[t,k,0,0]*mf[t,k,0,0] + (1.0-mf[t,k,0,0])*cf2[t,k,0,0]
    #return precip_frac[t,k,0,0]

def silhs_var(t,k,x,y):
    #return lh_precip_frac[t,k,0,0]
    return lh_rrm[t,k,0,0]

for k in range(0,nz):
    for t in range(time1,time2):
        clubb_var_plt[k] = clubb_var_plt[k] + clubb_var(t,k,0,0)
        silhs_var_plt[k] = silhs_var_plt[k] + silhs_var(t,k,0,0)
    clubb_var_plt[k] = clubb_var_plt[k] / (time2-time1)
    silhs_var_plt[k] = silhs_var_plt[k] / (time2-time1)

pl.plot(clubb_var_plt[:], range(0,nz), label="analytic")
pl.plot(silhs_var_plt[:], range(0,nz), label="SILHS")
pl.show()
pl.legend()

clubb_nc.close()
silhs_nc.close()