import netCDF4
import numpy as np
import matplotlib.pyplot as pl

time1 = 0
time2 = 864

clubb_var_str = 'rrm_accr'
silhs_var_str = 'lh_rrm_accr'

clubb_nc = netCDF4.Dataset('2Cat-Cld/silhs_1/rico_lh_zt.nc')

altitude = clubb_nc.variables['altitude']
altlow = 1
for althigh in range(2,len(altitude)):
    if altitude[althigh] > 4000.0:
        break

clubb_var_profile = np.average(clubb_nc.variables[clubb_var_str][:,altlow:althigh,0,0], axis=0)

for i in range(1,11): # Integers 1-10, inclusive
    silhs_nc_1 = netCDF4.Dataset('2Cat-Cld/silhs_'+str(i)+'/rico_lh_lh_zt.nc')
    silhs_nc_2 = netCDF4.Dataset('2Cat-CldPcp/silhs_'+str(i)+'/rico_lh_lh_zt.nc')
    silhs_nc_3 = netCDF4.Dataset('old-old/silhs_'+str(i)+'/rico_lh_lh_zt.nc')
    silhs_var_profile_1 = np.average(silhs_nc_1.variables[silhs_var_str][:,altlow:althigh,0,0], axis=0)
    silhs_var_profile_2 = np.average(silhs_nc_2.variables[silhs_var_str][:,altlow:althigh,0,0], axis=0)
    silhs_var_profile_3 = np.average(silhs_nc_3.variables[silhs_var_str][:,altlow:althigh,0,0], axis=0)

    fig = pl.figure()
    pl.plot(clubb_var_profile[:], altitude[altlow:althigh], 'k-', label='analytic')
    pl.plot(silhs_var_profile_1, altitude[altlow:althigh], 'b-', label='2Cat-Cld')
    pl.plot(silhs_var_profile_2, altitude[altlow:althigh], 'r-', label='2Cat-CldPcp')
    pl.plot(silhs_var_profile_3, altitude[altlow:althigh], 'm-', label='old-old')

    pl.legend()
    pl.savefig('plot_output/'+str(i)+'.svg')
