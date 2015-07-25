import netCDF4
import numpy

def pdf_p_at_k_lh_start(prefix, t):
    clubb_nc = netCDF4.Dataset(prefix+'_zt.nc')
    silhs_sfc_nc = netCDF4.Dataset(prefix+'_lh_sfc.nc')
    k = int(silhs_sfc_nc.variables['k_lh_start'][t,0,0,0])-1

    cf1 = clubb_nc.variables['cloud_frac_1'][t,k,0,0]
    cf2 = clubb_nc.variables['cloud_frac_2'][t,k,0,0]
    pf1 = clubb_nc.variables['precip_frac_1'][t,k,0,0]
    pf2 = clubb_nc.variables['precip_frac_2'][t,k,0,0]
    mf1 = clubb_nc.variables['mixt_frac'][t,k,0,0]
    mf2 = 1.0 - mf1

    p = numpy.empty(8)
    p[0] = cf1*pf1*mf1
    p[1] = cf2*pf2*mf2
    p[2] = (1-cf1)*pf1*mf1
    p[3] = (1-cf2)*pf2*mf2
    p[4] = cf1*(1-pf1)*mf1
    p[5] = cf2*(1-pf2)*mf2
    p[6] = (1-cf1)*(1-pf1)*mf1
    p[7] = (1-cf2)*(1-pf2)*mf2

    clubb_nc.close()
    silhs_sfc_nc.close()
    return p
