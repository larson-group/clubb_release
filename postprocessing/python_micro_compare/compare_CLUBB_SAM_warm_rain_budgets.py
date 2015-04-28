# $Id$
#
# compare_CLUBB_SAM_warm_rain_budgets

# Description:
#   Compares SAM and CLUBB rrm 'bulk' budgets. That is, turbulent advection,
#   sedimentation, etc. for Morrison Microphysics. Only plots the major terms.
import numpy as np
import matplotlib.pyplot as plt
import netCDF4

case = 'ADG1'
out_dir = '/home/weberjk/Ticket764/output/%s/'%(case)
sam_file = '/home/weberjk/Ticket764/input/input_fields/LBA.nc'
clubb_file = '%s/lba_zt.nc'%(out_dir)

z0 = 0 # [m]
z1 = 8000 

t0 = 189 # [min]
t1 = 249

t0_in_s = t0*60. # CLUBB's time is in seconds.
t1_in_s = t1*60.

out_name = '%s_%dmin_%dmin_warm_rain.png'%(case, t0, t1)

def pull_profiles(nc, varname, conversion):
    var = nc.variables[varname]
    var = np.squeeze(var)
    var = var*conversion
    return var

def return_mean_profiles(var, idx_t0, idx_t1, idx_z0, idx_z1):
    var = np.mean(var[idx_t0:idx_t1,idx_z0:idx_z1],axis=0)
    return var

print "Grab SAM's profiles"
nc = netCDF4.Dataset(sam_file)
sam_time = pull_profiles(nc, 'time',1.)
sam_z = pull_profiles(nc, 'z',1.)

idx_z0 = (np.abs(sam_z[:] - z0)).argmin()
idx_z1 = (np.abs(sam_z[:] - z1)).argmin()
sam_z = sam_z[idx_z0:idx_z1]

idx_t0 = (np.abs(sam_time[:] - t0)).argmin()
idx_t1 = (np.abs(sam_time[:] - t1)).argmin()
sam_rho = pull_profiles(nc, 'RHO',1.)

gdaym1_kgsm1 = (1.*10**-3) / 86400
sam_nr_conv = (1./sam_rho)*( (1.0*10**6)/86400.0 )

# Rain mass
sam_adv = return_mean_profiles(pull_profiles(nc, 'QRADV',gdaym1_kgsm1),idx_t0,idx_t1,idx_z0,idx_z1)
sam_diff = return_mean_profiles(pull_profiles(nc, 'QRDIFF',gdaym1_kgsm1),idx_t0,idx_t1,idx_z0,idx_z1)
sam_lsadv = return_mean_profiles(pull_profiles(nc, 'QRLSADV',gdaym1_kgsm1),idx_t0,idx_t1,idx_z0,idx_z1)
sam_mc = return_mean_profiles(pull_profiles(nc, 'QRMPHY',gdaym1_kgsm1),idx_t0,idx_t1,idx_z0,idx_z1)
sam_sed = return_mean_profiles(pull_profiles(nc, 'QRSED',gdaym1_kgsm1),idx_t0,idx_t1,idx_z0,idx_z1)
sam_sto = return_mean_profiles(pull_profiles(nc, 'QRSTO',gdaym1_kgsm1),idx_t0,idx_t1,idx_z0,idx_z1)

sam_nadv = return_mean_profiles(pull_profiles(nc, 'NRADV',sam_nr_conv),idx_t0,idx_t1,idx_z0,idx_z1)
sam_ndiff = return_mean_profiles(pull_profiles(nc, 'NRDIFF',sam_nr_conv),idx_t0,idx_t1,idx_z0,idx_z1)
sam_nlsadv = return_mean_profiles(pull_profiles(nc, 'NRLSADV',sam_nr_conv),idx_t0,idx_t1,idx_z0,idx_z1)
sam_nmc = return_mean_profiles(pull_profiles(nc, 'NRMPHY',sam_nr_conv),idx_t0,idx_t1,idx_z0,idx_z1)
sam_nsed = return_mean_profiles(pull_profiles(nc, 'NRSED',sam_nr_conv),idx_t0,idx_t1,idx_z0,idx_z1)
sam_nsto = return_mean_profiles(pull_profiles(nc, 'NRSTO',sam_nr_conv),idx_t0,idx_t1,idx_z0,idx_z1)
nc.close()

print "Grab CLUBB's profiles"
nc = netCDF4.Dataset(clubb_file)
clb_time = pull_profiles(nc, 'time',1.)
clb_z = pull_profiles(nc, 'altitude',1.)

idx_z0 = (np.abs(clb_z[:] - z0)).argmin()
idx_z1 = (np.abs(clb_z[:] - z1)).argmin()
clb_z = clb_z[idx_z0:idx_z1]

idx_t0 = (np.abs(clb_time[:] - t0_in_s)).argmin()
idx_t1 = (np.abs(clb_time[:] - t1_in_s)).argmin()

nc = netCDF4.Dataset(clubb_file)
clb_adv = return_mean_profiles(pull_profiles(nc, 'rrm_ta',1.),idx_t0,idx_t1,idx_z0,idx_z1)
clb_ts = return_mean_profiles(pull_profiles(nc, 'rrm_ts',1.),idx_t0,idx_t1,idx_z0,idx_z1)
clb_hf = return_mean_profiles(pull_profiles(nc, 'rrm_hf',1.),idx_t0,idx_t1,idx_z0,idx_z1)
clb_wvhf = return_mean_profiles(pull_profiles(nc, 'rrm_wvhf',1.),idx_t0,idx_t1,idx_z0,idx_z1)
clb_cl = return_mean_profiles(pull_profiles(nc, 'rrm_cl',1.),idx_t0,idx_t1,idx_z0,idx_z1)
clb_mc = return_mean_profiles(pull_profiles(nc, 'rrm_mc',1.),idx_t0,idx_t1,idx_z0,idx_z1)
clb_sed = return_mean_profiles(pull_profiles(nc, 'rrm_sd_morr',1.),idx_t0,idx_t1,idx_z0,idx_z1)
clb_sto = return_mean_profiles(pull_profiles(nc, 'rrm_bt',1.),idx_t0,idx_t1,idx_z0,idx_z1)

clb_mc = clb_mc - clb_sed
clb_sed = clb_sed +clb_ts
clb_clip = clb_hf + clb_wvhf + clb_cl

clb_nadv = return_mean_profiles(pull_profiles(nc, 'Nrm_ta',1.),idx_t0,idx_t1,idx_z0,idx_z1)
clb_nts = return_mean_profiles(pull_profiles(nc, 'Nrm_ts',1.),idx_t0,idx_t1,idx_z0,idx_z1)
clb_ncl = return_mean_profiles(pull_profiles(nc, 'Nrm_cl',1.),idx_t0,idx_t1,idx_z0,idx_z1)
clb_nmc = return_mean_profiles(pull_profiles(nc, 'Nrm_mc',1.),idx_t0,idx_t1,idx_z0,idx_z1)
clb_nsed = return_mean_profiles(pull_profiles(nc, 'NRSTEN',1.),idx_t0,idx_t1,idx_z0,idx_z1)
clb_nsto = return_mean_profiles(pull_profiles(nc, 'Nrm_bt',1.),idx_t0,idx_t1,idx_z0,idx_z1)

clb_nmc = clb_nmc - clb_nsed
clb_nsed = clb_nsed + clb_nts
nc.close()


print "Plot"
f, (ax1, ax2) = plt.subplots(1,2,sharey=True)
ax1.set_ylim([z0,z1])
#SAM
ax1.plot(sam_adv,sam_z, c='green', linewidth=2, ls='--')
ax1.plot(sam_mc,sam_z, c='blue', linewidth=2, ls='--')
ax1.plot(sam_sed,sam_z, c='red', linewidth=2, ls='--')
ax1.plot(sam_sto,sam_z, c='k', linewidth=2, ls='--')

#CLUBB
ax1.plot(clb_adv,clb_z, c='green', linewidth=2,label='Turbulent Advection' )
ax1.plot(clb_mc,clb_z, c='blue', linewidth=2 ,label='Microphysics')
ax1.plot(clb_sed,clb_z,c='red', linewidth=2 ,label='Sedimentation')
ax1.plot(clb_sto,clb_z,c='k', linewidth=2 ,label='Time Tendency')
ax1.legend()


#SAM
ax2.plot(sam_nadv,sam_z, c='green', linewidth=2, ls='--')
ax2.plot(sam_nmc,sam_z, c='blue', linewidth=2, ls='--')
ax2.plot(sam_nsed,sam_z,c='red', linewidth=2, ls='--')
ax2.plot(sam_nsto,sam_z, c='k', linewidth=2, ls='--')

#CLUBB
ax2.plot(clb_nadv,clb_z, c='green', linewidth=2)
ax2.plot(clb_nmc,clb_z, c='blue', linewidth=2)
ax2.plot(clb_nsed,clb_z,c='red', linewidth=2)
ax2.plot(clb_nsto,clb_z,c='k', linewidth=2 )

ax1.grid()
ax2.grid()

ax1.set_xlabel(r'$kg\ kg^{-1}\ s^{-1}$')
ax2.set_xlabel(r'$\#\ kg^{-1}\ s^{-1}$')
ax1.set_ylabel(r'$m$')
plt.savefig(out_dir+out_name)
plt.close()
