#
# Compare turbulent transport
#
import numpy as np
import matplotlib.pyplot as plt
import netCDF4

test = ['Control','NonLcl','const_lscale']
case = 'LBA'

sam_vars = ['covarnce_w_rr','covarnce_w_rg','covarnce_w_rs','covarnce_w_ri',
            'covarnce_w_Nr','covarnce_w_Ng','covarnce_w_Ns','covarnce_w_Ni']

#clubb_vars = ['rrm_ta','rgm_ta','rsm_ta','rim_ta',
#             'Nrm_ta','Ngm_ta','Nsm_ta','Nim_ta']

clubb_vars = ['wprrp','wprgp','wprsp','wprip',
             'wpNrp','wpNgp','wpNsp','wpNip']

out_dir = '/home/weberjk/precip_ta/output/'
sam_file = '/home/weberjk/precip_ta/input/input_fields/LBA.nc'

z0 = 0 # [m]
z1 = 18000 

t0 = 189 # [min]
t1 = 360

t0_in_s = t0*60. # CLUBB's time is in seconds.
t1_in_s = t1*60.

#----------------------------------------------------------------------------------------
# Functions
#----------------------------------------------------------------------------------------

def pull_profiles(nc, varname, conversion):
    # This function retrieves the profiles from SAM and performs a unit conversion
    # if nessesary
    var = nc.variables[varname]
    var = np.squeeze(var)
    var = var*conversion
    return var

def return_mean_profiles(var, idx_t0, idx_t1, idx_z0, idx_z1):
    # This function returns the mean profiles over a specified interval and altitude
    var = np.mean(var[idx_t0:idx_t1,idx_z0:idx_z1],axis=0)
    return var

#----------------------------------------------------------------------------------------
# Retrieve SAM's altitude, time, and flux profiles
#----------------------------------------------------------------------------------------
nc = netCDF4.Dataset(sam_file)
sam_z = pull_profiles(nc, 'z', 1.)
sam_t = pull_profiles(nc, 'time', 1.)
sam_rho = pull_profiles(nc, 'RHO', 1.) # For conversions
sam_nr_conv = (sam_rho*1000.)**-1 

idx_z0 = (np.abs(sam_z[:] - z0)).argmin()
idx_z1 = (np.abs(sam_z[:] - z1)).argmin()
sam_z = sam_z[idx_z0:idx_z1]

idx_t0 = (np.abs(sam_t[:] - t0)).argmin()
idx_t1 = (np.abs(sam_t[:] - t1)).argmin()

# Create array to hold all the profiles
sam_flux = np.empty( (len(sam_vars),len(sam_z) ) )

for j in np.arange(0,len(sam_vars)):
    print "Pulling SAM variable: %s"%(sam_vars[j])
    sam_flux[j,:] = ( return_mean_profiles(
                             pull_profiles(nc, sam_vars[j],1.),
                                      idx_t0, idx_t1, idx_z0, idx_z1) )
nc.close()

#----------------------------------------------------------------------------------------
# Retrieve CLUBB's altitude, time, and flux profiles
#----------------------------------------------------------------------------------------
clubb_file = '%s/lba_zm.nc'%(out_dir+test[0])
nc = netCDF4.Dataset(clubb_file)
clb_z = pull_profiles(nc, 'altitude', 1.)
clb_t = pull_profiles(nc, 'time', 1.)

idx_z0 = (np.abs(clb_z[:] - z0)).argmin()
idx_z1 = (np.abs(clb_z[:] - z1)).argmin()
clb_z = clb_z[idx_z0:idx_z1]

idx_t0 = (np.abs(clb_t[:] - t0_in_s)).argmin()
idx_t1 = (np.abs(clb_t[:] - t1_in_s)).argmin()

# Create an array to hold all the profiles for each test
clubb_flux = np.empty( ( len(test),len(clubb_vars),len(clb_z) ) )

for i in np.arange(0,len(test)):
    for j in np.arange(0,len(clubb_vars)):
        print "Pulling CLUBB variable: %s from case:%s"%(clubb_vars[j], test[i])
        clubb_file = '%s/lba_zm.nc'%(out_dir+test[i])
        nc = netCDF4.Dataset(clubb_file)
        clubb_flux[i,j,:] = ( return_mean_profiles(
                             pull_profiles(nc, clubb_vars[j],1.),
                                      idx_t0, idx_t1, idx_z0, idx_z1) )
nc.close()

f,((ax1,ax2,ax3,ax4),(ax5,ax6,ax7,ax8)) = plt.subplots(2,4,sharey=True)
# First, plot SAM data
ax1.plot(sam_flux[0],sam_z,lw=2,c='k')
ax1.grid()
ax1.set_ylim(0,max(sam_z))
ax2.plot(sam_flux[1],sam_z,lw=2,c='k')
ax2.grid()
ax3.plot(sam_flux[2],sam_z,lw=2,c='k')
ax3.grid()
ax4.plot(sam_flux[3],sam_z,lw=2,c='k')
ax4.grid()
ax5.plot(sam_flux[4],sam_z,lw=2,c='k')
ax5.grid()
ax6.plot(sam_flux[5],sam_z,lw=2,c='k')
ax6.grid()
ax7.plot(sam_flux[6],sam_z,lw=2,c='k')
ax7.grid()
ax8.plot(sam_flux[7],sam_z,lw=2,c='k')
ax8.grid()

for i in np.arange(0,len(test)):
    ax1.plot(clubb_flux[i,0,:],clb_z,label=test[i])
    ax2.plot(clubb_flux[i,1,:],clb_z)
    ax3.plot(clubb_flux[i,2,:],clb_z)
    ax4.plot(clubb_flux[i,3,:],clb_z)
    ax5.plot(clubb_flux[i,4,:],clb_z)
    ax6.plot(clubb_flux[i,5,:],clb_z)
    ax7.plot(clubb_flux[i,6,:],clb_z)
    ax8.plot(clubb_flux[i,7,:],clb_z) 

ax1.legend()
plt.show()