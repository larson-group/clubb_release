# $Id$
#
# compare_precip_profiles 
# 
# Description:
#   Compares the precipitation profiles for the Morrison Microphysics parameterization


# Import libraries
import numpy as np
import matplotlib.pyplot as plt
import netCDF4

# Point to CLUBB's 'output' directory and location of SAM's stat file 
out_dir = '/home/weberjk/precip_ta/output/'
sam_file = '/home/weberjk/precip_ta/input/input_fields/LBA.nc'

# Within CLUBB's 'output' directory, specify the names of each subdirectory
# containing different simulations. Each subdirectory's data will be overplotted
# and the name of the subdirectory used in the legend.
output_subdirectories = ['Control','NonLcl','const_lscale']

# For plotting, include what case is being plotted.
case = 'LBA'


# SAM and CLUBB vars. For each SAM variable, the CLUBB equivalent must be in the same 
# array position.
sam_vars = ['QR','QG','QS','QI','NR','NG','NS','NI']
clubb_vars = ['rrm','rgm','rsm','rim','Nrm','Ngm','Nsm','Nim']


z0 = 0 # Start height [m]
z1 = 18000 # end height

t0 = 189 # Start time [min]
t1 = 360 # end time

#----------------------------------------------------------------------------------------
# Functions
#----------------------------------------------------------------------------------------

def pull_profiles(nc, varname, conversion):
    # This function takes a netcdf object, variable name, and conversion factor as input.
    # It returns the converted variable
    var = nc.variables[varname]
    var = np.squeeze(var)
    var = var*conversion
    return var

def return_mean_profiles(var, idx_t0, idx_t1, idx_z0, idx_z1):
    # This function returns the mean profile of var, spanning the time indicies idx_t0 to
    # idx_t1 and height indicies idx_z0 to idx_z1.
    var = np.mean(var[idx_t0:idx_t1,idx_z0:idx_z1],axis=0)
    return var

#----------------------------------------------------------------------------------------
# Begin Code
#----------------------------------------------------------------------------------------

t0_in_s = t0*60. # CLUBB's time is in seconds.
t1_in_s = t1*60.

#----------------------------------------------------------------------------------------
# Retrieve SAM's altitude, time, and flux profiles
#----------------------------------------------------------------------------------------
nc = netCDF4.Dataset(sam_file)
sam_z = pull_profiles(nc, 'z', 1.)
sam_t = pull_profiles(nc, 'time', 1.)
sam_rho = pull_profiles(nc, 'RHO', 1.) # For conversions
sam_nr_conv = 100**3 * sam_rho**-1 # Convert SAM's #/cm^3 to #/kg
sam_rr_conv = 1000.**-1 # Convert SAM's g/kg to kg/kg

# Store the conversions corrosponding to SAM's indicies for each variable
convert = [sam_rr_conv,sam_rr_conv,sam_rr_conv,sam_rr_conv,
          sam_nr_conv,sam_nr_conv,sam_nr_conv,sam_nr_conv]

idx_z0 = (np.abs(sam_z[:] - z0)).argmin()
idx_z1 = (np.abs(sam_z[:] - z1)).argmin()
sam_z = sam_z[idx_z0:idx_z1]

idx_t0 = (np.abs(sam_t[:] - t0)).argmin()
idx_t1 = (np.abs(sam_t[:] - t1)).argmin()

# Create a structure that will hold all the profiles for each SAM variable
sam_mean = np.empty( (len(sam_vars),len(sam_z) ) )

# Loop through all the variables
for j in np.arange(0,len(sam_vars)):
    print "Pulling SAM variable: %s"%(sam_vars[j])
    sam_mean[j,:] = ( return_mean_profiles(
                             pull_profiles(nc, sam_vars[j],convert[j]),
                                      idx_t0, idx_t1, idx_z0, idx_z1) )
nc.close()

#----------------------------------------------------------------------------------------
# Retrieve CLUBB's altitude, time, and flux profiles
#----------------------------------------------------------------------------------------
clubb_file = '%s/lba_zt.nc'%(out_dir+test[0])
nc = netCDF4.Dataset(clubb_file)
clb_z = pull_profiles(nc, 'altitude', 1.)
clb_t = pull_profiles(nc, 'time', 1.)

idx_z0 = (np.abs(clb_z[:] - z0)).argmin()
idx_z1 = (np.abs(clb_z[:] - z1)).argmin()
clb_z = clb_z[idx_z0:idx_z1]

idx_t0 = (np.abs(clb_t[:] - t0_in_s)).argmin()
idx_t1 = (np.abs(clb_t[:] - t1_in_s)).argmin()

# Create an array to hold all the clubb profiles, this time, for each subdirectory
clubb_mean = np.empty( ( len(output_subdirectories),len(clubb_vars),len(clb_z) ) )

# Loop through each subdirectory
for i in np.arange(0,len(output_subdirectories)):
    clubb_file = '%s/lba_zt.nc'%(out_dir+test[i])
    nc = netCDF4.Dataset(clubb_file)

    # Loop through each variable
    for j in np.arange(0,len(clubb_vars)):
        print "Pulling CLUBB variable: %s from case:%s"%(clubb_vars[j], test[i])
        clubb_mean[i,j,:] = ( return_mean_profiles(
                             pull_profiles(nc, clubb_vars[j],1.),
                                      idx_t0, idx_t1, idx_z0, idx_z1) )
nc.close()

#----------------------------------------------------------------------------------------
# Plot
#-----------------------------------------------------------------------------------------

f,((ax1,ax2,ax3,ax4),(ax5,ax6,ax7,ax8)) = plt.subplots(2,4,sharey=True)
# First, plot SAM data
ax1.plot(sam_mean[0],sam_z,lw=2,c='k')
ax1.grid()
ax1.set_ylim(0,max(sam_z))
ax2.plot(sam_mean[1],sam_z,lw=2,c='k')
ax2.grid()
ax3.plot(sam_mean[2],sam_z,lw=2,c='k')
ax3.grid()
ax4.plot(sam_mean[3],sam_z,lw=2,c='k')
ax4.grid()
ax5.plot(sam_mean[4],sam_z,lw=2,c='k')
ax5.grid()
ax6.plot(sam_mean[5],sam_z,lw=2,c='k')
ax6.grid()
ax7.plot(sam_mean[6],sam_z,lw=2,c='k')
ax7.grid()
ax8.plot(sam_mean[7],sam_z,lw=2,c='k')
ax8.grid()

for i in np.arange(0,len(test)):
    ax1.plot(clubb_mean[i,0,:],clb_z,label=test[i])
    ax2.plot(clubb_mean[i,1,:],clb_z)
    ax3.plot(clubb_mean[i,2,:],clb_z)
    ax4.plot(clubb_mean[i,3,:],clb_z)
    ax5.plot(clubb_mean[i,4,:],clb_z)
    ax6.plot(clubb_mean[i,5,:],clb_z)
    ax7.plot(clubb_mean[i,6,:],clb_z)
    ax8.plot(clubb_mean[i,7,:],clb_z) 

ax1.legend()
plt.show()
