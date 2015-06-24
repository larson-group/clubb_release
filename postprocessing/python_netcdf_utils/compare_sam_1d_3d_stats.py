# $Id$
"""
compare_sam_1d_3d_stats

Description: Reads in 1D statistics and 3D files from SAM and checks for consistency 
            between thlp2 and thlp3. 


"""
# Import modules
import numpy as np
import matplotlib.pyplot as plt
import netCDF4

stat_file = '/home/weberjk/SAM_CLUBB/OUT_STAT/DYCOMS_RF02_16x16x32_mpi.nc'
threeD_file = '/home/weberjk/SAM_CLUBB/OUT_3D/DYCOMS_RF02_16x16x32_mpi_4_0000000600_micro.nc' 

moment = 3 # either 2 or 3
threeD_var_name = 'THL'
oneD_var_name = 'THEL%d'%(moment)

z0 = 775 # [m]

# Clear plots
plt.close('all')

#----------------------------------------------------------------------------------------
# Functions
#----------------------------------------------------------------------------------------
def pull_data(nc, varname):
    """
    Retrieves data and unit
    """
    
    var = nc.variables[varname]
    var_u = nc.variables[varname].units
    var = np.squeeze(var)
    return var, var_u
    
def calc_moment_at_level(slab,order):
    """
    User supplies horizontal slab and the order of the moment to be computed. Returns
    the statistical moment of that altitude.
    """

    # Similar to SAM. Creates a temp variable and continually adds to it.
    moment = 0.

    # Find the dimensions of the slab
    ny = np.size(slab,axis=0)
    nx = np.size(slab,axis=1)
    
    # I trust numpy to calculate the mean
    mean = np.mean(slab)
    
    # SAM's factor_xy
    factor_xy = 1. / (nx*ny)

    # Loop through the slab, finding the deviations^moment
    for i in np.arange(0,nx):
        for j in np.arange(0,ny):
            moment = moment + (slab[i,j] - mean)**order
            
    # Return similarly to SAM        
    return moment*factor_xy
    
#----------------------------------------------------------------------------------------
# Begin Code
#----------------------------------------------------------------------------------------

# Grab 3D data
nc = netCDF4.Dataset(threeD_file)                         # Open file object
z, z_u = pull_data(nc, 'z')                              # Get altitude
time_3D_file, time_3D_unit = pull_data(nc, 'time')        # Get time
threeD_var, threeD_var_u = pull_data(nc, threeD_var_name)  # Get variate
threeD_var = np.swapaxes(threeD_var,1,2)                  # Swap axes so in z,x,y coordinates
nc.close()                                              # Close file object



# Grab 1D data
nc = netCDF4.Dataset(stat_file)                          # Open file object
oneD_var, oneD_var_u = pull_data(nc, oneD_var_name)      # Get SAM's stat file var
time_1D_file, time_1D_unit = pull_data(nc, 'time')       # Get time
nc.close()                                             # Close file object

# 3D file stores time in days. The 1D file stores time in minutes
days_to_min = 1440

# Find the stat file index the matches the 3D output
idx_time = np.abs(time_1D_file - (time_3D_file*days_to_min)).argmin()
oneD_var = oneD_var[idx_time,:]

# Compute the statistic from the 3D file
threeD_moment = np.empty_like(z)
for k in np.arange(0,len(z)):
    threeD_moment[k] = calc_moment_at_level(threeD_var[k,:,:],moment)

# Plot profiles
f, ax = plt.subplots()
ax2 = ax.twiny()
ax.plot(oneD_var,z,lw=3,c='k',label='Stat File')
ax.plot(threeD_moment,z,c='r',label='3D-Chunk Statistic')
ax.set_xlabel('$\overline{K\'^{%d}}$'%(moment), color='r')
for tl in ax.get_xticklabels():
    tl.set_color('r')

ax2.plot(oneD_var-threeD_moment,z,c='b',lw=2,ls='--',label='1D-3D Diff')
ax2.set_xlabel('1D file - 3D stats for %s'%(threeD_var_name), color='b')
for tl in ax2.get_xticklabels():
    tl.set_color('b')
    
ax.set_ylabel('z [m]')

ax.legend()
ax.grid()
plt.show()


# For visual comparison, create plan view plot
x = np.arange(0,1600,100)   
y = np.arange(0,1600,100)  
    
idx_z0 = (np.abs(z[:] - z0)).argmin()
anom = threeD_var[idx_z0,:,:] - np.mean(threeD_var[idx_z0,:,:])
f, ax = plt.subplots()
im = ax.contourf(x,y,anom,cmap=plt.get_cmap('bwr'))
cbar = f.colorbar(im,ax=ax)
cbar.set_label(threeD_var_u)
ax.set_xlabel('x [m]')
ax.set_ylabel('y [m]')
ax.set_title('%s\' at %f m'%(threeD_var_name,z[idx_z0]))
plt.show()
