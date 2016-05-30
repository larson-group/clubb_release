# $Id: $
"""
 plotter.py 
 
 Description:
   Used to plot one case of CLUBB given a clubb output .nc file.
   Version .01, basic functionality only.  Requires editing within script
   to change between cases.
"""

# Import libraries
import numpy as np
import matplotlib.pyplot as plt
import netCDF4

# Point to .nc file within 'output' directory
clubb_ncfile = 'lba_zt.nc'

# Point to CLUBB's 'output' directory 
out_dir = '../../output'

# For plotting, include what case is being plotted.
case = 'LBA'

#CLUBB variables to plot.  Note, percentages are in 0.xx form due to nature of script.
clubb_vars = ['thlm',
'rtm',
'wpthlp2',
'wprtp2',
'cloud_frac',
'rcm',
'wp2_zt',
'wp3',
'thlp2_zt',
'rtp2_zt',
'rtpthlp_zt',
'wm',
'um',
'vm',
'upwp_zt',
'vpwp_zt',
'up2_zt',
'vp2_zt',
'wpthvp',
'Lscale',
'tau_zm',
'radht',
'precip_frac',
'ice_supersat_frac',
'rrm',
'Nrm',
'Ncm',
'Nc_in_cloud',
'rgm',
'Ngm',
'rim',
'Nim',
'rsm',
'Nsm',
'lwp',
'rwp',
'precip_rate_sfc',
'wp2_vert_avg',
'iwp',
'swp']

axis_titles = ['thlm [K]',
'rtm [kg/kg]',
'wpthlp [K m/s]',
'wprtp [(kg/kg) m/s]',
'cloud_frac [%]',
'rcm [kg/kg]',
'wp2 [m^2/s^2]',
'wp3 [m^3/s^3]',
'thlp2 [K^2]',
'rtp2 [(kg/kg)^2]',
'rtpthlp [(kg/kg) K]',
'wm [m/s]',
'um [m/s]',
'vm [m/s]',
'upwp [m^2/s^2]',
'vpwp [m^s/s^2]',
'up2 [m^s/s^2]',
'vp2 [m^2/s^2]',
'wpthvp [K m/s]',
'Lscale [m]',
'tau_zm [s]',
'radht [K/s]',
'precip_frac [%]',
'ice_supersat_frac [%]',
'rrm [kg/kg]',
'Nrm [num/kg]',
'Ncm [num/kg]',
'Nc_in_cloud [num/kg]',
'rgm [kg/kg]',
'Ngm [#/kg]',
'rim [kg/kg]',
'Nim [num/kg]',
'rsm [kg/kg]',
'Nsm [num/kg]',
'lwp [kg/m^2]',
'rwp [kg/m^2]',
'rain_rate_sfc [mm/day]',
'wp2 [m^2/s^2]',
'iwp [kg/m^2]',
'swp [kg/m^2]']

graph_titles = ['Liquid Water Potential Temperature, theta_l',
'Total Water Mixing Ratio, r_t',
'Turbulent Flux of theta_l',
'Turbulent Flux of r_t',
'Cloud Fraction',
'Cloud Water Mixing Ratio, r_c',
'Variance of w',
'Third-order Moment of w',
'Variance of theta_l',
'Variance of r_t',
'Covariance of r_t & theta_l',
'Vertical Wind Component, w (subsidence)',
'Zonal Wind Component, u',
'Meridional Wind Component, v',
'Covariance of u & w',
'Covariance of v & w',
'Variance of u wind',
'Variance of v wind',
'Buoyancy flux ',
'Mixing Length',
'Dissipation Time Scale',
'LW + SW radiative heating rate',
'Precip Fraction',
'Ice Fraction',
'Rain Water Mixing Ratio, r_r',
'Rain Drop Concentration, N_r',
'Cloud Droplet Concentration',
'In Cloud Droplet Concentration',
'Graupel Mixing Ratio',
'Graupel Number Concentration',
'Cloud Ice Mixing Ratio',
'Cloud Ice Concentration',
'Snow Mixing Ratio',
'Snow Number Concentration',
'Liquid Water Path',
'Rain Water Path',
'Surface rainfall rate',
'Density-Weighted Vertically Averaged wp2',
'Cloud Ice Water Path',
'Snow Water Path']


z0 = 0 # Start height [m]
z1 = 18000 # end height

t0 = 189 # Start time [min]
t1 = 360 # end time

#----------------------------------------------------------------------------------------
# Should not have to edit below this line
#----------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------
# Useful constants
#----------------------------------------------------------------------------------------
s_in_min = 60

#----------------------------------------------------------------------------------------
# Functions
#----------------------------------------------------------------------------------------

def pull_profiles(nc, varname, conversion):
    """
    Input:
      nc         --  Netcdf file object
      varname    --  Variable name string
      conversion --  Conversion factor

    Output:
      time x height array of the specified variable
    """

    var = nc.variables[varname]
    var = np.squeeze(var)
    var = var*conversion
    return var

def return_mean_profiles(var, idx_t0, idx_t1, idx_z0, idx_z1):
    """
    Input:
      var    -- time x height array of some property
      idx_t0 -- Index corrosponding to the beginning of the averaging interval
      idx_t1 -- Index corrosponding to the end of the averaging interval
      idx_z0 -- Index corrosponding to the lowest model level of the averaging interval
      idx_z1 -- Index corrosponding to the highest model level of the averaging interval

    Output:
      var    -- time averaged vertical profile of the specified variable
    """

    var = np.mean(var[idx_t0:idx_t1,idx_z0:idx_z1],axis=0)
    return var

#----------------------------------------------------------------------------------------
# Begin Code
#----------------------------------------------------------------------------------------

t0_in_s = t0*s_in_min # CLUBB's time is in seconds.
t1_in_s = t1*s_in_min

#----------------------------------------------------------------------------------------
# Retrieve CLUBB's altitude, time, and flux profiles
#----------------------------------------------------------------------------------------
clubb_file = '%s/%s'%(out_dir, clubb_ncfile)
nc = netCDF4.Dataset('/home/nieznan3/clubb/output/lba_zt.nc')
clb_z = pull_profiles(nc, 'altitude', 1.)
clb_t = pull_profiles(nc, 'time', 1.)

idx_z0 = (np.abs(clb_z[:] - z0)).argmin()
idx_z1 = (np.abs(clb_z[:] - z1)).argmin()
clb_z = clb_z[idx_z0:idx_z1]

idx_t0 = (np.abs(clb_t[:] - t0_in_s)).argmin()
idx_t1 = (np.abs(clb_t[:] - t1_in_s)).argmin()

# Create an array to hold all the clubb profiles, this time, for each subdirectory
clubb_mean = np.empty( ( len(clubb_vars),len(clb_z) ) )

# Loop through each variable
for j in np.arange(0,len(clubb_vars)):
    print "Pulling CLUBB variable: %s from case:"%(clubb_vars[j])
    clubb_mean[j,:] = ( return_mean_profiles(
                         pull_profiles(nc, clubb_vars[j],1.),
                                  idx_t0, idx_t1, idx_z0, idx_z1) )
nc.close()

#----------------------------------------------------------------------------------------
# Plot
#-----------------------------------------------------------------------------------------

f,((ax1,ax2,ax3,ax4),(ax5,ax6,ax7,ax8)) = plt.subplots(2,4,sharey=True)

#for i in np.arange(0,len(test)):
ax1.plot(clubb_mean[0,:],clb_z)
ax2.plot(clubb_mean[1,:],clb_z)
ax3.plot(clubb_mean[2,:],clb_z)
ax4.plot(clubb_mean[3,:],clb_z)
ax5.plot(clubb_mean[4,:],clb_z)
ax6.plot(clubb_mean[5,:],clb_z)
ax7.plot(clubb_mean[6,:],clb_z)
ax8.plot(clubb_mean[7,:],clb_z) 

ax1.legend()
plt.show()
