# $Id$
#
# compare_sam_clubb_Skewness

# Description:
#   Compares SAM and CLUBB skewness of w, rt, thl
#   As of 30-Mar 2015, only plots SAM's profiles.

import numpy as np
import matplotlib.pyplot as plt
import netCDF4
import sys

case = 'ADG1'
out_dir = '/home/weberjk/Ticket764/output/%s/'%(case)
sam_file = '/home/weberjk/Ticket764/input/input_fields/LBA.nc'
clubb_file = '%s/lba_zt.nc'%(out_dir)

Skw_denom_coef = 1.
w_tol = 2.*10**-2
thl_tol = 1.*10**-2
rt_tol = 1.*10**-8

z0 = 0 # [m]
z1 = 6000 

t0 = 189 # [min]
t1 = 249

t0_in_s = t0*60. # CLUBB's time is in seconds.
t1_in_s = t1*60.

out_name = '%dmin_%dmin_skw_compare.png'%(t0, t1)

def pull_profiles(nc, varname, conversion):
    # This function retrieves the profiles from SAM and supplies a unit conversion
    # if nessesary
    var = nc.variables[varname]
    var = np.squeeze(var)
    var = var*conversion
    return var

def return_mean_profiles(var, idx_t0, idx_t1, idx_z0, idx_z1):
    # This function returns the mean profiles over a specified interval and altitude
    var = np.mean(var[idx_t0:idx_t1,idx_z0:idx_z1],axis=0)
    return var
    
def return_skewness(xp3,xp2,tol):
    # This is CLUBB's formulation for Skewness. 
    return xp3 / ( (xp2 + Skw_denom_coef * (tol**2))**(3./2.) )

# Retrieve SAM's time and altitude arrays
nc = netCDF4.Dataset(sam_file)
sam_time = pull_profiles(nc, 'time',1.)
sam_z = pull_profiles(nc, 'z',1.)

# Slice arrays according to user's input
idx_z0 = (np.abs(sam_z[:] - z0)).argmin()
idx_z1 = (np.abs(sam_z[:] - z1)).argmin()
sam_z = sam_z[idx_z0:idx_z1]

idx_t0 = (np.abs(sam_time[:] - t0)).argmin()
idx_t1 = (np.abs(sam_time[:] - t1)).argmin()

# Retrieve SAM's mean cloud liquid water. To juxtapose with the skewness profiles 
sam_qcl = return_mean_profiles(pull_profiles(nc, 'QCL',1.),idx_t0,idx_t1,idx_z0,idx_z1)

# Retrieve the second and third order moments.
sam_thlp2 = pull_profiles(nc, 'THEL2',1.)
sam_thlp3 = pull_profiles(nc, 'THEL3',1.)
sam_Skthl = return_mean_profiles(return_skewness(sam_thlp3,sam_thlp2,thl_tol),idx_t0, idx_t1, idx_z0, idx_z1)

sam_rtp2 = pull_profiles(nc, 'QTO2',(1./1000.)**2)
# The conversion (1./1000.)**2 to get QTO3 to kg/kg is not a bug. In SAM's statistics.F90,
# to convert QTO3 to g/kg, it was multiplied by 1e6.
sam_rtp3 = pull_profiles(nc, 'QTO3',(1./1000.)**2)
sam_Skrt = return_mean_profiles(return_skewness(sam_rtp3,sam_rtp2,rt_tol),idx_t0, idx_t1, idx_z0, idx_z1)

sam_wp2 = pull_profiles(nc, 'W2',1.)
sam_wp3 = pull_profiles(nc, 'W3',1.)
sam_Skw = return_mean_profiles(return_skewness(sam_wp3,sam_wp2,w_tol),idx_t0, idx_t1, idx_z0, idx_z1)

nc.close()


#print "Plot"
f, ax1 = plt.subplots(1)
ax1.set_ylim([z0,z1])
ax2 = ax1.twiny()

ax1.plot(sam_Skw,sam_z,lw=2,label='$Skw$')
ax1.plot(sam_Skthl,sam_z,lw=2,label='$Sk_{\\theta_{l}}$')
ax1.plot(sam_Skrt,sam_z,lw=2,label='$Sk_{r_{t}}$')
ax1.set_xlabel('$Sk$')
ax1.legend()

ax2.plot(sam_qcl,sam_z,c='k',lw=2,label='$\overline{r_{c}}$')
ax2.set_xlabel('$g\ kg^{-1}$')

plt.grid()
plt.show()
plt.savefig(out_dir+out_name)
