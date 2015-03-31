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
sigma_tilde_w_sqd = .4 # For simplicity, fix as in Larson et al. (2002)
                      # We could make dependent on gamma 
num_beta_lines = 1000


beta = np.linspace(0,3,num_beta_lines)

z0 = 0 # [m]
z1 = 6000 

t0 = 189 # [min]
t1 = 249

t0_in_s = t0*60. # CLUBB's time is in seconds.
t1_in_s = t1*60.

out_name = '%dmin_%dmin'%(t0, t1)

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
    
def return_skewness(xp3,xp2,tol):
    # This is CLUBB's formulation for Skewness. 
    return xp3 / ( (xp2 + Skw_denom_coef * (tol**2))**(3./2.) )
    
def return_nrmlz_skw(Skw,sigma_tilde_w_sqd):
    # Normalized Skw. See appendix A of Larson and Golaz (2005) 
    return Skw*(1. - sigma_tilde_w_sqd)**(-3./2.)
    
def return_nrmlz_w_corr(wpxp,wp2,xp2,w_tol,x_tol,sigma_tilde_w_sqd):
    # Normalized corr_w_x. E.g. Equation 16 of Larson and Golaz (2005) 
    fac = (1. - sigma_tilde_w_sqd)**(-1./2.)
    corr = np.empty_like(wpxp)
    
    # Must iterate over time and space. Otherwise, max(w_tol**2,wp2) throws an error
    for t in np.arange(0,np.size(corr,axis=0)):
        for i in np.arange(0,np.size(corr,axis=1)):
            corr[t,i] = (wpxp[t,i] / 
            ( np.sqrt(max(w_tol**2,wp2[t,i])) * np.sqrt(max(x_tol**2,xp2[t,i])) ))

    return fac*corr
    
def prognose_skex(sam_Skw_hat,sam_nrmlz_corr_wx,beta):
    return sam_Skw_hat*sam_nrmlz_corr_wx*(beta+(1.-beta)*sam_nrmlz_corr_wx**2)

# Retrieve SAM's time, altitude, and air density arrays
nc = netCDF4.Dataset(sam_file)
sam_time = pull_profiles(nc, 'time',1.)
sam_z = pull_profiles(nc, 'z',1.)
rho = pull_profiles(nc, 'RHO',1.)

# Reshape arrays according to user's input
idx_z0 = (np.abs(sam_z[:] - z0)).argmin()
idx_z1 = (np.abs(sam_z[:] - z1)).argmin()
sam_z = sam_z[idx_z0:idx_z1]

idx_t0 = (np.abs(sam_time[:] - t0)).argmin()
idx_t1 = (np.abs(sam_time[:] - t1)).argmin()
sam_time = sam_time[idx_t0:idx_t1]

# Retrieve SAM's mean cloud liquid water. To juxtapose with the skewness profiles 
sam_qcl = return_mean_profiles(pull_profiles(nc, 'QCL',1.),idx_t0,idx_t1,idx_z0,idx_z1)


#----------------------------------------------------------------------------------------
# Retrieve the second and third order moments.
#----------------------------------------------------------------------------------------

#Theta_l
sam_thlp2 = pull_profiles(nc, 'THEL2',1.)
sam_thlp3 = pull_profiles(nc, 'THEL3',1.)
sam_wpthlp = pull_profiles(nc, 'TLFLUX',(rho * 1004.)**-1)

#rt
sam_rtp2 = pull_profiles(nc, 'QT2',(1./1000.)**2)
# The conversion (1./1000.)**2 to get QTO3 to kg/kg is not a bug. 
# In SAM's statistics.F90, to convert QTO3 to g/kg, it was multiplied by 1e6.
sam_rtp3 = pull_profiles(nc, 'QTO3',(1./1000.)**2)
sam_wprtp = pull_profiles(nc, 'QTFLUX',(rho * 2.5104*10**6)**-1)

#w
sam_wp2 = pull_profiles(nc, 'W2',1.)
sam_wp3 = pull_profiles(nc, 'W3',1.)

nc.close()


#----------------------------------------------------------------------------------------
# Compute normalized Skewness of w and the normalized correlation of (w and thl) and
# (w and rt). 
sam_Skw = return_skewness(sam_wp3,sam_wp2,w_tol)
sam_Skw_hat = return_nrmlz_skw(sam_Skw,sigma_tilde_w_sqd)
sam_nrmlz_corr_wthl = return_nrmlz_w_corr(sam_wpthlp,sam_wp2,sam_thlp2,
                                        w_tol,thl_tol,sigma_tilde_w_sqd)
sam_nrmlz_corr_wrt = return_nrmlz_w_corr(sam_wprtp,sam_wp2,sam_rtp2,
                                        w_tol,rt_tol,sigma_tilde_w_sqd)
                                        

# Initialize arrays for the prognosed skewness. The first contains the time,height,
# and beta coordinates. The second will be the time mean profile.
prog_Skthl = np.empty( (np.size(sam_Skw,axis=0),
                       np.size(sam_Skw,axis=1),
                       num_beta_lines) )
prog_Skthlm = np.empty( (len(sam_z),num_beta_lines) )

prog_Skrt = np.empty( (np.size(sam_Skw,axis=0),
                       np.size(sam_Skw,axis=1),
                       num_beta_lines) )
prog_Skrtm = np.empty( (len(sam_z),num_beta_lines) )

# Loop through each value of beta.
for i in np.arange(0,num_beta_lines):
    prog_Skthl[:,:,i] = prognose_skex(sam_Skw_hat,sam_nrmlz_corr_wthl,beta[i])
    prog_Skrt[:,:,i] = prognose_skex(sam_Skw_hat,sam_nrmlz_corr_wrt,beta[i])

# Time mean of each 'beta-line'.
for i in np.arange(0,num_beta_lines):
    prog_Skthlm[:,i] = return_mean_profiles(prog_Skthl[:,:,i],
                                          idx_t0, idx_t1, idx_z0, idx_z1)
    prog_Skrtm[:,i] = return_mean_profiles(prog_Skrt[:,:,i],
                                          idx_t0, idx_t1, idx_z0, idx_z1)

#----------------------------------------------------------------------------------------
# Return mean profiles of skewness of w, rt, and thl
#----------------------------------------------------------------------------------------
sam_Skw = (return_mean_profiles(sam_Skw,
           idx_t0, idx_t1, idx_z0, idx_z1))
sam_Skthl = (return_mean_profiles(return_skewness(sam_thlp3,sam_thlp2,thl_tol),
            idx_t0, idx_t1, idx_z0, idx_z1))
sam_Skrt = (return_mean_profiles(return_skewness(sam_rtp3,sam_rtp2,rt_tol),
           idx_t0, idx_t1, idx_z0, idx_z1))


#----------------------------------------------------------------------------------------
# Plot for different values of beta, Skthl
#----------------------------------------------------------------------------------------
print "Plot"
# Create 'throwaway' plot for the colorbar
Z = [[0,0,],[0,0]]
jet = plt.cm.jet
levels = beta
CS3 = plt.contourf(Z,levels,cmap=jet)
plt.clf()

# Thl

idx = np.linspace(0,1,num_beta_lines)
fig, ax = plt.subplots()
ax.set_color_cycle(jet(idx))

for i in range(1, num_beta_lines):
    ax.plot(prog_Skthlm[:,i],sam_z)
ax.plot(sam_Skthl,sam_z,lw=2,c='k',label='Sam')
ax.legend()
ax.set_xlabel('$Sk_{\\theta_{l}}$')
ax.set_ylabel('m')
plt.grid()
plt.colorbar(CS3,orientation='horizontal').set_label('$\\beta$')
plt.savefig(out_dir+out_name+'_prog_Skthl.png')
plt.close()

# rt
fig, ax = plt.subplots()
ax.set_color_cycle(jet(idx))

for i in range(1, num_beta_lines):
    ax.plot(prog_Skrtm[:,i],sam_z)
ax.plot(sam_Skrt,sam_z,lw=2,c='k',label='Sam')
ax.legend()
ax.set_xlabel('$Sk_{r_{t}}$')
ax.set_ylabel('m')
plt.grid()
plt.colorbar(CS3,orientation='horizontal').set_label('$\\beta$')
plt.savefig(out_dir+out_name+'_prog_Skrt.png')

# Plot profiles of SAM's skewness
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
plt.savefig(out_dir+out_name+'_skw_compare.png')
plt.close()
