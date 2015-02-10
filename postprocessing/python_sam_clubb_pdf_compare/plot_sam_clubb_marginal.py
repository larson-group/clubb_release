# $Id$

#----------------------------------------------------------------------------------------
# compare_CLUBB_SAM_marginal.py
# 
# Description:
#   Plots and compares SAM's and CLUBB's marginal for some variable 

import numpy as np
import matplotlib.pyplot as plt
import netCDF4

#----------------------------------------------------------------------------------------    
# User input
#----------------------------------------------------------------------------------------

# Directory for SAM's 3D files 
sam_3D_dir = '/home/weberjk/Plot_CRM/SAM_LBA/Data/3D/1km/'

# Note: SAM outputs 3D files in NSTEPS, so it is dependent on timestep. Therefore, 
#   this script must do some gymnastics to ensure we are comparing the correct times.
#   Only supply 'a' filename upto the NSTEP counter.  
sam_3D_file = 'LBA_128kmx1kmx128_1km_Morrison_64_'
# e.g.        LBA_128kmx1kmx128_1km_Morrison_64_0000000150_micro.nc

sam_dt = 6. # [s] SAM's model timestep, not output frequency. 

# SAM's stat file 
sam_stat_file = '/home/weberjk/Plot_CRM/SAM_LBA/Data/LBA.nc'

# Directory for CLUBB's sample points file
clubb_3D_dir = '/home/weberjk/Summarize_Progress/CLUBB/output/prescribe_thl_thlp2_rt_rtp2/'

# CLUBB's file names
clb_nl_file = clubb_3D_dir + 'lba_nl_lh_sample_points_2D.nc'
clb_u_file = clubb_3D_dir + 'lba_u_lh_sample_points_2D.nc'
clubb_stat_file = clubb_3D_dir + 'lba_zt.nc' # The statistics file and sample points
                                           # files are output from CLUBB into the 
                                           # same directory by default
                                           

# Variable name, time, and height
SAM_hist_var_name = 'CHI'   # SAM's name for variable
CLUBB_hist_var_name = 'chi' # CLUBB's name for variable

out_dir = '/home/weberjk/Summarize_Progress/CLUBB/output/prescribe_thl_thlp2_rt_rtp2/' # Figure output directory
out_name = 'hist_chi.png'                                # and and name

plot_var_name = '$\chi$'         # How the variable looks when it is plotted
plot_var_units = '$kg\ kg^{-1}$' # The plotted units
                      
z0 = 1800                  # Analysis height [m]
t0_in_min = 225            # Analysis time [min] 
                            # Note: This exact time must exist for SAM's output. 
                            # CLUBB will use the closest available.

#----------------------------------------------------------------------------------------
# Should not have to edit below this line
#----------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------
# Definitions
#----------------------------------------------------------------------------------------
def pull_samples(nc, varname):
    var = nc.variables[varname]
    var_u = nc.variables[varname].units
    var = np.squeeze(var)
    return var, var_u

def reshape_sam_slice(var, idx_z0):
    var = np.squeeze(var)
    var = np.swapaxes(var,2,0)
    var = var[:,:,idx_z0]
    var = np.reshape(var,(np.size(var,0)*np.size(var,1)))
    return var
    
def reshape_clubb_samples(var, idx_z0, idx_t0):
    var = var[idx_t0,idx_z0,:]
    return var


# Find the correct SAM 3D file.
t0_in_s = t0_in_min*60.

# Find which NSTEP corrosponds to our desired t0
t0_filename = str(int((t0_in_s)/sam_dt))

# SAM's nstep counter is a fixed 10 digits long.
sam_nstep_counter='0000000000'

# Find how many digits the desire NSTEP takes up
place = len(t0_filename)

# From the beginning of the nstep counter, find how many places until we insert our 
# NSTEP
end_zeros = len(sam_nstep_counter) - place

# Combine SAM's nstep counter and our desired NSTEP
t0_filename = sam_nstep_counter[0:end_zeros] + t0_filename

# Paste the strings together and call it a filename.
sam_3D_file = sam_3D_dir+sam_3D_file+t0_filename+'_micro.nc'


print "Extract SAM Slice"
nc = netCDF4.Dataset(sam_3D_file)
sam_z, z_u = pull_samples(nc, 'z')
sam_chi, chi_u = pull_samples(nc, SAM_hist_var_name)
nc.close()    

print "Find Height Indicies for SAM Slice"
idx_z0 = (np.abs(sam_z[:] - z0)).argmin()
sam_chi = reshape_sam_slice(sam_chi,  idx_z0)

print "Extract mixt_frac"
nc = netCDF4.Dataset(clubb_stat_file)
time, time_u = pull_samples(nc, 'time')
clb_z, z_u = pull_samples(nc, 'altitude')
mixt_frac, mixt_frac_u = pull_samples(nc, 'mixt_frac')
nc.close()

print "Find Time and Height Indicies for CLUBB's mixt_frac"
idx_t0 = (np.abs(time[:] - t0_in_s)).argmin()
idx_z0 = (np.abs(clb_z[:] - z0)).argmin()
mixt_frac = mixt_frac[idx_t0,idx_z0] 

print "Extract CLUBB Data"
nc = netCDF4.Dataset(clb_nl_file)
clb_z, z_u = pull_samples(nc, 'altitude')
clb_time, clb_time_u = pull_samples(nc, 'time')
clb_chi, chi_u = pull_samples(nc, CLUBB_hist_var_name)
nc.close()

nc = netCDF4.Dataset(clb_u_file)
dp1, dp1_u = pull_samples(nc, 'dp1')
nc.close()

idx_t0 = (np.abs(clb_time[:] - t0_in_s)).argmin()
idx_z0 = (np.abs(clb_z[:] - z0)).argmin()

clb_chi = reshape_clubb_samples(clb_chi, idx_z0,idx_t0)
dp1 = reshape_clubb_samples(dp1, idx_z0,idx_t0)

# Seperate clb_chi into components 1 and 2. Where dp1 < mixt_frac, then the sample 
# is from component 1. 
idx_comp_1 = np.where(dp1 < mixt_frac)
idx_comp_2 = np.where(dp1 >= mixt_frac)
clb_chi_1 = clb_chi[idx_comp_1]
clb_chi_2 = clb_chi[idx_comp_2]

max_chi = max(max(clb_chi),max(sam_chi))
min_chi = min(min(clb_chi),min(sam_chi))
H, sam_bins = np.histogram(sam_chi, bins = 100, range=(min_chi,max_chi), density=True)

# Plot histograms
f, ax1 = plt.subplots()
plt.hist([clb_chi_1,clb_chi_2], bins=sam_bins, stacked='True', color = ('red','blue'), normed=True, label=('Comp.-1','Comp.-2'))
plt.hist(sam_chi, bins=sam_bins, normed=True, color = 'k', lw = 3, histtype='step', label='SAM')
plt.grid()
plt.legend()
ax1.set_xlabel(plot_var_units)
ax1.set_ylabel('Density')
ax1.set_title('Marginal of %s at z=%d m and t=%d min.'%(plot_var_name, z0, t0_in_min))
plt.savefig(out_dir + out_name)
plt.close()
