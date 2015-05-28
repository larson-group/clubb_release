# $Id$

#----------------------------------------------------------------------------------------
"""
 compare_CLUBB_SAM_Morr_warm_rain_processes.py
 
 Description:
   Plots and compares SAM's and CLUBB's Morrison microphysics warm rain processes
"""

# Import libraries
import numpy as np
import matplotlib.pyplot as plt
import netCDF4

# Files in which the tendencies are located
sam_file = '/home/weberjk/Summarize_Progress/CLUBB/input/input_fields/LBA.nc'
clubb_file = '/home/weberjk/Summarize_Progress/CLUBB/output/prescribe_thl_thlp2_rt_rtp2_wp2_wp3/lba_zt.nc'

# Specify a directory and name for output
savedir = '/home/weberjk/Summarize_Progress/'
save_filename = 'LBA_warmrain.png'

z0 = 0    # Height start [m]
z1 = 6000 # Height end

t0 = 189  # Averaging interval start [min]
t1 = 249  # Averaging interval end


#----------------------------------------------------------------------------------------
# Should not have to edit below this line. 
#----------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------
# Useful constants 
#----------------------------------------------------------------------------------------
s_in_min = 60

#----------------------------------------------------------------------------------------
# Functions 
#----------------------------------------------------------------------------------------
def pull_profiles(nc, varname):
    """
    Input:
      nc         --  Netcdf file object
      varname    --  Variable name string

    Output:
      time x height array of the specified variable
    """

    var = nc.variables[varname]
    var = np.squeeze(var)
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
# Begin code 
#----------------------------------------------------------------------------------------

t0_in_s = t0*s_in_min # CLUBB's time is in seconds.
t1_in_s = t1*s_in_min

print "Grab SAM's profiles"
nc = netCDF4.Dataset(sam_file)
sam_time = pull_profiles(nc, 'time')
sam_z = pull_profiles(nc, 'z')

# Find the indices closest to z0, and z1
idx_z0 = (np.abs(sam_z[:] - z0)).argmin()
idx_z1 = (np.abs(sam_z[:] - z1)).argmin()
sam_z = sam_z[idx_z0:idx_z1]

# Find the indices closest to t0, and t1
idx_t0 = (np.abs(sam_time[:] - t0)).argmin()
idx_t1 = (np.abs(sam_time[:] - t1)).argmin()

# SAM's
# Warm-rain processes for rrm. Evaporation (PRE), Accretion (PRA), and Autoconversion (PRC)
sam_pre = return_mean_profiles(pull_profiles(nc, 'PRE'),idx_t0,idx_t1,idx_z0,idx_z1)
sam_pra = return_mean_profiles(pull_profiles(nc, 'PRA'),idx_t0,idx_t1,idx_z0,idx_z1)
sam_prc = return_mean_profiles(pull_profiles(nc, 'PRC'),idx_t0,idx_t1,idx_z0,idx_z1)

# Warm-rain processes for Nrm. Autoconversion (NPRC1), Self-Collection (NRAGG), 
# and Evaporation(NSUBR).
sam_nprc = return_mean_profiles(pull_profiles(nc, 'NPRC1'),idx_t0,idx_t1,idx_z0,idx_z1)
sam_nragg = return_mean_profiles(pull_profiles(nc, 'NRAGG'),idx_t0,idx_t1,idx_z0,idx_z1)
sam_nsubr = return_mean_profiles(pull_profiles(nc, 'NSUBR'),idx_t0,idx_t1,idx_z0,idx_z1)

print "Grab CLUBB's profiles"
nc = netCDF4.Dataset(clubb_file)
clb_time = pull_profiles(nc, 'time')
clb_z = pull_profiles(nc, 'altitude')

idx_z0 = (np.abs(clb_z[:] - z0)).argmin()
idx_z1 = (np.abs(clb_z[:] - z1)).argmin()
clb_z = clb_z[idx_z0:idx_z1]

idx_t0 = (np.abs(clb_time[:] - t0_in_s)).argmin()
idx_t1 = (np.abs(clb_time[:] - t1_in_s)).argmin()

# Rain mass
clb_pre = return_mean_profiles(pull_profiles(nc, 'PRE'),idx_t0,idx_t1,idx_z0,idx_z1)
clb_pra = return_mean_profiles(pull_profiles(nc, 'PRA'),idx_t0,idx_t1,idx_z0,idx_z1)
clb_prc = return_mean_profiles(pull_profiles(nc, 'PRC'),idx_t0,idx_t1,idx_z0,idx_z1)

# Rain number
clb_nprc = return_mean_profiles(pull_profiles(nc, 'NPRC1'),idx_t0,idx_t1,idx_z0,idx_z1)
clb_nragg = return_mean_profiles(pull_profiles(nc, 'NRAGG'),idx_t0,idx_t1,idx_z0,idx_z1)
clb_nsubr = return_mean_profiles(pull_profiles(nc, 'NSUBR'),idx_t0,idx_t1,idx_z0,idx_z1)

sam_rr_sum = sam_pre[:] + sam_pra[:] + sam_prc[:]
sam_Nr_sum = sam_nprc[:] + sam_nragg[:] + sam_nsubr[:]

clb_rr_sum = clb_pre[:] + clb_pra[:] + clb_prc[:]
clb_Nr_sum = clb_nprc[:] + clb_nragg[:] + clb_nsubr[:]


print "Plot"
f, (ax1, ax2) = plt.subplots(1,2,sharey=True)
ax1.set_ylim([z0,z1])
#SAM
ax1.plot(sam_prc,sam_z, c='green', linewidth=2, ls='--')
ax1.plot(sam_pra,sam_z, c='blue', linewidth=2, ls='--')
ax1.plot(sam_pre,sam_z, c='red', linewidth=2, ls='--')
ax1.plot(sam_rr_sum,sam_z, c='k', linewidth=2, ls='--')

#CLUBB
ax1.plot(clb_prc,clb_z, c='green', linewidth=2,label='Auto. Conv.' )
ax1.plot(clb_pra,clb_z, c='blue', linewidth=2 ,label='Accr.')
ax1.plot(clb_pre,clb_z,c='red', linewidth=2 ,label='Evap.')
ax1.plot(clb_rr_sum,clb_z,c='k', linewidth=2 ,label='CLUBB Sum')
ax1.legend()
#SAM
ax2.plot(sam_nprc,sam_z, c='green', linewidth=2, ls='--')
ax2.plot(sam_nragg,sam_z, c='purple', linewidth=2, ls='--')
ax2.plot(sam_nsubr,sam_z,c='red', linewidth=2, ls='--')
ax2.plot(sam_Nr_sum,sam_z, c='k', linewidth=2, ls='--')

#CLUBB
ax2.plot(clb_nprc,clb_z, c='green', linewidth=2, label = 'Auto. Conv.')
ax2.plot(clb_nragg,clb_z, c='purple', linewidth=2, label = 'Self. Coll.')
ax2.plot(clb_nsubr,clb_z,c='red', linewidth=2, label = 'Evap.')
ax2.plot(clb_Nr_sum,clb_z,c='k', linewidth=2 ,label='CLUBB Sum')
ax2.legend()

ax1.grid()
ax2.grid()

ax1.set_title('$\overline{r_{r}}$')
ax2.set_title('$\overline{N_{r}}$')

ax1.set_xlabel(r'$kg\ kg^{-1}\ s^{-1}$')
ax2.set_xlabel(r'$\#\ kg^{-1}\ s^{-1}$')
ax1.set_ylabel('$m$')
f.suptitle('Warm-rain budgets averaged over t=%d to %d min.'%(t0,t1))
f.savefig(savedir+save_filename)
plt.close()
