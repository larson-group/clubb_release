import matplotlib
import matplotlib.pyplot as pyplot
from netCDF4 import Dataset
import numpy as np

matplotlib.rcParams['font.size'] = 16
figsize = (8,5)
figsizethin = (8,3)

# oscillation vs non-oscillation

data = Dataset('rico_godunov-scalarwp3_dt-0.0234375_refine-8_tfinal-3600_zt.nc')
wp3_ta_stable = data['wp3_ta'][:]
z = data['altitude'][:]
data = Dataset('rico_godunov-scalarwp3_dt-0.234375_refine-8_tfinal-3600_zt.nc')
wp3_ta_unstable = data['wp3_ta'][:]

f, ax = pyplot.subplots(figsize=figsize)
ax.plot(wp3_ta_unstable[14,:,0,0], z, ':b', label='t=15min (dt=dt0)')
ax.plot(wp3_ta_stable[14,:,0,0], z, '-b', label='t=15min (dt=dt0/10)', linewidth=2)
ax.plot(wp3_ta_unstable[16,:,0,0], z, ':r', label='t=17min (dt=dt0)')
ax.plot(wp3_ta_stable[16,:,0,0], z, '-r', label='t=17min (dt=dt0/10)', linewidth=2)
ax.plot(wp3_ta_unstable[18,:,0,0], z, ':g', label='t=19min (dt=dt0)')
ax.plot(wp3_ta_stable[18,:,0,0], z, '-g', label='t=19min (dt=dt0/10)', linewidth=2)
ax.legend()
ax.set_xlabel('wp3_ta', fontsize='large')
ax.set_ylabel('altitude (m)', fontsize='large')
ax.set_ylim(550,900)
f.tight_layout()
f.savefig('stable_unstable.png')

# smooth tau and unsmooth tau

data = Dataset('rico_godunov-scalarwp3_dt-0.05_refine-7_tfinal-1200_dt_output-1_zm.nc')
tau_original = data['tau_zm'][:]
z = data['altitude'][:]
data = Dataset('rico_godunov-scalarwp3_dt-0.05_refine-7_bvf-smooth_tfinal-1200_dt_output-1_zm.nc')
tau_smooth = data['tau_zm'][:]

f, ax = pyplot.subplots(2,2,figsize=figsize, sharey=True)
ax[0,0].plot(tau_original[1152,:,0,0], z, label='original BVF')
ax[0,0].plot(tau_smooth[1152,:,0,0], z, label='smoothed BVF')
ax[0,0].set_xlabel('tau at t=1153s', fontsize='large')
ax[0,0].set_ylabel('altitude (m)', fontsize='large')
ax[0,0].set_ylim(1120,1130)
ax[0,1].plot(tau_original[1154,:,0,0], z, label='original BVF')
ax[0,1].plot(tau_smooth[1154,:,0,0], z, label='smoothed BVF')
ax[0,1].set_xlabel('tau at t=1155s', fontsize='large')
ax[1,0].plot(tau_original[1156,:,0,0], z, label='original BVF')
ax[1,0].plot(tau_smooth[1156,:,0,0], z, label='smoothed BVF')
ax[1,0].set_xlabel('tau at t=1157s', fontsize='large')
ax[1,0].set_ylabel('altitude (m)', fontsize='large')
ax[1,0].set_ylim(1120,1130)
ax[1,1].plot(tau_original[1158,:,0,0], z, label='original BVF')
ax[1,1].plot(tau_smooth[1158,:,0,0], z, label='smoothed BVF')
ax[1,1].set_xlabel('tau at t=1159s', fontsize='large')
f.tight_layout()
f.savefig('smooth_unsmooth.png')

# convergence

refineLevels = [0, 1, 2, 3, 4, 5, 6, 8]
dtValuesU = [60, 30, 15, 7.5, 3.75, 1.875, 0.9375, 0.234375]
dtValuesS = [4, 2, 1, 0.5, 0.25, 0.125, 0.0625, 0.015625]
tIndex = 59 

varName = 'thlm'
f, ax = pyplot.subplots(figsize=figsizethin)

for dtValues, label in ((dtValuesU, 'unstable dt/dz'), (dtValuesS, 'stable dt/dz')):

  dataList = []
  for i,refineLevel in enumerate(refineLevels):
    dt = dtValues[i]
    dataList.append(Dataset('rico_godunov-scalarwp3_dt-{}_refine-{}_bvf-smooth_tfinal-3600_zt.nc'.format(dt, refineLevel)))
  
  n = refineLevels[-1] - refineLevels[0]
  varRef = dataList[-1][varName][:,range(0,len(dataList[-1]['altitude'][:]),2**n),0,0]
  
  varList = []
  for n in range(len(dataList)-1):
    data = dataList[n]
    varList.append(data[varName][:,range(0,len(data['altitude'][:]),2**n),0,0])
  
  z = dataList[0]['altitude'][1:]
  
  rmsList = []
  dzList = []
  for n,var in enumerate(varList):
    dzList.append(1.0/2**n)
    rms = np.sqrt(np.mean((var[tIndex,:]-varRef[tIndex,:])**2))
    rmsList.append(rms)
  order = np.log(rmsList[-2]/rmsList[-1])/np.log(2)
  ax.loglog(dzList, rmsList, '-o', label=label)

ax.set_ylabel('RMSE', fontsize='large')
ax.set_xlabel('normalized grid spacing', fontsize='large')
ax.set_ylim(5e-4,8e-2)
ax.legend(loc='best')
f.tight_layout()
f.savefig('convergence.png')

pyplot.show()
