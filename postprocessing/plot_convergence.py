import matplotlib.pyplot as pyplot
import numpy as np
from netCDF4 import Dataset

refineLevels = [0, 1, 2, 3, 4, 5, 6, 8]
#dtValues = [60, 30, 15, 7.5, 3.75, 1.875, 0.9375, 0.234375]
#idtValues = [10, 5, 2.5, 1.25, 0.625, 0.3125, 0.15625, 0.0390625]
dtValues = [6, 3, 1.5, 0.75, 0.375, 0.1875, 0.09375, 0.0234375]
#dtValues = [4, 2, 1, 0.5, 0.25, 0.125, 0.0625, 0.015625]
tIndices = [4, 9, 14, 19, 24, 29, 34, 39, 44, 49, 54, 59] 

varName = 'thlm'
suffix = 'zt'

if (suffix == 'zm'):
  def startIndex(n):
    return 0
else:
  def startIndex(n):
    return 2**n

dataList = []
for i,refineLevel in enumerate(refineLevels):
  dt = dtValues[i]
#  dataList.append(Dataset('rico_godunov-scalarwp3_dt-{}_refine-{}_tfinal-3600_{}.nc'.format(dt, refineLevel, suffix)))
  dataList.append(Dataset('rico_dt-{}_refine-{}_tfinal-3600_{}.nc'.format(dt, refineLevel, suffix)))

n = refineLevels[-1] - refineLevels[0]
varRef = dataList[-1][varName][:,range(startIndex(n),len(dataList[-1]['altitude'][:]),2**n),0,0]

varList = []
for n in range(len(dataList)-1):
  data = dataList[n]
  varList.append(data[varName][:,range(startIndex(n),len(data['altitude'][:]),2**n),0,0])

z = dataList[0]['altitude'][1:]
f, ax = pyplot.subplots()
for tIndex in tIndices:
  rmsList = []
  dzList = []
  for n,var in enumerate(varList):
    dzList.append(1.0/2**n)
    rms = np.sqrt(np.mean((var[tIndex,:]-varRef[tIndex,:])**2))
    rmsList.append(rms)
  order = np.log(rmsList[-2]/rmsList[-1])/np.log(2)
  ax.loglog(dzList, rmsList, '-o', label='t={}s ({:.2f})'.format(dataList[0]['time'][tIndex], order))

ax.set_title('Convergence of {} at with $dt_0$={}'.format(varName,dtValues[0]), fontsize='x-large')
ax.set_ylabel('RMSE', fontsize='x-large')
ax.set_xlabel('normalized grid spacing', fontsize='x-large')
ax.set_ylim(4e-5,5e-2)
ax.axis('equal')
ax.legend(loc='best',fontsize='large')

f.tight_layout()
pyplot.show()

