import matplotlib.pyplot as pyplot
from netCDF4 import Dataset
import numpy as np
import sys

pyplot.rcParams.update({'font.size': 18})

M = len(sys.argv) - 1

if (M < 1):
  print("usage: python plot_hovmoller.py <run 1> [run 2] ...")
  sys.exit(1)

f, ax = pyplot.subplots(M, 1, sharex=True, figsize=(8,4*M))
if (M == 1):
  ax = [ax,]
for m in range(M):
  nc_file = sys.argv[1+m] + '_zm.nc'
  data = Dataset(nc_file)
  t = data['time'][:]
  z = data['altitude'][:]
  wp2 = data['wp2'][:]
  pcr = ax[m].pcolor(t/60.0, z, np.transpose(wp2[:,:,0,0]), cmap='gist_heat_r', vmin=0, vmax=1.0)
  f.colorbar(pcr, ax=ax[m])
  ax[m].set_ylabel(sys.argv[1+m])
  ax[m].grid(True, which='major', axis='x')

pyplot.setp(ax, xlim=(0,1000), ylim=(0,5000))
#pyplot.grid(True, which='major', axis='x')
ax[0].set_title('wp2')
ax[M-1].set_xlabel('time (min)')
f.tight_layout()
pyplot.show()
