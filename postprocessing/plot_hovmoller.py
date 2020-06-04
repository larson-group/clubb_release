import matplotlib.pyplot as pyplot
import matplotlib.colors as colors
from netCDF4 import Dataset
import numpy as np
import sys

pyplot.rcParams.update({'font.size': 18})

M = len(sys.argv) - 1

if (M < 1):
  print("usage: python plot_hovmoller.py <run 1> [run 2] ...")
  sys.exit(1)

for (suffix, var) in (('zm','wp2'),('zt','wp3')):
  f, ax = pyplot.subplots(M, 1, sharex=True, figsize=(8,4*M))
  if (M == 1):
    ax = [ax,]
  for m in range(M):
    nc_file = sys.argv[1+m] + '_{}.nc'.format(suffix)
    data = Dataset(nc_file)
    t = data['time'][:]
    z = data['altitude'][:]
    val = data[var][:]
    #norm = colors.TwoSlopeNorm(vmin=min(0.0,np.amin(val))-1e-12, vcenter=0.0, vmax=np.amax(val))
    valmax = min(2.0, np.amax(val))
    pcr = ax[m].contourf(t/60.0, z, np.transpose(val[:,:,0,0]), cmap='gist_stern_r', levels=np.linspace(0.0,valmax,10))
    f.colorbar(pcr, ax=ax[m])
    ax[m].set_ylabel(sys.argv[1+m])
    ax[m].grid(True, which='major', axis='x')
  pyplot.setp(ax, xlim=(1000,1440), ylim=(0,5000))
  ax[0].set_title(var)
  ax[M-1].set_xlabel('time (min)')
  f.tight_layout()

pyplot.show()
