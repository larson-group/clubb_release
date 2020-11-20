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

varList = (('wp2',0.0,2.0,'zm'), ('wp3',-0.1,0.1,'zt'))
varList = (('wp3',-0.1,0.1,'zt'), ('wp3_bt',None,None,'zt'), ('wp3_ma',None,None,'zt'), ('wp3_ta',None,None,'zt'), ('wp3_tp',None,None,'zt'), ('wp3_ac',None,None,'zt'), ('wp3_bp1',None,None,'zt'), ('wp3_bp2',None,None,'zt'), ('wp3_pr1',None,None,'zt'), ('wp3_pr2',None,None,'zt'), ('wp3_pr3',None,None,'zt'), ('wp3_dp1',None,None,'zt'), ('wp3_sdmp',None,None,'zt'), ('wp3_cl',None,None,'zt'), ('wp3_splat',None,None,'zt'))

for (var, valmin, valmax, suffix) in varList:
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
    if (valmin is None):
      valmin = np.amin(val)
    if (valmax is None):
      valmax = np.amax(val)
    if (valmax <= valmin):
      continue
    print(var)
    pcr = ax[m].contourf(t/60.0, z, np.transpose(val[:,:,0,0]), cmap='gist_stern_r', levels=np.linspace(valmin,valmax,10))
    f.colorbar(pcr, ax=ax[m])
    ax[m].set_title(sys.argv[1+m])
    ax[m].set_ylabel(var)
    ax[m].grid(True, which='major', axis='x')
  pyplot.setp(ax, ylim=(0,5000))
  ax[M-1].set_xlabel('time (min)')
  f.tight_layout()

pyplot.show()
