import matplotlib.pyplot as pyplot
from netCDF4 import Dataset
import numpy as np
import sys

pyplot.rcParams.update({'font.size': 18})

def plot_tv(t, values, ax, label):

  N = np.shape(values)[0]
  M = np.shape(values)[1]
  tv = np.sum(np.abs(values[:, 1:M, 0, 0] - values[:, 0:(M-1), 0, 0]), axis=1)
  ax.plot(t, tv, label=label, linewidth=2)

if (__name__ == '__main__'):

  M = len(sys.argv) - 1

  if (M < 2):
    print("usage: python plot_tv.py <variable> <run 1> [run 2] ...")
    sys.exit(1)

  var_name = sys.argv[1]
  f, ax = pyplot.subplots()
  for m in range(2,M+1):
    nc_file = sys.argv[m] + '_zm.nc'
    data = Dataset(nc_file)
    if (var_name not in data.variables):
      nc_file = sys.argv[m] + '_zt.nc'
      data = Dataset(nc_file)
    t = data['time'][:]
    var = data[var_name][:]
    plot_tv(t/60.0, var, ax, nc_file)

  ax.set_title('Total Variation of {}'.format(var_name))
  ax.set_xlabel('time (min)')
  ax.set_ylabel('total variation')
  ax.legend(loc='best')
  f.tight_layout()
  pyplot.show()
