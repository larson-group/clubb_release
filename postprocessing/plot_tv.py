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

  if (M < 1):
    print("usage: python plot_tv.py <run 1> [run 2] ...")
    sys.exit(1)

  f, ax = pyplot.subplots()
  for m in range(M):
    nc_file = sys.argv[1+m] + '_zm.nc'
    data = Dataset(nc_file)
    t = data['time'][:]
    wp2 = data['wp2'][:]
    plot_tv(t/60.0, wp2, ax, nc_file)

  ax.set_title('Total Variation of wp2')
  ax.set_xlabel('time (min)')
  ax.set_ylabel('total variation')
  ax.legend(loc='best')
  f.tight_layout()
  pyplot.show()
