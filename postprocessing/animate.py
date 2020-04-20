import matplotlib.animation as animation
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import sys
from netCDF4 import Dataset

matplotlib.rcParams.update({'font.size': 16})


if (__name__ == '__main__'):

  timestep = 60.0

  if (len(sys.argv) < 3):
    print("usage: python animate.py zm_file zt_file <start_index> <stop_index> <output_file>")
    sys.exit(1) 

  zm_file = sys.argv[1]
  zm_data = Dataset(zm_file)
  zt_file = sys.argv[2]
  zt_data = Dataset(zt_file)

  if (len(sys.argv) > 3):
    n0 = int(sys.argv[3])
  else:
    n0 = 0

  if (len(sys.argv) > 4):
    nMax = int(sys.argv[4])
  else:
    nMax = -1

  if (len(sys.argv) > 5):
    outputFile = sys.argv[5]
  else:
    outputFile = None

  def add_to_axis(datum, z, name, ax, dict, linespec='-'):
    line, = ax.plot(datum[n0,:,0,0], z, linespec, label=name)
    dict[line] = datum

  zm = zm_data['altitude'][:]
  zt = zt_data['altitude'][:]

  M = 3
  f, ax = plt.subplots(1,M, figsize=(5*M,6),sharey=True)
  data = {}
  # create (wp2(t+dt)-wp2(t))/dt
  dwp2dt = (zm_data['wp2'][1:,:,:,:]-zm_data['wp2'][0:-1,:,:,:])/timestep
  # wp2
  add_to_axis(zm_data['wp2'][:], zm, 'wp2', ax[0], data)
  add_to_axis(zt_data['cloud_frac'][:], zt, 'cld_frac', ax[0], data)
  add_to_axis(dwp2dt, zm, 'dwp2/dt', ax[0], data, linespec='--k')
  ax[0].set_xlim(-0.005,1.0)
  # significant quantities
  add_to_axis(zm_data['wp2_ta'][:], zm, 'ta', ax[1], data)
  add_to_axis(zm_data['wp2_pr1'][:], zm, 'pr1', ax[1], data)
  add_to_axis(zm_data['wp2_pr2'][:], zm, 'pr2', ax[1], data)
  add_to_axis(zm_data['wp2_bp'][:], zm, 'bp', ax[1], data)
  add_to_axis(dwp2dt, zm, 'dwp2/dt', ax[1], data, linespec='--k')
  ax[1].set_xlim(-0.005,0.005)
  # insignificant_quantities
  add_to_axis(zm_data['wp2_ma'][:], zm, 'ma', ax[2], data)
  add_to_axis(zm_data['wp2_ac'][:], zm, 'ac', ax[2], data)
  add_to_axis(zm_data['wp2_pr3'][:], zm, 'pr3', ax[2], data)
  add_to_axis(zm_data['wp2_dp1'][:], zm, 'dp1', ax[2], data)
  add_to_axis(zm_data['wp2_dp2'][:], zm, 'dp2', ax[2], data)
  add_to_axis(dwp2dt, zm, 'dwp2/dt', ax[2], data, linespec='--k')
  ax[2].set_xlim(-0.005,0.005)
  # numerical adjustments
  if 0:
    add_to_axis(zm_data['wp2_sdmp'][:], zm, 'sdmp', ax[3], data)
    add_to_axis(zm_data['wp2_cl'][:], zm, 'cl', ax[3], data)
    add_to_axis(zm_data['wp2_pd'][:], zm, 'pd', ax[3], data)
    add_to_axis(dwp2dt, zm, 'dwp2/dt', ax[3], data, linespec='--k')
    ax[3].set_xlim(-0.005,0.005)

  time = zm_data['time'][:]
  if (nMax > 0):
    N = nMax-1
  else:
    N = len(time)-1
  t_final = time[N]

  ax[0].set_ylabel('altitude (m)', fontsize='large')
  for axis in ax:
    axis.legend(fontsize='small')
    axis.set_ylim(0.0,5000.0)
    t = time[n0]
    axis.set_xlabel('t = %d min of %d min' % (t/60.0, t_final/60.0), fontsize='large')

  f.tight_layout()

  def update(i):
    t = time[i]
    for axis in ax:
      axis.set_xlabel('t = %d min of %d min' % (t/60.0, t_final/60.0), fontsize='large')
    n = 0
    while (time[n+1] < t):
      n = n+1

      for line in data:
        line.set_xdata(data[line][n,:,0,0])

  ani = animation.FuncAnimation(f, update, frames=range(n0,N), interval=1)
  if (outputFile is None):
    plt.show()
  else:
    ani.save(outputFile, fps=5)
