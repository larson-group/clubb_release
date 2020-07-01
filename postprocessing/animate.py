import matplotlib.animation as animation
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import sys
from netCDF4 import Dataset

matplotlib.rcParams.update({'font.size': 16})


if (__name__ == '__main__'):

  if (len(sys.argv) < 2):
    print("usage: python animate.py run_name  <start_index> <stop_index> <animation_step> <output_file>")
    sys.exit(1)

  zm_file = sys.argv[1] + '_zm.nc'
  zm_data = Dataset(zm_file)
  zt_file = sys.argv[1] + '_zt.nc'
  zt_data = Dataset(zt_file)

  if (len(sys.argv) > 2):
    n0 = int(sys.argv[2])
  else:
    n0 = 0

  if (len(sys.argv) > 3):
    nMax = int(sys.argv[3])
  else:
    nMax = -1

  if (len(sys.argv) > 4):
    step = int(sys.argv[4])
  else:
    step = 1

  if (len(sys.argv) > 5):
    outputFile = sys.argv[5]
  else:
    outputFile = None

  def add_to_axis(datum, z, name, ax, dict, linespec='-'):
    line, = ax.plot(datum[n0,:,0,0], z, linespec, label=name)
    dict[line] = datum

  zm = zm_data['altitude'][:]
  zt = zt_data['altitude'][:]
  time = zm_data['time'][:]
  dt = time[1]-time[0]

  M = 4
  f, ax = plt.subplots(2, M, figsize=(5*M,12),sharey=True)
  data = {}
  # create (wp2(t+dt)-wp2(t))/dt & (wp3(t+dt)-wp3(t))/dt
  dwp2dt = (zm_data['wp2'][1:,:,:,:]-zm_data['wp2'][0:-1,:,:,:])/dt
  dwp3dt = (zt_data['wp3'][1:,:,:,:]-zt_data['wp3'][0:-1,:,:,:])/dt

  # ***** wp2 *****
  cax = ax[0,0]
  add_to_axis(zm_data['wp2'][:], zm, 'wp2', cax, data)
  add_to_axis(zt_data['cloud_frac'][:], zt, 'cld_frac', cax, data)
  add_to_axis(dwp2dt, zm, 'dwp2/dt', cax, data, linespec='--k')
  cax.set_xlim(-0.005,5.0)
  # significant quantities
  cax = ax[0,1]
  add_to_axis(zm_data['wp2_ta'][:], zm, 'ta', cax, data)
  add_to_axis(zm_data['wp2_pr1'][:], zm, 'pr1', cax, data)
  add_to_axis(zm_data['wp2_pr2'][:], zm, 'pr2', cax, data)
  add_to_axis(zm_data['wp2_bp'][:], zm, 'bp', cax, data)
  add_to_axis(dwp2dt, zm, 'dwp2/dt', cax, data, linespec='--k')
  cax.set_xlim(-0.05,0.05)
  # insignificant_quantities
  cax = ax[0,2]
  add_to_axis(zm_data['wp2_ma'][:], zm, 'ma', cax, data)
  add_to_axis(zm_data['wp2_ac'][:], zm, 'ac', cax, data)
  add_to_axis(zm_data['wp2_pr3'][:], zm, 'pr3', cax, data)
  add_to_axis(zm_data['wp2_dp1'][:], zm, 'dp1', cax, data)
  add_to_axis(zm_data['wp2_dp2'][:], zm, 'dp2', cax, data)
  add_to_axis(dwp2dt, zm, 'dwp2/dt', cax, data, linespec='--k')
  cax.set_xlim(-0.05,0.05)
  # numerical adjustments
  cax = ax[0,3]
  add_to_axis(zm_data['wp2_sdmp'][:], zm, 'sdmp', cax, data)
  add_to_axis(zm_data['wp2_cl'][:], zm, 'cl', cax, data)
  add_to_axis(zm_data['wp2_pd'][:], zm, 'pd', cax, data)
  add_to_axis(dwp2dt, zm, 'dwp2/dt', cax, data, linespec='--k')
  cax.set_xlim(-0.05,0.05)

  # ***** wp3 *****
  cax = ax[1,0]
  add_to_axis(zt_data['wp3'][:], zt, 'wp3', cax, data)
  add_to_axis(zt_data['wp3_on_wp2_zt'][:]/10.0, zt, '0.1*wp3/wp2 (zt)', cax, data)
  add_to_axis(zm_data['wp3_on_wp2'][:]/10.0, zm, '0.1*wp3/wp2 (zm)', cax, data)
  #add_to_axis(zt_data['cloud_frac'][:], zt, 'cld_frac', cax, data)
  add_to_axis(dwp3dt, zt, 'dwp3/dt', cax, data, linespec='--k')
  cax.set_xlim(-0.1,0.1)
  # significant quantities
  cax = ax[1,1]
  add_to_axis(zt_data['wp3_ta'][:], zt, 'ta', cax, data)
  add_to_axis(zt_data['wp3_pr1'][:], zt, 'pr1', cax, data)
  add_to_axis(zt_data['wp3_pr2'][:], zt, 'pr2', cax, data)
  add_to_axis(zt_data['wp3_bp1'][:], zt, 'bp1', cax, data)
  add_to_axis(zt_data['wp3_tp'][:], zt, 'tp', cax, data)
  add_to_axis(dwp3dt, zt, 'dwp3/dt', cax, data, linespec='--k')
  cax.set_xlim(-0.05,0.05)
  # insignificant_quantities
  cax = ax[1,2]
  add_to_axis(zt_data['wp3_ma'][:], zt, 'ma', cax, data)
  add_to_axis(zt_data['wp3_ac'][:], zt, 'ac', cax, data)
  add_to_axis(zt_data['wp3_pr3'][:], zt, 'pr3', cax, data)
  add_to_axis(zt_data['wp3_dp1'][:], zt, 'dp1', cax, data)
  add_to_axis(zt_data['wp3_bp2'][:], zt, 'bp2', cax, data)
  add_to_axis(dwp3dt, zt, 'dwp3/dt', cax, data, linespec='--k')
  cax.set_xlim(-0.05,0.05)
  # numerical adjustments
  cax = ax[1,3]
  add_to_axis(zt_data['wp3_sdmp'][:], zt, 'sdmp', cax, data)
  add_to_axis(zt_data['wp3_cl'][:], zt, 'cl', cax, data)
  add_to_axis(zt_data['wp3_splat'][:], zt, 'splat', cax, data)
  add_to_axis(dwp3dt, zt, 'dwp3/dt', cax, data, linespec='--k')
  cax.set_xlim(-0.05,0.05)

  if (nMax > 0):
    N = nMax-1
  else:
    N = len(time)-1
  t_final = time[N]

  ax[0,0].set_ylabel('altitude (m)', fontsize='large')
  ax[1,0].set_ylabel('altitude (m)', fontsize='large')
  for m in range(M):
    ax[0,m].legend(fontsize='small')
    ax[1,m].legend(fontsize='small')
    ax[0,m].set_ylim(0.0,5000.0)
    ax[1,m].set_ylim(0.0,5000.0)
    t = time[n0]
    ax[1,m].set_xlabel('t = %d min of %d min' % (t/60.0, t_final/60.0), fontsize='large')

  f.tight_layout()

  def update(i):
    t = time[i]
    for m in range(M):
      ax[1,m].set_xlabel('t = %d min of %d min' % (t/60.0, t_final/60.0), fontsize='large')
    n = 0
    while (time[n+1] < t):
      n = n+1

      for line in data:
        line.set_xdata(data[line][n,:,0,0])

  ani = animation.FuncAnimation(f, update, frames=range(n0,N,step), interval=1)
  if (outputFile is None):
    plt.show()
  else:
    ani.save(outputFile, fps=5)
