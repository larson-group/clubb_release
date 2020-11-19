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
  # create (thlm(t+dt)-thlm(t))/dt & (thlmwp(t+dt)-thlmw(t))/dt
  dthlmdt = (zt_data['thlm'][1:,:,:,:]-zt_data['thlm'][0:-1,:,:,:])/dt
  dwpthlpdt = (zm_data['wpthlp'][1:,:,:,:]-zm_data['wpthlp'][0:-1,:,:,:])/dt

  # ***** thlm *****
  cax = ax[0,0]
  add_to_axis(zt_data['thlm'][:]-300.0, zt, 'thlm-300', cax, data)
  add_to_axis(zt_data['cloud_frac'][:], zt, 'cld_frac', cax, data)
  add_to_axis(dthlmdt, zt, 'dthlm/dt', cax, data, linespec='--k')
  cax.set_xlim(-5,7)
  # significant quantities
  cax = ax[0,1]
  add_to_axis(zt_data['thlm_ta'][:], zt, 'ta', cax, data)
  add_to_axis(dthlmdt, zt, 'dthlm/dt', cax, data, linespec='--k')
  cax.set_xlim(-0.015,0.015)
  # insignificant_quantities
  cax = ax[0,2]
  add_to_axis(zt_data['thlm_ma'][:], zt, 'ma', cax, data)
  add_to_axis(zt_data['thlm_forcing'][:], zt, 'forcing', cax, data)
  cax.set_xlim(-0.0001,0.0001)
  # numerical adjustments
  cax = ax[0,3]
  add_to_axis(zt_data['thlm_sdmp'][:], zt, 'sdmp', cax, data)
  add_to_axis(zt_data['thlm_cl'][:], zt, 'cl', cax, data)
  add_to_axis(zt_data['thlm_mfl'][:], zt, 'mfl', cax, data)
  add_to_axis(zt_data['thlm_tacl'][:], zt, 'tacl', cax, data)
  cax.set_xlim(-0.0001,0.0001)

  # ***** wpthlp *****
  cax = ax[1,0]
  add_to_axis(zm_data['wpthlp'][:], zm, 'wpthlp', cax, data)
#  #add_to_axis(zt_data['cloud_frac'][:], zt, 'cld_frac', cax, data)
  add_to_axis(dwpthlpdt, zm, 'dwpthlp/dt', cax, data, linespec='--k')
  cax.set_xlim(-2.0,2.0)
  # significant quantities
  cax = ax[1,1]
  add_to_axis(zm_data['wpthlp_tp'][:], zm, 'tp', cax, data)
  add_to_axis(zm_data['wpthlp_pr1'][:], zm, 'pr1', cax, data)
  add_to_axis(zm_data['wpthlp_pr3'][:], zm, 'pr3', cax, data)
  add_to_axis(zm_data['wpthlp_bp'][:], zm, 'bp', cax, data)
  add_to_axis(zm_data['wpthlp_dp1'][:], zm, 'dp1', cax, data)
  add_to_axis(dwpthlpdt, zm, 'dwpthlp/dt', cax, data, linespec='--k')
  cax.set_xlim(-0.05,0.05)
  # insignificant_quantities
  cax = ax[1,2]
  add_to_axis(zm_data['wpthlp_ma'][:], zm, 'ma', cax, data)
  add_to_axis(zm_data['wpthlp_ta'][:], zm, 'ta', cax, data)
  add_to_axis(zm_data['wpthlp_ac'][:], zm, 'ac', cax, data)
  add_to_axis(zm_data['wpthlp_pr2'][:], zm, 'pr2', cax, data)
  add_to_axis(zm_data['wpthlp_forcing'][:], zm, 'forcing', cax, data)
#  add_to_axis(dwp3dt, zt, 'dwp3/dt', cax, data, linespec='--k')
  cax.set_xlim(-0.01,0.01)
  # numerical adjustments
  cax = ax[1,3]
  add_to_axis(zm_data['wpthlp_mc'][:], zm, 'mc', cax, data)
#  add_to_axis(zt_data['wp3_cl'][:], zt, 'cl', cax, data)
#  add_to_axis(zt_data['wp3_splat'][:], zt, 'splat', cax, data)
#  add_to_axis(dwp3dt, zt, 'dwp3/dt', cax, data, linespec='--k')
#  cax.set_xlim(-0.05,0.05)
#
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
    ax[0,m].set_ylim(0.0,2000.0)
    ax[1,m].set_ylim(0.0,2000.0)
    t = time[n0]
    #ax[1,m].set_xlabel('t = %d min of %d min' % (t/60.0, t_final/60.0), fontsize='large')
    ax[1,m].set_xlabel('t = %d s of %d s' % (t, t_final), fontsize='large')

  f.tight_layout()

  def update(i):
    t = time[i]
    for m in range(M):
      #ax[1,m].set_xlabel('t = %d min of %d min' % (t/60.0, t_final/60.0), fontsize='large')
      ax[1,m].set_xlabel('t = %d s of %d s' % (t, t_final), fontsize='large')
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
