import matplotlib.pyplot as pyplot
import numpy as np
from netCDF4 import Dataset
from scipy.interpolate import interp1d
import sys

def odd_extension(f):
  N = len(f)-1
  fe = np.empty(2*N)
  fe[0:N] = -f[range(N,0,-1)]
  fe[N:2*N] = f[1:N+1]
  return fe
 

if (len(sys.argv) < 3):
  sys.exit('usage: python {} <field name> <zt dataset>'.format(sys.argv[0]))

field_name = sys.argv[1]
dataset_name = sys.argv[2]

data = Dataset(dataset_name)
field = data[field_name][:]
time = data['time'][:]
z = data['altitude'][:]
 
f, ax = pyplot.subplots(3, 2, figsize=(8, 9), constrained_layout=True)
ax[0,0].set_title('all wave amplitudes', fontsize='large')
ax[0,1].set_title('grid-scale wave amplitudes', fontsize='large')
ax[2,0].set_xlabel('wavelength (m)', fontsize='large')
ax[2,1].set_xlabel('wavelength (m)', fontsize='large')


max_amp = 0.0
max_time = 0
for n in range(len(time)):

  # create odd extensions
  z_ext = odd_extension(z)
  field_ext = odd_extension(field[n,:,0,0])
  
  # interpolate to uniform grid
  dz = np.min(z[1:]-z[0:-1])
  zi = np.arange(z_ext[0], z_ext[-1], dz)
  interpolant = interp1d(z_ext, field_ext, kind='cubic')
  field_interp = interpolant(zi)
  
  # FFT
  K = len(zi)/2
  field_fft = np.fft.fft(field_interp)/len(field_interp)
  
  # check for max grid-scale amplitude
  k = (z_ext[-1] - z_ext[0])/np.arange(1,K)
  k_grid = 2.0*np.max(z[1:]-z[0:-1])
  indices = np.where(k <= 2.0*k_grid)[0]
  field_fft_gridscale = field_fft[1:K][indices]
  current_max_amp = np.max(np.abs(field_fft_gridscale))

  # plot if appropriate
  i = -1
  if (n == 0):
    i = 0
  elif (current_max_amp > max_amp):
    max_amp = current_max_amp
    max_time = n
    i = 1
    ax[1,0].clear()
    ax[1,1].clear()
  elif (n == len(time)-1):
    i = 2
  
  if (i > -1):
    k = (z_ext[-1] - z_ext[0])/np.arange(1,K)
    ax[i,0].plot(k, np.abs(field_fft[1:K]), ':.')
    ax[i,1].plot(k, np.abs(field_fft[1:K]), ':.')
    ax[i,0].set_ylabel('at time {}'.format(time[n]), fontsize='large')
    ax[i,1].set_xlim(k[-1], k[indices[0]])
    ax[i,1].set_ylim(0.0, current_max_amp)
 
f2, ax2 = pyplot.subplots(constrained_layout=True)
ax2.plot(field[0,:,0,0], z, label='t={}'.format(time[0]))
ax2.plot(field[max_time,:,0,0], z, label='t={}'.format(time[max_time]))
ax2.plot(field[-1,:,0,0], z, label='t={}'.format(time[-1]))
ax2.set_ylim(0,2000)
ax2.legend()
ax2.set_title(dataset_name, fontsize='large')
ax2.set_xlabel(field_name, fontsize='large')
ax2.set_ylabel('altitude (m)', fontsize='large')

f.suptitle(dataset_name)
pyplot.show()

