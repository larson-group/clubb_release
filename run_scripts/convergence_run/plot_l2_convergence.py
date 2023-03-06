import numpy as np
from cycler import cycler
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
from matplotlib.font_manager import FontProperties
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from matplotlib import rc
import seaborn as sns

case         = "CASE_NAME"
expnam       = "EXP_NAME"
#refinment factor and time step sizes 
#refineLevels = [0, 1, 2, 3, 4, 5, 7]
refineLevels = [REFINE_LEVELS]
#coresponded time steps
#dtValues     = [4, 2, 1, 0.5, 0.25, 0.125, 0.03125]
dtValues     = [DTIME_VALUES]
#corresponded plotting time index for 1-h, 2-h, 3-h, 4-h, 5-h, 6-h 
#tIndices     = [5, 11, 17, 23, 29, 35]
tIndices     = [DTIME_INDEX]

Datadir      = '../output'
Filename     = Datadir + '/' + case + '_dt-{}_refine-{}_output_name-'+expnam+'_{}.nc'
Figstr       = case+'_l2_cnvg_'+expnam
tstart       = 0

color_cycle  = ['#377eb8', '#ff7f00', '#4daf4a',
                '#f781bf', '#a65628', '#984ea3',
                '#999999', '#e41a1c', '#dede00',
                '#76b7b2', '#121211', '#FF00FF']
rc('font',**{'family':'serif','serif':['Times New Roman'],'style':'normal','weight':'normal'})
#rc('text',usetex=True)
font_size = 20

##variable names and output levels 
vnmList      = ['wp3','wp2','rtm','thlm','rtp2','thlp2','rtpthlp', 'um', 'vm', 'upwp', 'vpwp', 'wprtp', 'wpthlp', 'up2', 'vp2']
sufList      = [ 'zt', 'zm', 'zt',  'zt',  'zm',   'zm',     'zm', 'zt', 'zt',   'zm',   'zm',    'zm',     'zm',  'zm',  'zm']

for iv in range(len(vnmList)-1):

  varName = vnmList[iv]
  suffix  = sufList[iv]

  if (suffix == 'zm'):
    def startIndex(n):
      return 0
  else:
    def startIndex(n):
      return 2**n

  dataList = []
  for i,refineLevel in enumerate(refineLevels):
    dt = dtValues[i]
    dataList.append(Dataset(Filename.format(dt, refineLevel, suffix)))

  n = refineLevels[-1] - refineLevels[0]
  varRef = dataList[-1][varName][:,range(startIndex(n),len(dataList[-1]['altitude'][:]),2**n),0,0]

  varList = []
  for n in range(len(dataList)-1):
    data = dataList[n]
    varList.append(data[varName][:,range(startIndex(n),len(data['altitude'][:]),2**n),0,0])

  z = dataList[0]['altitude'][1:]
  f, ax = plt.subplots(figsize=(8.0,6.6),nrows=1,ncols=1)
  ax.set_prop_cycle('color', color_cycle)

  for tIndex in tIndices:
    rmsList = []
    dzList = []
    for n,var in enumerate(varList):
      dzList.append(1.0/2**n)
      rms = np.sqrt(np.mean((var[tIndex,:]-varRef[tIndex,:])**2))
      rmsList.append(rms)

    # These are the colors that will be used in the plot
    order = np.log(rmsList[-2]/rmsList[-1])/np.log(2)
    ylog = np.log10(rmsList)
    xlog = np.log10(dzList)
    ax.plot(xlog,ylog, '-o', linewidth=2.0, markersize=8.0,
               markeredgecolor='black', markeredgewidth=0.5,
               label='{:.0f}{} ({:.2f})'.format(
                   (dataList[0]['time'][tIndex])/3600.0,'h',order))

  ax.set_title('Convergence of {} ( {} )'.format(varName,case.upper()),fontsize=font_size*1.1,loc='center',x=0.5,y=1.01)
  ax.set_ylabel('{}'.format('log10 (RMSE)'), fontsize=font_size)
  #normalized grid spacing(z^*)
  ax.set_xlabel('{}'.format('log10 (Normalized grid spacing)'), fontsize=font_size)
  #set legend 
  leg = ax.legend(loc='best',fontsize = font_size*0.94,labelspacing=0.24,markerscale=0.8,
            handlelength=1.0,handletextpad=0.5,handleheight=0.8) 
  # set the linewidth of each legend object
  for legobj in leg.legendHandles:
      legobj.set_linewidth(2.0)

  for axis in ['top','bottom','left','right']:
      ax.spines[axis].set_linewidth(0.6)
 
  ax.xaxis.set_minor_locator(MultipleLocator(0.25))
  ax.xaxis.set_major_locator(MultipleLocator(0.5))
  ax.yaxis.set_minor_locator(MultipleLocator(0.25))
  ax.yaxis.set_major_locator(MultipleLocator(0.5))

  ax.tick_params(which='both',  direction='in', pad=6, 
                 width=0.5, labelsize=font_size*0.95, right=True, top=True)
  ax.tick_params(which='major', length=4.5)
  ax.tick_params(which='minor', length=3.0) #, color='b')
  ax.set_xlim(-2.0, 1.0)
  ax.set_aspect('equal')

  #add reference line with slope 1 
  xref = ax.get_xlim()
  yrex = ax.get_ylim()
  yref = xref + np.mean(yrex) - np.mean(xref)
  ax.plot(xref,yref,dashes=[6, 2], linewidth=1.2, color='black')

  #plt.tight_layout(pad=1.5)
  plt.savefig('figure/convergence_{}_{}.png'.format(Figstr,varName))

