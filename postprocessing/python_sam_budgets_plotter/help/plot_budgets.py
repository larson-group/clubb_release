"""
 plot_budgets

 Description:
   help module to grap the data of an .nc file and to plot budgets
"""

import numpy as np
import matplotlib.pyplot as plt
import logging
from matplotlib.ticker import ScalarFormatter as stick

#-------------------------------------------------------------------------------
#   L O G G E R
#-------------------------------------------------------------------------------
logger = logging.getLogger('plot_budgets')
#logger.setLevel(logging.DEBUG)


#-------------------------------------------------------------------------------
#   D E F I N I T I O N S
#-------------------------------------------------------------------------------
#styles = ['-']
styles = ['--','-.',':']
ticks = stick()
ticks.set_powerlimits((-3,3))
legend_pos = [2,1]
#-------------------------------------------------------------------------------
#   F U N C T I O N S
#-------------------------------------------------------------------------------
def plot_budgets(budgets_data, level, xLabel, yLabel, title, name, linewidth = 1, color = 'nipy_spectral',pdf=None):
    logger.info('plot_budgets')
    """
    Plots a plot with budgets
    Input:
      budgets_data   --  list of budgets (names (index 0) and values (index 3))
      level          --  levels
      xLabel         --  label of the x axis
      yLabel         --  label of the y axis
      title          --  title of the plot
      name           --  name of the file
      linewidth      --  linewidth of the budgets
      color          --  name of colormap
      pdf            --  axes object for output in pdf
      xax            --  switch to plot on different x-scale (0=bottom, 1=top)
    """
    # clear the plot
    #plt.clf()
    fig = plt.figure(figsize=(6,10))
    ax = fig.add_subplot(111)
    if any([b[3]==1 for b in budgets_data]):
        axes = [ax, ax.twiny()]
        titlepos = 1.04
    else:
        axes = [ax]
        titlepos = 1.01
    
    # set axis labels and title
    ax.set_xlabel(xLabel)
    ax.set_ylabel(yLabel)
    ax.set_title(title, fontsize=18,y=titlepos)
    
    # show grid
    #ax.grid(True, which='both')
    ax.axvline(x=0, color='k',ls="--",alpha=.5)
    ax.axhline(y=0, color='k',ls="--",alpha=.5)
    if pdf is not None:
        pdf[0].axhline(y=0, color='k',ls="--",alpha=.5)
        pdf[0].axvline(x=0, color='k',ls="--",alpha=.5)
    
    # color of lines
    cmap = plt.get_cmap(color)
    # getting as many colors as lines to plot
    counter_colors = 0
    for budget in budgets_data:
        if budget[1]:
            counter_colors += 1
    colors = cmap(np.linspace(0, 1.0, counter_colors))
    
    # plot each variable
    j = 0
    lines = []
    pdf_lines = []
    for i in range(len(budgets_data)):
        logger.debug('dimension of %s: %s', budgets_data[i][0], len(budgets_data[i]))
        if budgets_data[i][1]:
            # if it is a help variable, like BUOY e.g., the variable should not be plotted. It is included in B+P variables
            axes[budgets_data[i][3]].plot(budgets_data[i][4][2:], level[2:], label=budgets_data[i][0], color=colors[j], linewidth=linewidth, ls=styles[j%len(styles)])
            if (pdf is not None):
                #pdf.plot(budgets_data[i][3][2:], level[2:], label=budgets_data[i][0], color=colors[j], linewidth=linewidth, ls=styles[j%len(styles)])
                pdf[budgets_data[i][3]].plot(budgets_data[i][4][2:], level[2:], label=budgets_data[i][0], color=colors[j], linewidth=linewidth, ls=styles[j%len(styles)])
            j += 1
            #limit = max(absbudgets_data[i][3][1:])*5
            #xlimits = ax.get_xlim()
            #limit = max(abs(xlimits[0]), abs(xlimits[1]))
            #ax.set_xlim(-limit, limit)
    
    # x axis should be symmetric
    ticklist = []
    for i in range(len(axes)):
        xlimits = axes[i].get_xlim()
        limit = max(abs(xlimits[0]), abs(xlimits[1]))
        #e = math.floor(math.log10(limit))
        #limit = round(limit/10**e)
        axes[i].set_xlim(-limit, limit)
        axes[i].xaxis.set_major_formatter(ticks)
        axes[i].grid(True)
        axes[i].legend(loc=legend_pos[i], prop={'size':8})
        ticklist.append(axes[i].get_xticks())
        if pdf is not None:
            pdf[i].set_xlim(-limit, limit)
            pdf[i].xaxis.set_major_formatter(ticks)
            pdf[i].grid(True)
            pdf[i].legend(loc=legend_pos[i], prop={'size':8})
    
    if len(axes)>1:
        ls = map(len,ticklist)
        l = max(ls)
        for i in range(len(axes)):
            axes[i].set_xticks(np.linspace(ticklist[i][0],ticklist[i][-1],l))
            if pdf is not None:
                pdf[i].set_xticks(np.linspace(ticklist[i][0],ticklist[i][-1],l))
    
    # plot the graphs
    #labels = [l.get_label() for l in lines]
    fig.savefig(name)
    plt.close(fig)
    if pdf is not None:
        #labels = [l.get_label() for l in pdf_lines]
        pdf[0].set_xlabel(xLabel)
        pdf[0].set_ylabel(yLabel)
        pdf[0].set_title(title, fontsize=18,y=titlepos)
        #pdf[0].grid(True, which='both')

def get_budgets_from_nc(nc, varname, conversion, n, t):
    logger.info('get_budgets_from_nc:%s', varname)
    """
    Input:
      nc         --  Netcdf file object
      varname    --  Variable name string
      conversion --  Conversion factor
      n          --  amount of level
      t          --  amount of timesteps
      n and t are used, if the variable cannot be found

    Output:
      time x height array of the specified variable, scaled by conversion factor
    """
    
    keys = nc.variables.keys()
    if varname in keys:
        logger.debug('%s is in keys', varname)
        var = nc.variables[varname]
        var = np.squeeze(var)
        var = var*conversion
    else:
        logger.debug('%s is not in keys', varname)
        var = np.zeros(shape=(n,t)) - 1000.
        
    return var
                
def mean_profiles(var, idx_t0, idx_t1, idx_z0, idx_z1):
    logger.info('mean_profiles')
    """
    Input:
      var    -- time x height array of some property
      idx_t0 -- Index corrosponding to the beginning of the averaging interval
      idx_t1 -- Index corrosponding to the end of the averaging interval
      idx_z0 -- Index corrosponding to the lowest model level of the averaging interval
      idx_z1 -- Index corrosponding to the highest model level of the averaging interval

    Output:
      var    -- time averaged vertical profile of the specified variable
    """
    
    var2 = np.nanmean(var[idx_t0:idx_t1,idx_z0:idx_z1],axis=0)
    return var2
    
def get_units(nc, varname):
    logger.info('get_units:%s', varname)
    """
    Input:
      nc         --  Netcdf file object
      varname    --  Variable name string
    Output:
      unit as string
    """
    
    keys = nc.variables.keys()
    if varname in keys:
        logger.debug('%s is in keys', varname)
        unit = nc.variables[varname].units
    else:
        logger.debug('%s is not in keys', varname)
        unit = "nm"
    
    return unit
    
def get_long_name(nc, varname):
    logger.info('get_long_name:%s', varname)
    """
    Input:
      nc         --  Netcdf file object
      varname    --  Variable name string
    Output:
      long_name as string
    """
    
    keys = nc.variables.keys()
    if varname in keys:
        logger.debug('%s is in keys', varname)
        long_name = nc.variables[varname].long_name
    else:
        logger.debug('%s is not in keys', varname)
        long_name = "nm"
    
    return long_name
    
    