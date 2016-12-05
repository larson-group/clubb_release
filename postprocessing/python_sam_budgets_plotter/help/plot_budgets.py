"""
 plot_budgets

 Description:
   help module to grap the data of an .nc file and to plot budgets
"""

import numpy as np
import matplotlib.pyplot as plt
import logging

#-------------------------------------------------------------------------------
#   L O G G E R
#-------------------------------------------------------------------------------
logger = logging.getLogger('plot_budgets')
#logger.setLevel(logging.DEBUG)

#-------------------------------------------------------------------------------
#   F U N C T I O N S
#-------------------------------------------------------------------------------
def plot_budgets(budgets_data, level, xLabel, yLabel, title, name, linewidth = 2, color = 'nipy_spectral'):
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
    """
    # clear the plot
    plt.clf()
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    # set axis labels and title
    ax.set_xlabel(xLabel)
    ax.set_ylabel(yLabel)
    ax.set_title(title, fontsize=18)
    
    # show grid
    ax.grid(True, which='both')
    ax.axhline(y=0, color='k')
    ax.axvline(x=0, color='k')
    
    # color of lines
    cmap = plt.get_cmap(color)
    colors = cmap(np.linspace(0, 1.0, len(budgets_data)))
    
    # plot each variable
    for i in range(len(budgets_data)):
        logger.debug('dimension of %s: %s', budgets_data[i][0], len(budgets_data[i]))
        if budgets_data[i][1]:
        # if it is a help variable, like BUOY e.g., the variable should not be plotted. It is included in B+P variables
            ax.plot(budgets_data[i][3], level, label=budgets_data[i][0], color=colors[i], linewidth=linewidth)
    
    # x axis should be symmetric
    xlimits = ax.get_xlim()
    limit = max(abs(xlimits[0]), abs(xlimits[1]))
    ax.set_xlim(-limit, limit)
    
    # plot the graphs
    plt.legend(loc=1,prop={'size':8})
    plt.savefig(name)

def get_budgets_from_nc(nc, varname, conversion, n, t):
    logger.info('get_budgets_from_nc')
    """
    Input:
      nc         --  Netcdf file object
      varname    --  Variable name string
      conversion --  Conversion factor
      n          --  amount of level
      t          --  amount of timesteps
      n and t are used, if the variable cannot be found

    Output:
      time x height array of the specified variable
    """
    
    keys = nc.variables.keys()
    if varname in keys:
        logger.debug(varname)
        var = nc.variables[varname]        
        var = np.squeeze(var)
        var = var*conversion
    else:
        var = np.zeros(shape=(n,t)) - 10000000.
        
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

    var = np.mean(var[idx_t0:idx_t1,idx_z0:idx_z1],axis=0)
    return var