"""
 plotNC

 Description:
   plots the QTADV and QTDIFF budgets of Sam
"""

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import logging

#-------------------------------------------------------------------------------
#   C O N S T A N T S
#-------------------------------------------------------------------------------
DAY = 24
HOUR = 3600
KG = 1000
FACTOR = 1. / (DAY * HOUR * KG)

#-------------------------------------------------------------------------------
#   L O G G E R
#-------------------------------------------------------------------------------
logger = logging.getLogger('plot_budgets')
#logger.setLevel(logging.DEBUG)

#-------------------------------------------------------------------------------
#   F U N C T I O N S
#-------------------------------------------------------------------------------
def plot_many_budgets(budgets_data, plot_data, level, yLabel, plot_name, linewidth):
    logger.info('plot_many_budgets')
    
    for plot in plot_data:
        data = []
        names = []
        for i in range(len(budgets_data)):
            if plot[0] == budgets_data[i][0] and budgets_data[i][2]:
                data.append(budgets_data[i][4])
                names.append(budgets_data[i][1])
        plot_budgets(data, level, plot[2], yLabel, names, plot[1], plot_name + str(plot[0]) + ".jpg", linewidth)

def plot_budgets(budgets_data, level, xLabel, yLabel, budgets_name, title, name, linewidth):
    logger.info('plot_budgets')
    """
    Plots a plot with budgets
    Input:
      budgets_data   --  list of budgets
      level          --  levels
      xLabel         --  label of the x axis
      yLabel         --  label of the y axis
      budgets_name   --  list of names for the legend
      title          --  title of the plot
      name           --  name of the file
    """
    plt.clf()
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    ax.set_xlabel(xLabel)
    ax.set_ylabel(yLabel)
    ax.set_title(title, fontsize=18)
    
    ax.grid(True, which='both')
    ax.axhline(y=0, color='k')
    ax.axvline(x=0, color='k')
    
    cmap = plt.get_cmap('jet')
    colors = cmap(np.linspace(0, 1.0, len(budgets_name)))
    
    for i in range(len(budgets_data)):
        logger.debug('dimension of %s: %s', budgets_name[i], len(budgets_data[i]))
        ax.plot(budgets_data[i], level, label=budgets_name[i], color=colors[i], linewidth=linewidth)
        xlimits = ax.get_xlim()
        limit = max(abs(xlimits[0]), abs(xlimits[1]))
        ax.set_xlim(-limit, limit)
    
    plt.legend(loc=1,prop={'size':8})
    plt.savefig(name)

def get_budgets_from_nc(nc, varname, conversion, n, t):
    logger.info('get_budgets_from_nc')
    """
    Input:
      nc         --  Netcdf file object
      varname    --  Variable name string
      conversion --  Conversion factor

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
        logger.warning("Could not find the variable %s", varname)
        var = np.zeros(shape=(n,t)) - 10000000.
        
    
    #var = nc.variables[varname]        
    #var = np.squeeze(var)
    #var = var*conversion
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

def get_all_budgets_mean(nc, varnames, conversions, idx_t0, idx_t1, idx_z0, idx_z1):
    logger.info('get_all_budgets_mean')
    n = len(varnames)
    budgets_data = []
    for i in range(n):
        budgets_data.append(mean_profiles(get_budgets_from_nc(nc, varnames[i], conversions[i]), idx_t0, idx_t1, idx_z0, idx_z1))
    return budgets_data