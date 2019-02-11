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
logger = logging.getLogger('plotgen.help.pb')
#logger.setLevel(logging.INFO)
logger.setLevel(logging.DEBUG)
#logger.setLevel(logging.CRITICAL)


#-------------------------------------------------------------------------------
#   D E F I N I T I O N S
#-------------------------------------------------------------------------------
#styles = ['-']
styles = ['-','--','-.',':']
#ticks = stick(useMathText=True)
#ticks.set_powerlimits((-3,3))
pow_lim = 1
legend_pos = [2,1]
legend_title = ['bottom', 'top']
# Create 9 maximally distinguishable colors TODO!!
#             red       blue      green     purple    brown     black     grey      orange    magenta
color_arr = ['#e41a1c','#377eb8','#4daf4a','#984ea3','#a65628','#000000','#999999','#ffa600','#f750bf']

## Comparison styles
comp_style = {
    'clubb' : {
        'color' : 'red',
        'lw'    : 2,
        'ls'    : '--',
        'label' : 'new CLUBB'
        },
    'sam'   : {
        'color' : 'black',
        'lw'    : 5,
        'ls'    : '-',
        'label' : 'SAM-LES'
        },
    'old'   : {
        'color' : 'green',
        'lw'    : 2,
        'ls'    : ':',
        'label' : 'old CLUBB'
        },
    }

fontsizes = {
    'labels' : 25,
    'ticks' : 20,
    'title' : 30,
    'legend' : 18,
    }

#-------------------------------------------------------------------------------
#   F U N C T I O N S
#-------------------------------------------------------------------------------
def plot_budgets(budgets_data, level, xLabel, yLabel, title, name, lw = 5, grid = True,  color = 'nipy_spectral', pdf=None):
    logger.info('plot_budgets')
    """
    This function is deprecated and has been replaced by plot_profiles below.
    Plots a plot with budgets
    Input:
      budgets_data   --  list of budgets (names (index 0) and values (index 3))
      level          --  levels
      xLabel         --  label of the x axis
      yLabel         --  label of the y axis
      title          --  title of the plot
      name           --  name of the file
      lw             --  linewidth of the budgets
      color          --  name of colormap
      pdf            --  axes object for output in pdf
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
    ax.set_title(title, fontsize=18, y=titlepos)
    
    # show grid
    #ax.grid(True, which='both')
    ax.axvline(x=0, color='k',ls="--", alpha=.5)
    ax.axhline(y=0, color='k',ls="--", alpha=.5)
    if pdf is not None:
        pdf[0].axhline(y=0, color='k',ls="--", alpha=.5)
        pdf[0].axvline(x=0, color='k',ls="--", alpha=.5)
    
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
        logger.debug('dimension of %s: %s', budgets_data[i][0], len(budgets_data[i][3]))
        if budgets_data[i][1]:
            # if it is a help variable, like BUOY e.g., the variable should not be plotted. It is included in B+P variables
            axes[budgets_data[i][3]].plot(budgets_data[i][4][2:], level[2:], label=budgets_data[i][0], color=colors[j], lw=lw, ls=styles[j%len(styles)])
            if (pdf is not None):
                #pdf.plot(budgets_data[i][3][2:], level[2:], label=budgets_data[i][0], color=colors[j], lw=lw, ls=styles[j%len(styles)])
                pdf[budgets_data[i][3]].plot(budgets_data[i][4][2:], level[2:], label=budgets_data[i][0], color=colors[j], lw=lw, ls=styles[j%len(styles)])
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
        ticks = stick(useMathText=True)
        ticks.set_powerlimits((-pow_lim,pow_lim))
        axes[i].xaxis.set_major_formatter(ticks)
        axes[i].grid(grid)
        axes[i].legend(loc=legend_pos[i], prop={'size':8})
        ticklist.append(axes[i].get_xticks())
        if pdf is not None:
            pdf[i].set_xlim(-limit, limit)
            pdf[i].xaxis.set_major_formatter(ticks)
            pdf[i].grid(grid)
            pdf[i].legend(loc=legend_pos[i], prop={'size':8})
    
    if len(axes)>1:
        ls = map(len,ticklist)
        l = max(ls)
        for i in range(len(axes)):
            axes[i].set_xticks(np.linspace(ticklist[i][0],ticklist[i][-1],l))
            if pdf is not None:
                pdf[i].set_xticks(np.linspace(ticklist[i][0],ticklist[i][-1],l))
    # plot the graphs
    fig.savefig(name)
    plt.close(fig)
    if pdf is not None:
        pdf[0].set_xlabel(xLabel)
        pdf[0].set_ylabel(yLabel)
        pdf[0].set_title(title, fontsize=18, y=titlepos)
    


def plot_profiles(data, level, xLabel, yLabel, title, name, startLevel = 0, lw = 5, grid = True, color = 'nipy_spectral', pdf=None):
    """
    This function replaced plot_budgets, so it should be used from now on.
    Plots a plot with budgets
    TODO: Change line label
    Input:
    data           --  list of data per plot line
                       Structure of data:
                       Index | Description
                       ______________________________________________
                       0   | plot label
                       1   | switch: - True  = show line in plot
                           |         - False = do not show line
                       2   | plot axis (0: left, 1: right)
                       3   | data of variable (numpy array)
    level          --  height levels
    xLabel         --  label of the x axis
    yLabel         --  label of the y axis
    title          --  title of the plot
    name           --  name of the file
    startLevel     --  Height level at which to plot lines, values below this level will be ignored
    lw             --  linewidth of the plot lines
    color          --  name of colormap (deprecated, as color sequence is hard coded at the momemt)
    pdf            --  PdfPages object
    """
    logger.info('plot_profiles')
    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111)
    #logger.debug([b[2] for b in data])
    if any([b[2]==1 for b in data]):
        axes = [ax, ax.twiny()]
        titlepos = 1.04
    else:
        axes = [ax]
        titlepos = 1.01
            
    # set axis labels and title
    ax.set_xlabel(xLabel, fontsize=fontsizes['labels'])
    ax.set_ylabel(yLabel, fontsize=fontsizes['labels'])
    ax.set_title(title, fontsize=fontsizes['title'], y=titlepos)
    
    # show grid
    #ax.grid(grid, which='both')
    ax.axvline(x=0, color='k',ls="--", alpha=.5)
    ax.axhline(y=0, color='k',ls="--", alpha=.5)
        
    # color of lines
    cmap = plt.get_cmap(color)
    # getting as many colors as lines to plot (deprecated)
    counter_colors = 0
    for budget in data:
        if budget[1]:
            counter_colors += 1
    colors = cmap(np.linspace(0, 1.0, counter_colors))
    
    # loop over lines in plot
    j = 0
    for i in range(len(data)):
        if data[i][3] is None and data[i][1]:
            logger.debug('Add dummy line for %s', data[i][0])
            # plot dummy line, even necessary?
            #axes[data[i][2]].plot([])
            # increment line counter j
            j+=1
        else:
            logger.debug('dimension of %s: %s', data[i][0], len(data[i][3]))
            # if it is a help variable, like BUOY e.g., the variable should not be plotted. It is included in B+P variables
            if data[i][1]:
                # plot line
                axes[data[i][2]].plot(data[i][3][startLevel:], level[startLevel:], label=data[i][0], color=color_arr[j%len(color_arr)], lw=lw, ls=styles[j%len(styles)])
                # Change color for next line
                j += 1
            
    # x axis should be symmetric
    # List to save ticks of each axis
    ticklist = []
    for i in range(len(axes)):
        xlimits = axes[i].get_xlim()
        # Calculate absolute maximum of x values
        limit = max(abs(xlimits[0]), abs(xlimits[1]))
        # Round to next bigger reasonable number, not needed?
        #e = math.floor(math.log10(limit))
        #limit = round(limit/10**e)
        # Set x limits symmetrically
        axes[i].set_xlim(-limit, limit)
        # Set tick label format to scalar formatter with length sensitive format (switch to scientific format when a set order of magnitude is reached)
        ticks = stick(useMathText=True)
        ticks.set_powerlimits((-pow_lim,pow_lim))
        axes[i].xaxis.set_major_formatter(ticks)
        axes[i].grid(grid, which='both')
        # Generate legend
        axes[i].legend(loc=legend_pos[i], prop={'size': fontsizes['legend']}, title=legend_title[i] if len(axes)>1 else None)
        # Add ticks to list
        ticklist.append(axes[i].get_xticks())

    # For multiple axes adjust ticks to coincide
    if len(axes)>1:
        ls = map(len,ticklist)
        l = max(ls)
        for i in range(len(axes)):
            axes[i].set_xticks(np.linspace(ticklist[i][0],ticklist[i][-1],l))
    
    # Increase fontsize for ticklabels
    for i in range(len(axes)):
        axes[i].tick_params(axis='both', labelsize=fontsizes['ticks'])
        axes[i].xaxis.get_offset_text().set_size(fontsizes['ticks'])
    fig.canvas.draw_idle()
                
    # Save plot
    fig.savefig(name)
    if pdf:
        pdf.savefig(fig)
    # Close figure
    plt.close(fig)


def plot_comparison(data_clubb, data_sam, level_clubb, level_sam, xLabel, yLabel, title, name, startLevel = 0, grid = False, pdf=None, plot_old_clubb=False, data_old=None, level_old=None):
    """
    Plots a plot with budgets
    Input:
    data_clubb      --  list of CLUBB data per plot line
    data_sam        --  list of SAM data per plot line
    Structure of data:
    Index | Description
    ______________________________________________
    0   | plot label
    1   | switch: - True  = show line in plot
        |         - False = do not show line
    2   | plot axis (0: left, 1: right)
    3   | data of variable (numpy array)
    level_clubb    --  height levels for CLUBB model
    level_sam      --  height levels for SAM model
    xLabel         --  label of the x axis
    yLabel         --  label of the y axis
    title          --  title of the plot
    name           --  name of the file
    startLevel     --  
    grid           --  switch: True = shpw grid for both axes, False=do not show grid (default)
    pdf            --  PdfPages object
    """
    logger.info('plot_profiles')
    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111)
    # A 2nd axis should not be necessary here, as there will only be one line for SAM and CLUBB each per plot
    #if any([b[2]==1 for b in data]):
        #axes = [ax, ax.twiny()]
        #titlepos = 1.04
    #else:
        #axes = [ax]
    titlepos = 1.01
    
    # set axis labels and title
    ax.set_xlabel(xLabel, fontsize=fontsizes['labels'])
    ax.set_ylabel(yLabel, fontsize=fontsizes['labels'])
    ax.set_title(title, fontsize=fontsizes['title'], y=titlepos)

    # show grid
    ax.grid(grid)
    ax.axvline(x=0, color='k',ls="--",alpha=.5)
    ax.axhline(y=0, color='k',ls="--",alpha=.5)
    
    # loop over lines in plot
    for i in range(len(data_clubb)):
        logger.debug('dimension of CLUBB data %s: %s', data_clubb[i][0], len(data_clubb[i][3]))
        logger.debug('dimension of SAM data %s: %s', data_sam[i][0], len(data_sam[i][3]))
        if plot_old_clubb:
            logger.debug('dimension of old CLUBB data %s: %s', data_old[i][0], len(data_old[i][3]))
        # if it is a help variable, like BUOY e.g., the variable should not be plotted. It is included in B+P variables
        if data_clubb[i][1]:
            ## plot lines
            # SAM
            d = comp_style['sam']
            ax.plot(data_sam[i][3][startLevel:], level_sam[startLevel:], label=d['label'], color=d['color'], lw=d['lw'], ls=d['ls'])
            # CLUBB
            d = comp_style['clubb']
            ax.plot(data_clubb[i][3][startLevel:], level_clubb[startLevel:], label=d['label'], color=d['color'], lw=d['lw'], ls=d['ls'])
            logger.debug("Data limits: CLUBB: (%f, %f); SAM: (%f, %f)", data_clubb[i][3].min(), data_clubb[i][3].max(), data_sam[i][3].min(), data_sam[i][3].max())
            xmin = min(data_sam[i][3].min(), data_clubb[i][3].min())
            xmax = max(data_sam[i][3].max(), data_clubb[i][3].max())
            # old CLUBB
            if plot_old_clubb:
                d = comp_style['old']
                ax.plot(data_old[i][3][startLevel:], level_old[startLevel:], label=d['label'], color=d['color'], lw=d['lw'], ls=d['ls'])
                xmin = min(xmin, data_old[i][3].min())
                xmax = max(xmax, data_old[i][3].max())
            xmin -= .1*(xmax-xmin)
            xmax += .1*(xmax-xmin)
            logger.debug("Set xlim to (%f, %f)", xmin, xmax)
            if xmax-xmin >0:
                ax.set_xlim(xmin, xmax)
    # Set tick label format to scalar formatter with length sensitive format (switch to scientific format when a set order of magnitude is reached)
    ticks = stick(useMathText=True)
    ticks.set_powerlimits((-pow_lim,pow_lim))
    ax.xaxis.set_major_formatter(ticks)
    # Increase fontsize for ticklabels
    ax.tick_params(axis='both', labelsize=fontsizes['ticks'])
    ax.xaxis.get_offset_text().set_size(fontsizes['ticks'])
    fig.canvas.draw_idle()
    # Add legend
    ax.legend(loc=0, prop={'size': fontsizes['legend']})

    # Save plot
    fig.savefig(name)
    if pdf:
        pdf.savefig(fig)
    # Close figure
    plt.close(fig)


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

def get_var_from_nc(nc, varname, conversion, n, t):
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