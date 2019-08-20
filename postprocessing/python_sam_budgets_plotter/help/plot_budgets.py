"""
 plot_budgets

 Description:
   help module to grap the data of an .nc file and to plot budgets
"""

import numpy as np
import matplotlib.pyplot as plt
import logging
from matplotlib.ticker import ScalarFormatter as stick
# Imports used for moving power offset
import types
import matplotlib.transforms as mpt
from plot_defs import *

#-------------------------------------------------------------------------------
#   L O G G E R
#-------------------------------------------------------------------------------
logger = logging.getLogger('plotgen.help.pb')
#logger.setLevel(logging.INFO)
logger.setLevel(logging.DEBUG)
#logger.setLevel(logging.CRITICAL)

#-------------------------------------------------------------------------------
#   F U N C T I O N S
#-------------------------------------------------------------------------------
# DEPRECATED!!!
def plot_budgets(budgets_data, level, xLabel, yLabel, title, name, lw=5, grid=True,  color='nipy_spectral', pdf=None):
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
        ticks.set_powerlimits(pow_lims)
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
    
# END plot_budgets (DEPRECATED!)

def plot_profiles(data, level, xLabel, yLabel, title, name, textEntry="", textPos=(0,0), startLevel=0, lw=5, grid=True, color='nipy_spectral', centering=False, pdf=None):
    """
    NOTE: This function replaced plot_budgets, so it should be used from now on.
    Creates a height profile plot with one line per entry in data
    Input:
    data           --  list of data per plot line
                       Structure of data:
                       Index | Description
                       ______________________________________________
                       0   | plot label
                       1   | switch: - True  = show line in plot
                           |         - False = do not show line
                       2   | netcdf variable name
                       3   | plot axis (0: left, 1: right)
                       4   | data of variable (numpy array)
    level          --  height levels
    xLabel         --  label of the x axis
    yLabel         --  label of the y axis
    title          --  title of the plot
    name           --  name of the file
    textEntry      --  Additional text to be put into plot
    textPos        --  Position of additional text in data coordinates
    startLevel     --  Height level at which to plot lines, values below this level will be ignored (Redundant with startHeight entry in case file)
    lw             --  linewidth of the plot lines
    color          --  name of colormap (deprecated, as color sequence is hard coded at the momemt)
    pdf            --  PdfPages object
    """
    logger.info('plot_profiles')
    fig = plt.figure(figsize=profile_figsize)
    ax = fig.add_subplot(111)
    #logger.debug([b[2] for b in data])
    if any([b[3]==1 for b in data]):
        axes = [ax, ax.twiny()]
        for i, el in enumerate(axes):
            el.xaxis._update_offset_text_position = types.MethodType(x_update_offset_text_position, el.xaxis)
        axes[1].xaxis.tick_top()
        axes[1].xaxis.offset_text_position = 'top'
        titlepos = 1.06
        centering=True
    else:
        axes = [ax]
        titlepos = 1.01
            
    # set axis labels, title, and additional text
    ax.set_xlabel(xLabel, fontsize=fontsizes['labels'])
    ax.set_ylabel(yLabel, fontsize=fontsizes['labels'])
    ax.set_title(title, fontsize=fontsizes['title'], y=titlepos)
    #ax.text(0.1,0.9,'c)',fontsize=30,transform=ax.transAxes)
    ax.text(textPos[0], textPos[1], textEntry, fontsize=30, transform=ax.transAxes)
    
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
    lims = np.full((len(axes),2), fill_value=np.nan)
    #logger.debug(lims)
    for i in range(len(data)):
        if data[i][4] is None and data[i][1]:
            logger.debug('Add dummy line for %s', data[i][0])
            # plot dummy line, even necessary?
            #axes[data[i][2]].plot([])
            # increment line counter j
            j+=1
        else:
            logger.debug('dimension of %s: %s', data[i][0], len(data[i][4]))
            # if it is a help variable, like BUOY e.g., the variable should not be plotted. It is included in B+P variables
            if data[i][1]:
                # plot line
                axes[data[i][3]].plot(data[i][4][startLevel:], level[startLevel:], label=data[i][0], color=color_arr[j%len(color_arr)], lw=lw, ls=styles[j%len(styles)])
                if not np.all(np.isnan(data[i][4])):
                    lims[data[i][3],0] = np.nanmin( (np.nanmin(data[i][4]), lims[data[i][3],0]) )
                    lims[data[i][3],1] = np.nanmax( (np.nanmax(data[i][4]), lims[data[i][3],1]) )
                # Change color for next line
                j += 1

    #logger.debug(lims)
    # x axis should be symmetric
    # List to save ticks of each axis
    logger.info('Plot formatting')
    ticklist = []
    for i in range(len(axes)):
        #logger.debug('Adjusting plot limits')
        if centering:
            xlimits = axes[i].get_xlim()
            # Calculate absolute maximum of x values
            limit = max(abs(xlimits[0]), abs(xlimits[1]))
            # Round to next bigger reasonable number, not needed?
            #e = math.floor(math.log10(limit))
            #limit = round(limit/10**e)
            # Set x limits symmetrically
            axes[i].set_xlim(-limit, limit)
        else:
            margin = np.abs(lims[i,1]-lims[i,0])*.1
            if not np.any(np.isnan(lims[i])):
                axes[i].set_xlim(lims[i]+np.array((-margin,margin)))
        # Set tick label format to scalar formatter with length sensitive format (switch to scientific format when a set order of magnitude is reached)
        #logger.debug('Setting tick labels to scientific notation')
        ticks = stick(useMathText=True)
        ticks.set_powerlimits(pow_lims)
        axes[i].xaxis.set_major_formatter(ticks)
        axes[i].grid(grid, which='both')
        # Generate legend
        axes[i].legend(loc=legend_pos[i], prop={'size': fontsizes['legend']}, title=legend_title[i] if len(axes)>1 else None)
        # Add ticks to list
        ticklist.append(axes[i].get_xticks())

    # For multiple axes adjust ticks to coincide
    if len(axes)>1:
        #logger.debug('Matching ticks for multiple x axes')
        ls = map(len,ticklist)
        l = max(ls)
        for i in range(len(axes)):
            axes[i].set_xticks(np.linspace(ticklist[i][0],ticklist[i][-1],l))
    
    # Increase fontsize for ticklabels
    for i in range(len(axes)):
        axes[i].tick_params(axis='both', labelsize=fontsizes['ticks'])
        axes[i].xaxis.get_offset_text().set_size(fontsizes['ticks'])
    #logger.debug('Redrawing figure')
    fig.canvas.draw_idle()
                
    # Save plot
    #logger.debug('Saving to jpg')
    fig.savefig(name)
    if pdf:
        #logger.debug('Saving to pdf')
        pdf.savefig(fig)
    # Close figure
    plt.close(fig)


def plot_comparison(data_clubb, data_sam, level_clubb, level_sam, xLabel, yLabel, title, name, textEntry="", textPos=(0,0), startLevel=0, grid=False, pdf=None, plot_old_clubb=False, data_old=None, level_old=None):
    """
    Plots a plot with budgets
    Input:
    data_clubb     --  list of CLUBB data per plot line
    data_sam       --  list of SAM data per plot line
    data           --  list of data per plot line
                        Structure of data:
                        Index | Description
                        ______________________________________________
                        0   | plot label
                        1   | switch: - True  = show line in plot
                            |         - False = do not show line
                        2   | netcdf variable name
                        3   | plot axis (0: left, 1: right)
                        4   | data of variable (numpy array)
    level_clubb    --  height levels for CLUBB model
    level_sam      --  height levels for SAM model
    xLabel         --  label of the x axis
    yLabel         --  label of the y axis
    title          --  title of the plot
    name           --  name of the file
    textEntry      --  Additional text to be put into plot
    textPos        --  Position of additional text in data coordinates
    startLevel     --  Height level at which to plot lines, values below this level will be ignored (Redundant with startHeight entry in case file)
    grid           --  switch: True = show grid for both axes, False=do not show grid (default)
    pdf            --  PdfPages object
    """
    logger.info('plot_profiles')
    fig = plt.figure(figsize=profile_figsize)
    ax = fig.add_subplot(111)
    # A 2nd axis should not be necessary here, as there will only be one line for SAM and CLUBB each per plot
    #if any([b[2]==1 for b in data]):
        #axes = [ax, ax.twiny()]
        #titlepos = 1.04
    #else:
        #axes = [ax]
    titlepos = 1.01
    
    # set axis labels, title, and additional text
    ax.set_xlabel(xLabel, fontsize=fontsizes['labels'])
    ax.set_ylabel(yLabel, fontsize=fontsizes['labels'])
    ax.set_title(title, fontsize=fontsizes['title'], y=titlepos)
    ax.text(textPos[0], textPos[1], textEntry, fontsize=30, transform=ax.transAxes, zorder=100)

    # show grid
    ax.grid(grid)
    ax.axvline(x=0, color='k',ls="--",alpha=.5)
    ax.axhline(y=0, color='k',ls="--",alpha=.5)
    
    # loop over lines in plot
    for i in range(len(data_clubb)):
        logger.debug('dimension of CLUBB data %s: %s', data_clubb[i][0], len(data_clubb[i][4]))
        logger.debug('dimension of SAM data %s: %s', data_sam[i][0], len(data_sam[i][4]))
        if plot_old_clubb:
            logger.debug('dimension of old CLUBB data %s: %s', data_old[i][0], len(data_old[i][4]))
        # if it is a help variable, like BUOY e.g., the variable should not be plotted. It is included in B+P variables
        if data_clubb[i][1]:
            ## plot lines
            # SAM
            d = comp_style['sam']
            ax.plot(data_sam[i][4][startLevel:], level_sam[startLevel:], label=d['label'], color=d['color'], lw=d['lw'], ls=d['ls'])
            # CLUBB
            d = comp_style['clubb']
            ax.plot(data_clubb[i][4][startLevel:], level_clubb[startLevel:], label=d['label'], color=d['color'], lw=d['lw'], ls=d['ls'])
            logger.debug("Data limits: CLUBB: (%f, %f); SAM: (%f, %f)", data_clubb[i][4].min(), data_clubb[i][4].max(), data_sam[i][4].min(), data_sam[i][4].max())
            xmin = min(data_sam[i][4].min(), data_clubb[i][4].min())
            xmax = max(data_sam[i][4].max(), data_clubb[i][4].max())
            # old CLUBB
            if plot_old_clubb:
                d = comp_style['old']
                ax.plot(data_old[i][4][startLevel:], level_old[startLevel:], label=d['label'], color=d['color'], lw=d['lw'], ls=d['ls'])
                xmin = min(xmin, data_old[i][4].min())
                xmax = max(xmax, data_old[i][4].max())
            xmin -= .1*(xmax-xmin)
            xmax += .1*(xmax-xmin)
            logger.debug("Set xlim to (%f, %f)", xmin, xmax)
            if xmax-xmin >0:
                ax.set_xlim(xmin, xmax)
    # Set tick label format to scalar formatter with length sensitive format (switch to scientific format when a set order of magnitude is reached)
    ticks = stick(useMathText=True)
    ticks.set_powerlimits(pow_lims)
    ax.xaxis.set_major_formatter(ticks)
    # Increase fontsize for ticklabels
    ax.tick_params(axis='both', labelsize=fontsizes['ticks'])
    ax.xaxis.get_offset_text().set_size(fontsizes['ticks'])
    fig.canvas.draw_idle()
    # Add legend
    ax.legend(loc=legend_pos[0], prop={'size': fontsizes['legend']})

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
        var = np.zeros(shape=(t,n)) - 1000.
        
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
        var = np.zeros(shape=(t,n)) - 1000.
        
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
    

# As the tick_top function of Axis currently does not move the power offset, this is a fix until it is officially implemented in pyplot
def x_update_offset_text_position(self, bboxes, bboxes2):
    logger.info('x_update_offset_text_position')
    x, y = self.offsetText.get_position()
    logger.debug((x,y))
    
    if self.offset_text_position == 'bottom':
        logger.debug('bottom')
        if bboxes:
            bbox = mpt.Bbox.union(bboxes)
        else:
            bbox = self.axes.bbox
        y = bbox.ymin - self.OFFSETTEXTPAD * self.figure.dpi / 72.0
        self.offsetText.set(va='top', ha='right')

    else:
        logger.debug('top')
        if bboxes2:
            bbox = mpt.Bbox.union(bboxes2)
        else:
            bbox = self.axes.bbox
        y = bbox.ymax + self.OFFSETTEXTPAD * self.figure.dpi / 72.0 - .75 * self.figure.dpi
        x = x * (1.007)
        self.offsetText.set(va='bottom', ha='left')
        
    logger.debug((x,y))
    self.offsetText.set_position((x, y))