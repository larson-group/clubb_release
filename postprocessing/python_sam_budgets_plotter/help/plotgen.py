
# General imports
import os
import sys
import re
import logging
from datetime import datetime as dt
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as mpp
import matplotlib.animation as manim
from matplotlib.backends.backend_pdf import PdfPages
from netCDF4 import Dataset
from PyPDF2 import PdfFileWriter as pdfw, PdfFileReader as pdfr

# plotgen imports
import plot_budgets as pb
import OutputWriter as ow
from plot_defs import *

#-------------------------------------------------------------------------------
#   L O G G E R
#-------------------------------------------------------------------------------
#formatter = logging.Formatter('%(asctime)s:%(levelname)s:%(message)s')
#fileHandler = logging.FileHandler('plotgen.log')
#streamHandler = logging.StreamHandler()
#fileHandler.setFormatter(formatter)
#streamHandler.setFormatter(formatter)
logger = logging.getLogger('plotgen.help')
#logger.setLevel(logging.INFO)
logger.setLevel(logging.DEBUG)
#logger.setLevel(logging.CRITICAL)


#-------------------------------------------------------------------------------
#   HELPER FUNCTIONS
#-------------------------------------------------------------------------------
def init_plotgen(plots, cf):
    """
    bla
    """
    logger.info("init_plotgen")
    date = dt.now().strftime(date_file_format)
    
    out_dir = cf.out_dir.format(date=date)
    jpg_dir = os.path.join(out_dir, 'jpg')
    logger.info("Creating jpg directories: %s", jpg_dir)
    makeDirectory(jpg_dir)
    plot_case_name = cf.plot_case_name.format(type=plots.name, date=date, plot='{plot}')
    logger.info("Generating jpg template name: %s", plot_case_name)
    
    out_pdf = os.path.join(out_dir, cf.out_pdf.format(type=plots.name, date=date))
    
    return out_dir, jpg_dir, plot_case_name, out_pdf

def load_nc(plots, cf, old_clubb=False):
    """
    bla
    """
    logger.info('load_nc')
    nc_list = []
    if 'clubb_zm' in plots.nc_files:
        logger.info('Loading clubb_zm netcdf')
        try:
            nc_list.append(Dataset(cf.clubb_zm_file,'r'))
        except IOError as e:
            logger.error('The file {}, specified as clubb_zm_file in the {} case file, could not be opened: {}'.format(e.filename, cf.case, e.message))
            sys.exit()
    if 'clubb_zt' in plots.nc_files:
        logger.info('Loading clubb_zt netcdf')
        try:
            nc_list.append(Dataset(cf.clubb_zt_file,'r'))
        except IOError as e:
            logger.error('The file {}, specified as clubb_zt_file in the {} case file, could not be opened: {}'.format(e.filename, cf.case, e.message))
            sys.exit()
    if 'sam' in plots.nc_files:
        logger.info('Loading sam netcdf')
        try:
            nc_list.append(Dataset(cf.sam_file,'r'))
        except IOError as e:
            logger.error('The file {}, specified as sam_file in the {} case file, could not be opened: {}'.format(e.filename, cf.case, e.message))
            sys.exit()
    if 'sam_3d' in plots.nc_files:
        logger.info('Loading sam_3d netcdf')
        try:
            nc_list.append(Dataset(cf.sam_3d_file,'r'))
        except IOError as e:
            logger.error('The file {}, specified as sam_3d_file in the {} case file, could not be opened: {}'.format(e.filename, cf.case, e.message))
            sys.exit()
    if old_clubb:
        logger.info('Loading netcdf for old CLUBB data.')
        try:
            if plots.sortPlots_zm:
                nc_list.append(Dataset(cf.old_clubb_zm_file, 'r'))
            if plots.sortPlots_zt:
                nc_list.append(Dataset(cf.old_clubb_zt_file, 'r'))
        except IOError as e:
            logger.error('The file {}, specified in the {} case file, could not be opened: {}'.format(e.filename, cf.case, e.message))
            sys.exit()
    return nc_list

def get_t_dim(nc, create=False, dt=None):
    """
    TODO:   - implement nsave3Dstart/end
            - match time scaling to consistent units of minutes/seconds
    bla
    """
    logger.info('get_t_dim')
    t = pb.get_var_from_nc(nc, 'time', 1, 0, 0)
    if create:
        if dt is None:
            logger.error('No timestep length specified.')
            sys.exit()
        t = (np.arange(len(t))+1) * dt
    return t

def get_h_dim(nc, model, create=False, cf=None):
    """
    bla
    """
    logger.info('get_h_dim')
    if model=='clubb':
        h = pb.get_var_from_nc(nc, 'altitude', 1, 0, 0)
    else:
        h = pb.get_var_from_nc(nc, 'z', 1, 0, 0)
    #logger.debug(h)
    if create:
        if cf is None:
            logger.error('No case file specified.')
            sys.exit()
        # Get entries from grd file and extrapolate
        try:
            grd = np.genfromtxt(cf.sam_grd)
        except IOError as ioe:
            logger.error('The file {}, specified as sam_grd in the {} case file, could not be opened: {}'.format(ioe.filename, cf.case, ioe.message))
            sys.exit()
        nh = len(h)
        # Generate height levels from grd
        h = np.zeros(nh)
        # Input grd entries
        h[:len(grd)] = grd
        inc = grd[-1]-grd[-2]
        # Fill in rest
        h[len(grd):] = grd[-1] + (np.arange(nh-len(grd))+1) * inc
    return h

def get_values_from_prm(cf):
    """
    bla
    """
    logger.info('get_values_from_prm')
    try:
        with open(cf.sam_prm,'r') as f:
            prm = f.read()
    except IOError as ioe:
        logger.error('The file {}, specified as sam_prm in the {} case file, could not be opened: {}'.format(ioe.filename, cf.case, ioe.message))
        sys.exit()
    prm_vars = {}
    logger.debug(prm)
    for var in prm_patterns:
        logger.debug(var)
        try:
            logger.debug(re.search(prm_patterns[var], prm).group(1))
            prm_vars[var] = float(re.search(prm_patterns[var], prm).group(1))
        except AttributeError as ve:
            logger.error('In file {}, no entry for {} could be found.',cf.sam_prm, var)
            sys.exit()
    return prm_vars


def get_all_variables(nc, lines, plotLabels, nh, nt, t0, t1, h0, h1, filler=0):
    """
    TODO: Include pb.get_units in output structure
    Get variables from netcdf and store in dictionary with the following structure:
    dict entry: 'variable name' : [<plot info>]
    1 entry = 1 plot
    length of list = max number of lines in plot
    plot info structure:
     Index | Description
    ______________________________________________
       0   | plot label
       1   | switch: - True  = show line in plot
           |         - False = do not show line
       2   | plot axis (0: left, 1: right)
       3   | data of variable (numpy array)
    
    List structure in plots.lines:
    lines = [<plot1>,<plot2>,[...]]
    List of all plots, where every plot descriptor has the following structure:
    plot = [<line1>,<line2>,[...]]
    List of lines/variables in plot, where every line has the following structure:
    line = [<label>, <switch>, <expression/var>, <conversion>, <axis>]
    line info structure:
     Index | Description
    ______________________________________________
       0   | line label: Identifier string
       1   | switch: - True  = show line in plot
           |         - False = do not show line
       2   | varname: variable name in nc file
           | expression: expression of other variables in plot
       3   | conversion: applies only to var; conversion factor, applied when reading data from nc file
       4   | axis: 0 = left, 1 = right; axis location parameter, used to plot variables of different scales in one plot
    
    Input:
      nc            --
      lines         --
      plotLabels     --
      nh            --
      nt            --
      t0, t1        --
      h0, h1        --
    Output:
      bla
    """
    logger.info('get_all_variables')
    logger.debug(plotLabels)
    # Initialize list of all plots
    plot_data = []
    # Rename to make some of the expressions in SAM standalone work
    n = nh
    for i,l in enumerate(lines):
        # Initialize list of all lines in plot l
        plot_lines = []
        # Initialize list of function lines in plot l
        functions = []
        logger.debug('Iterate over lines')
        for j,var in enumerate(l):
            if var[2] is None:
                logger.debug("Add dummy line")
                plot_lines.append([var[0], var[1], var[4], None])
            elif not isFunction(var[2]):
                logger.debug(plotLabels[i])
                logger.debug(var[2])
                logger.debug(var[0])
                long_name = pb.get_long_name(nc, var[2])
                logger.debug(long_name)
                logger.info('Plot %s: Variable %s, label %s, long_name %s', plotLabels[i], var[2], var[0], long_name)
                # Fetch data
                data = pb.get_var_from_nc(nc, var[2], var[3], nh, nt)
                logger.debug(data.shape)
                # NAN test: -1000 is the value assigned for invalid data
                if np.any(np.isnan(data)) or np.any(data <= -1000):
                    # if there are invalid data points in the variable replace by filler value
                    logger.warning("Invalid data in variable %s of plot %s.", var[0], plotLabels[i])
                    data = np.where(np.logical_or(np.isnan(data), data<=-9000), filler, data)
                # Average over given time indices
                logger.debug("%d>%d?",t1,t0)
                logger.debug(t1>t0)
                if t1>t0:
                    data = pb.mean_profiles(data, t0, t1, h0, h1)
                # Save to plots, var[2] not needed?
                #plot_lines.append([var[0], var[1], var[2], var[4], data])
                plot_lines.append([var[0], var[1], var[4], data])
            else:
                plot_lines.append(None)
                functions.append([var[0], var[1], var[2], var[4], j])
        logger.debug('Iterate over functions')
        for func in functions:
            expression = func[2]
            logger.debug('Calculate %s', func[0])
            # Insert plot_lines entries into function string
            all_nan = True # Check if all values in expression are NAN
            for k in range(len(plot_lines)):
                if plot_lines[k] is not None and plot_lines[k][3] is not None:
                    if not np.all(np.isnan(plot_lines[k][3])):
                        all_nan = False
                    # Replace variables in expression with their respective entries from plot_lines and replace NAN values with zero
                    expression = expression.replace(plot_lines[k][0], 'np.where(np.isnan(plot_lines['+str(k)+'][3]),0,plot_lines['+str(k)+'][3])')
            # if all entries are NAN, return array containing only NANs
            if all_nan:
                data = np.full(nh, np.nan)
            else:
                try:
                    data = eval(expression)
                except Exception as e:
                    logger.error('Expression %s of line %s in plot %s could not be evaluated: %s', expression, plotLabels[i], func[0], e.message)
                    data = np.full(nh, np.nan)
            logger.debug(func[4])
            plot_lines[func[4]] = [func[0], func[1], func[3], data]
        plot_data.append(plot_lines)
    return dict(zip(plotLabels, plot_data))


def plot_default(plots, cf, data, h, centering):
    """
    bla
    TODO: pass startLevel etc. to plot routine (style file?)
    """
    logger.info('plot_default')
    # Initialise names
    out_dir, jpg_dir, plot_case_name, out_pdf = init_plotgen(plots, cf)
    # Create  output pdf
    pdf = PdfPages(out_pdf)
    # Loop over plots
    imageNames = []
    for i,plot_label in enumerate(plots.sortPlots):
        logger.info('Generating plot %s', plot_label)
        title, units = plots.plotNames[i]
        # pb.get_units() needs nc
        name = plot_case_name.format(plot=plot_label)
        imageNames.append(name)
        pb.plot_profiles(data[plot_label], h, units, cf.yLabel, title, os.path.join(jpg_dir,name), startLevel=0, lw=cf.lw, grid=False, centering=centering, pdf=pdf)
    pdf.close()
    
    # write html page
    logger.info("Write HTML page")
    index = out_dir + 'index.html'
    mode = 'Splotgen'
    logger.debug("HTML header template: "+cf.headerText)
    headerText = cf.headerText.format(type=plots.header)
    logger.debug("HTML header: "+headerText)
    ow.writeNavPage(out_dir, headerText)
    ow.writePlotsPage(out_dir, headerText, mode, imageNames)
    ow.writeIndex(index, mode)
    
    ## Add bookmarks to pdf output:
    logger.info("Add bookmarks to pdf output")
    writer = pdfw() # create pdf output class
    with open(out_pdf, 'rb') as inpdf: # open input file
        reader = pdfr(inpdf) # open pdf input class
        writer.appendPagesFromReader(reader) # copy all pages to output
        for i in range(writer.getNumPages()):
            writer.addBookmark(plots.sortPlots[i], i, parent=None) # add bookmark
        writer.setPageMode("/UseOutlines") # make pdf open bookmarks overview
        with open(os.path.join(out_dir,'tmp.pdf'),'wb') as tmp_pdf: # open output pdf file
            writer.write(tmp_pdf) # write pdf content to output file
    # file streams closed
    os.rename(os.path.join(out_dir,'tmp.pdf'), out_pdf) # rename newly created pdf


def plot_3d(plots, cf, data, h, prm_vars, fps=1, gif=False):
    """
    Generate horizontal output showing cloud outlines and wind speeds
    TODO: Special handling of RICO, because of huge grid -> split up grid in multiple subgrids?
    """
    logger.info('plot_3d')
    # Initialise names
    out_dir, jpg_dir, plot_case_name, out_pdf = init_plotgen(plots, cf)
    # Unpack arrays
    qn_3d = data['qn']
    u_3d = data['u']
    v_3d = data['v']
    w_3d = data['w']
    logger.debug(qn_3d.shape)
    logger.debug(u_3d.shape)
    logger.debug(v_3d.shape)
    logger.debug(w_3d.shape)
    # Define height level limits
    h_limits = np.array([0])
    logger.debug(h)
    for i in range(len(h)):
        h_limits = np.append(h_limits, 2*h[i]-h_limits[i])
    logger.debug(h_limits)
    h_extent = np.diff(h_limits)
    logger.debug(h_extent)
    # Get shape of arrays, here (z,x,y)
    dims = qn_3d.shape
    # Create meshgrid for 2d plot
    logger.info('Create meshgrid')
    grid_1d = np.arange(dims[1])*cf.dxy
    xgrid, ygrid = np.meshgrid(grid_1d, grid_1d)
    logger.info('Setup color maps')
    # Create normalization for projection onto colormap
    norm = mpl.colors.Normalize(vmin=w_3d.min(), vmax=w_3d.max())
    # Calculate cloud fraction for shading of height indicator
    cloud_frac_cmap = mpl.cm.get_cmap(plots.cloud_frac_cmap)
    cloud_frac = (qn_3d>0).mean(axis=(1,2))
    cloud_frac = cloud_frac/cloud_frac.max()
    colors = cloud_frac_cmap(cloud_frac)
    # Get quiver cmap
    quiver_cmap = mpl.cm.get_cmap(plots.quiver_cmap)
    # Define base figsize
    figsize = np.array(mpp.rcParams['figure.figsize'])*4
    logger.info('Setup mp4 output')
    # Initialize ffmpeg
    Writer = manim.writers['ffmpeg']
    writer = Writer(fps=fps, metadata={'artist':'Steffen Domke'})
    # Generate title
    title = plots.title_template.format(case = cf.case,
                                        x = dims[1],
                                        y = dims[2],
                                        z = dims[0],
                                        dx = prm_vars['dx'],
                                        dz = h_extent.mean(),
                                        t = (prm_vars['nsave3Dstart']+prm_vars['nsave3Dend'])*.5*prm_vars['dt'],
                                        h = '{h}',
                                        wt = '{wt}')
    # Define update function for mp4/gif creation
    def update_plot(i, wt):
        logger.debug('Generating frame %d',i)
        # Clear figures
        fig.clear()
        # Prepare canvas
        ax1, ax2 = fig.subplots(1, 2, gridspec_kw={'width_ratios':[50,1]})
        wt_title = wt.replace('_',' ')
        fig.suptitle(title.format(h=h[i], wt=wt_title), fontsize=20)
        #logger.debug("Set ylim to (%f, %f)", -prm_vars['dx'], dims[1]*prm_vars['dx'])
        ax1.set_ylim(-prm_vars['dx'], dims[1]*prm_vars['dx'])
        ax2.set_xlim(0,1)
        ax2.set_xticks([])
        ax2.set_ylim(0, h_limits[-1])
        ax1.set_title('Cloud cover and {} vectors'.format(wt_title), loc ='left', fontsize=20)
        ax2.set_title('height\nlevel')
        ax2.set_ylabel('Cloud fraction', fontsize=16)
        # Get data at height level i
        qn = qn_3d[i]
        u = u_3d[i]
        v = v_3d[i]
        w = w_3d[i]
        # Plot cloud contour fields
        cs = ax1.contourf(xgrid, ygrid, qn>0, levels=[0,.5,1], colors=['white','skyblue'])
        # Plot wind fields
        if 'total' in wt:
            q = ax1.quiver(xgrid, ygrid, u, v, w, angles='xy', minshaft=2, minlength=0, scale=quiver_scale_factor[cf.case], cmap=quiver_cmap, norm=norm)
        else:
            # Calculate deviation from horizontal mean
            u_dev = u - u.mean()
            v_dev = v - v.mean()
            q = ax1.quiver(xgrid, ygrid, u_dev, v_dev, w, angles='xy', minshaft=2, minlength=0, scale=quiver_scale_factor[cf.case]/5, cmap=quiver_cmap, norm=norm)
        
        # Add colorbar to ax1
        cb = fig.colorbar(q, ax=ax1, aspect=50)
        cb.ax.set_ylabel('Vertical wind velocity', fontsize=16)
        # Create wind vector field
        wind_field = np.stack((u,v))
        # Calculate mean horizontal wind strength
        mean_speed = np.linalg.norm(wind_field, axis=0).mean()
        # Calculate mean wind direction
        mean_wind = wind_field.mean(axis=(1,2))
        mean_dir = np.arccos(mean_wind[0]/np.linalg.norm(mean_wind))*180/np.pi
        if mean_wind[1]<0:
            mean_dir = 360-mean_dir
            # Generate quiver key
            ax1.quiverkey(q, .8, 1.01, mean_speed, 'Mean wind: {:.1f}'.format(mean_speed)+r'$\frac{m}{s}$', angle=mean_dir, labelsep=.5, labelpos='E', fontproperties={'size':20})
            # Fill second subplot
            # Generate Rectangles
            for idx in range(len(h)):
                ax2.add_patch(mpl.patches.Rectangle((0, h_limits[idx]), 1, h_extent[idx], color=colors[idx]))
                #ax2_dev.add_patch(mpl.patches.Rectangle((0, h_limits[idx]), 1, h_extent[idx], color=colors[idx]))
            # Add height indicator
            ax2.axhline(h[i], 0, 1, color='black')
            ax2.arrow(2.5, h[i], -.01, 0, color='black', head_width=h_extent.mean(), head_length=1).set_clip_on(False)
            # Generate output names
            #logger.debug(plot_case_name)
            plot_case_name_frame = plot_case_name.format(wt=wt, plot=int(h[i]))
            #logger.debug(os.path.join(out_dir, plot_case_name.format(plot=int(h[i]))))
            for el in ax1.get_xticklabels()+ax1.get_yticklabels():
                el.set_fontsize(20)
            # Reduce blank space between subplots (using fig.subplots_adjust, because tight_layout isnt working properly)
            #fig.tight_layout()
            fig.subplots_adjust(wspace=0)
            # Save as png
            #logger.debug(os.path.join(jpg_dir, plot_case_name_frame))
            fig.savefig(os.path.join(jpg_dir, plot_case_name_frame))
            # Save to pdf
            pdf.savefig(fig)
            #logger.debug("frame %d saved",i)
            fig.canvas.draw_idle()
    logger.info('Generate ouput')
    for wt in plots.wind_types:
        logger.info(wt)
        # Create figures for output
        fig = mpp.figure(figsize=figsize*figure_scale[cf.case])
        logger.info("Creating output pdf: %s", out_pdf)
        pdf = PdfPages(out_pdf.format(wt=wt))
        logger.info('Generate output')
        ani = manim.FuncAnimation(fig, update_plot, frames=len(h), fargs=(wt,), interval=1000./fps, repeat=False)
        logger.info('Save files')
        logger.debug('%d frames',len(h))
        if gif:
            ani.save(os.path.join(out_dir, plot_case_name.format(wt=wt, plot='mov'))+'.gif', writer='imagemagick', fps=30)
        else:
            ani.save(os.path.join(out_dir, plot_case_name.format(wt=wt, plot='mov'))+'.mp4', fps=30)
        pdf.close()
    

def plot_comparison(plots, cf, data_clubb, data_sam, h_clubb, h_sam, plot_old_clubb=False, data_old=None, h_old=None):
    """
    bla
    Added functionality to include old CLUBB data in plots. TODO: Deuglify code
    TODO: get units from clubb long name
    """
    logger.info('plot_comparison')
    # Initialise names
    out_dir, jpg_dir, plot_case_name, out_pdf = init_plotgen(plots, cf)
    # Create  output pdf
    pdf = PdfPages(out_pdf)
    # Loop over plots
    for i, plot_label in enumerate(plots.sortPlots):
        logger.info('Generating plot %s', plot_label)
        title, units = plots.plotNames[i]
        pb.plot_comparison(data_clubb[plot_label], data_sam[plot_label], h_clubb, h_sam, units, cf.yLabel, title, os.path.join(jpg_dir, plot_case_name.format(plot=plot_label)), startLevel=0, grid=False, pdf=pdf, plot_old_clubb=plot_old_clubb, data_old=data_old[plot_label], level_old=h_old)
    pdf.close()
    # TODO: Generate html
    
    ## Add bookmarks to pdf output:
    logger.info("Add bookmarks to pdf output")
    writer = pdfw() # create pdf output class
    with open(out_pdf, 'rb') as inpdf: # open input file
        reader = pdfr(inpdf) # open pdf input class
        writer.appendPagesFromReader(reader) # copy all pages to output
        for i in range(writer.getNumPages()):
            writer.addBookmark(plots.sortPlots[i], i, parent=None) # add bookmark
        writer.setPageMode("/UseOutlines") # make pdf open bookmarks overview
        with open(os.path.join(out_dir,'tmp.pdf'),'wb') as tmp_pdf: # open output pdf file
            writer.write(tmp_pdf) # write pdf content to output file
    # file streams closed
    os.rename(os.path.join(out_dir,'tmp.pdf'), out_pdf) # rename newly created pdf



#-------------------------------------------------------------------------------
#   MAIN PLOT FUNCTIONS
#-------------------------------------------------------------------------------

def plotgen_default(plots, cf):
    """
    bla
    This plotting routine plots height profiles of data from either the SAM or CLUBB simulations
    """
    logger.info("plotgen_default")
    ncs = load_nc(plots, cf)
    if 'clubb' in plots.name:
        # Separate CLUBB nc files
        if plots.sortPlots_zm and plots.sortPlots_zt:
            nc_zm, nc_zt = ncs
            nc = nc_zm
        elif plots.sortPlots_zm:
            nc_zm, = ncs
            nc = nc_zm
            nc_zt = None
        else:
            nc_zt, = ncs
            nc = nc_zt
            nc_zm = None
    else:
        nc_sam, = ncs
        nc = nc_sam
    logger.info('Fetching dimension variables')
    t = get_t_dim(nc)
    # Distinguish between CLUBB and SAM, as z dimension variables have different names
    if 'clubb' in plots.name:
        h = get_h_dim(nc, model='clubb')
    else:
        h = get_h_dim(nc_sam, model='sam')
    logger.info('Find nearest levels to case setup')
    idx_h0 = (np.abs(h - cf.startHeight)).argmin()
    idx_h1 = (np.abs(h - cf.endHeight)).argmin()+1
    # Times for CLUBB are given in seconds, for SAM in minutes
    if 'clubb' in plots.name:
        idx_t0 = (np.abs(t - cf.startTime*60)).argmin()
        idx_t1 = (np.abs(t - cf.endTime*60)).argmin()+1
    else:
        idx_t0 = (np.abs(t - cf.startTime)).argmin()
        idx_t1 = (np.abs(t - cf.endTime)).argmin()+1
    h = h[idx_h0:idx_h1]
    if 'clubb' in plots.name:
        # Initialize data variables as empty dicts
        data_zm = {}
        data_zt = {}
        if nc_zm:
            logger.info("Fetching clubb_zm data")
            data_zm = get_all_variables(nc_zm, plots.lines_zm, plots.sortPlots_zm, len(h), len(t), idx_t0, idx_t1, idx_h0, idx_h1, filler=plots.filler)
        if nc_zt:
            logger.info("Fetching clubb_zt data")
            data_zt = get_all_variables(nc_zt, plots.lines_zt, plots.sortPlots_zt, len(h), len(t), idx_t0, idx_t1, idx_h0, idx_h1, filler=plots.filler)
        # Combine both dictionaries
        data = data_zm.copy()
        data.update(data_zt)
    else:
        logger.info("Fetching sam data")
        data = get_all_variables(nc_sam, plots.lines, plots.sortPlots, len(h), len(t), idx_t0, idx_t1, idx_h0, idx_h1, filler=plots.filler)
    logger.info("Create plots")
    # Center plots if plotting budgets
    logger.debug("Centering=%s", str('budget' in plots.name))
    plot_default(plots, cf, data, h, centering='budget' in plots.name)
    #logger.debug(mpp.rcParams['text.usetex'])


def plotgen_3d(plots, cf):
    """
    bla
    TODO: specify mp4 or gif outout type
    """
    logger.info("plotgen_3d")
    nc, = load_nc(plots, cf)
    logger.info('Fetching information from prm file')
    prm_vars = get_values_from_prm(cf)
    logger.info('Fetching dimension variables')
    # t not necessarily needed for momentary plot TODO
    #t = get_t_dim(nc, create=True, dt=prm_vars['dt'])
    h = get_h_dim(nc, model='sam', create=True, cf=cf)
    idx_h0 = (np.abs(h - cf.startHeight)).argmin()
    idx_h1 = (np.abs(h - cf.endHeight)).argmin()+1
    h = h[idx_h0:idx_h1]
    logger.info("Fetching sam_3d data")
    data = get_all_variables(nc, plots.lines, plots.sortPlots, len(h), 1, 0, 0, idx_h0, idx_h1, filler=plots.filler)
    data = {key:value[0][-1]  for (key,value) in data.iteritems()}
    logger.debug(data)
    logger.info("Create plots")
    plot_3d(plots, cf, data, h, prm_vars)
    
    
def plotgen_comparison(plots, cf):
    """
    bla
    Added fugly code to include old CLUBB data in comparison plots. TODO: Find better way to implement this.
    """
    logger.info("plotgen_comparison")
    old_clubb = None
    for n in range(ntrials):
        user_input = raw_input('Include old CLUBB data?(y/n) -> ').lower()
        if user_input in affirmatives:
            old_clubb = True
            break
        elif user_input in negatives:
            old_clubb = False
            break
        else:
            logger.info('Trial #%d/%d. Invalid input. Please type y to include old CLUBB data and n otherwise.', n+1, ntrials)
    if old_clubb is None:
        logger.error('Too many tries. Ending program...')
        sys.exit()
    if old_clubb:
        logger.info('Using old CLUBB data.')
    ncs = load_nc(plots, cf, old_clubb)
    # Separate nc files
    # Get nc files for old CLUBB data
    nc_zm_old = None
    nc_zt_old = None
    if old_clubb:
        if plots.sortPlots_zt:
            nc_zt_old = ncs.pop()
        if plots.sortPlots_zm:
            nc_zm_old = ncs.pop()
            nc_clubb_old = nc_zm_old
        else:
            nc_clubb_old = nc_zt_old
    if plots.sortPlots_zm and plots.sortPlots_zt:
        nc_zm, nc_zt, nc_sam = ncs
        nc_clubb = nc_zm
    elif plots.sortPlots_zm:
        nc_zm, nc_sam = ncs
        nc_clubb = nc_zm
        nc_zt = None
    else:
        nc_zt, nc_sam = ncs
        nc_clubb = nc_zt
        nc_zm = None
    logger.info('Fetching dimension variables')
    logger.info('New CLUBB')
    t_clubb = get_t_dim(nc_clubb)
    h_clubb = get_h_dim(nc_clubb, model='clubb')
    idx_h0_clubb = (np.abs(h_clubb - cf.startHeight)).argmin()
    idx_h1_clubb = (np.abs(h_clubb - cf.endHeight)).argmin()+1
    # Times for CLUBB are given in seconds, for SAM in minutes
    idx_t0_clubb = (np.abs(t_clubb - cf.startTime*60)).argmin()
    idx_t1_clubb = (np.abs(t_clubb - cf.endTime*60)).argmin()+1
    h_clubb = h_clubb[idx_h0_clubb:idx_h1_clubb]
    logger.info('SAM')
    t_sam = get_t_dim(nc_sam)
    h_sam = get_h_dim(nc_sam, model='sam')
    idx_h0_sam = (np.abs(h_sam - cf.startHeight)).argmin()
    idx_h1_sam = (np.abs(h_sam - cf.endHeight)).argmin()+1
    idx_t0_sam = (np.abs(t_sam - cf.startTime)).argmin()
    idx_t1_sam = (np.abs(t_sam - cf.endTime)).argmin()+1
    h_sam = h_sam[idx_h0_sam:idx_h1_sam]
    if old_clubb:
        t_clubb_old = get_t_dim(nc_clubb_old)
        h_clubb_old = get_h_dim(nc_clubb_old, model='clubb')
        idx_h0_clubb_old = (np.abs(h_clubb_old - cf.startHeight)).argmin()
        idx_h1_clubb_old = (np.abs(h_clubb_old - cf.endHeight)).argmin()+1
        # Times for CLUBB are given in seconds, for SAM in minutes
        idx_t0_clubb_old = (np.abs(t_clubb_old - cf.startTime*60)).argmin()
        idx_t1_clubb_old = (np.abs(t_clubb_old - cf.endTime*60)).argmin()+1
        h_clubb_old = h_clubb_old[idx_h0_clubb_old:idx_h1_clubb_old]
        logger.debug("h0_clubb_old = %d",idx_h0_clubb)
        logger.debug("h1_clubb_old = %d",idx_h1_clubb)
        logger.debug("t0_clubb_old = %d",idx_t0_clubb)
        logger.debug("t1_clubb_old = %d",idx_t1_clubb)
        logger.debug(h_clubb_old.shape)
        logger.debug(h_clubb_old[idx_h0_clubb_old:idx_h1_clubb_old])
    logger.debug("h0_clubb = %d",idx_h0_clubb)
    logger.debug("h1_clubb = %d",idx_h1_clubb)
    logger.debug("t0_clubb = %d",idx_t0_clubb)
    logger.debug("t1_clubb = %d",idx_t1_clubb)
    logger.debug("h0_sam = %d",idx_h0_sam)
    logger.debug("h1_sam = %d",idx_h1_sam)
    logger.debug("t0_sam = %d",idx_t0_sam)
    logger.debug("t1_sam = %d",idx_t1_sam)
    logger.debug(h_clubb.shape)
    logger.debug(h_clubb[idx_h0_clubb:idx_h1_clubb])
    logger.debug(h_sam.shape)
    logger.debug(h_sam[idx_h0_sam:idx_h1_sam])
    # TODO: distinguish between zm and zt data!!
    logger.info("Fetching SAM data")
    data_sam = get_all_variables(nc_sam, plots.lines_sam, plots.sortPlots, len(h_sam), len(t_sam), idx_t0_sam, idx_t1_sam, idx_h0_sam, idx_h1_sam, filler=plots.filler)
    #logger.debug(data_sam)
    logger.info("Fetching CLUBB data")
    data_zm = None
    data_zt = None
    if plots.sortPlots_zm:
        logger.info("Fetching clubb_zm data")
        data_zm = get_all_variables(nc_zm, plots.lines_zm, plots.sortPlots_zm, len(h_clubb), len(t_clubb), idx_t0_clubb, idx_t1_clubb, idx_h0_clubb, idx_h1_clubb, filler=plots.filler)
        #logger.debug(data_clubb_zm)
    if plots.sortPlots_zt:
        logger.info("Fetching clubb_zt data")
        data_zt = get_all_variables(nc_zt, plots.lines_zt, plots.sortPlots_zt, len(h_clubb), len(t_clubb), idx_t0_clubb, idx_t1_clubb, idx_h0_clubb, idx_h1_clubb, filler=plots.filler)
        #logger.debug(data_clubb_zt)
    # combine both dictionaries
    data_clubb = data_zm.copy()
    data_clubb.update(data_zt)
    logger.info("Fetching old CLUBB data")
    if old_clubb:
        data_zm_old = None
        data_zt_old = None
        if plots.sortPlots_zm:
            logger.info("Fetching old clubb_zm data")
            data_zm_old = get_all_variables(nc_zm_old, plots.lines_zm, plots.sortPlots_zm, len(h_clubb_old), len(t_clubb_old), idx_t0_clubb_old, idx_t1_clubb_old, idx_h0_clubb_old, idx_h1_clubb_old, filler=plots.filler)
            logger.debug(data_zm_old)
        if plots.sortPlots_zt:
            logger.info("Fetching old clubb_zt data")
            data_zt_old = get_all_variables(nc_zt_old, plots.lines_zt, plots.sortPlots_zt, len(h_clubb_old), len(t_clubb_old), idx_t0_clubb_old, idx_t1_clubb_old, idx_h0_clubb_old, idx_h1_clubb_old, filler=plots.filler)
            logger.debug(data_zt_old)
        # combine both dictionaries
        data_clubb_old = data_zm_old.copy()
        data_clubb_old.update(data_zt_old)
        logger.debug(data_clubb_old)
    logger.info("Create plots")
    plot_comparison(plots, cf, data_clubb, data_sam, h_clubb, h_sam, plot_old_clubb=old_clubb, data_old=data_clubb_old, h_old=h_clubb_old)


# Compile plot functions in dict
plot_dict = {
    "budget" : plotgen_default,
    "standalone" : plotgen_default,
    "3d" : plotgen_3d,
    "comparison" : plotgen_comparison,
    }