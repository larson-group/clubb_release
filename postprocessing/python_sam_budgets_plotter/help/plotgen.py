
# General imports
import os
import sys
import re
import logging
from datetime import datetime as dt
import numpy as np
from scipy.ndimage.morphology import binary_dilation as bin_dilation
import matplotlib as mpl
from matplotlib import pyplot as mpp
import matplotlib.animation as manim
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import ScalarFormatter as stick
from mpl_toolkits.axes_grid1 import make_axes_locatable
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
                else:
                    data = data[h0:h1]
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
        pb.plot_profiles(data[plot_label], h, units, cf.yLabel, title, os.path.join(jpg_dir,name), startLevel=plots.startLevel, lw=cf.lw, grid=False, centering=centering, pdf=pdf)
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


def plot_3d(plots, cf, data, h, h_limits, h_extent, prm_vars, fps=2, gif=False):
    """
    Generate horizontal output showing cloud outlines and wind speeds
    TODO: Special handling of RICO, because of huge grid -> split up grid in multiple subgrids?
    """
    logger.info('plot_3d')
    # Initialise names
    out_dir, jpg_dir, plot_case_name, out_pdf = init_plotgen(plots, cf)
    
    ### Get data
    logger.info('Unpack data and perform calculations')
    ## Unpack arrays
    qn_3d = data['qn']
    u_3d = data['u']
    v_3d = data['v']
    w_3d = data['w']
    ## Calculate needed values
    # Create cloud mask
    cld_mask = qn_3d>cld_lim
    # Create dilated cloud mask with cross structure element
    # Create structure element
    dot = np.zeros((3,3))
    dot[1,1]=1
    mid = np.ones((3,3))
    for i in range(4):
        mid[i/2*2,i%2*2]=0
    struct = np.stack((dot,mid,dot))
    # Apply dilation
    cld_mask_dilated = bin_dilation(cld_mask, struct)
    # Profile of mean wind component
    um = u_3d.mean(axis=(1,2))
    vm = v_3d.mean(axis=(1,2))
    wm = w_3d.mean(axis=(1,2))
    # 3d field of u'w' and v'w'
    up = u_3d-um.reshape((um.size,1,1))
    vp = v_3d-vm.reshape((vm.size,1,1))
    wp = w_3d-wm.reshape((wm.size,1,1))
    uw = up*wp
    vw = vp*wp
    # Profile of horizontal wind strength
    wind_strength = np.linalg.norm(np.stack((um,vm)), axis=0)
    # Profile of wind direction
    mean_dir = np.arccos(um/wind_strength)*180/np.pi
    mean_dir = np.where(vm<0, 360-mean_dir, mean_dir)
    # DEBUGGING
    logger.debug("qn_3d.shape=%s",str(qn_3d.shape))
    logger.debug("u_3d.shape=%s",str(u_3d.shape))
    logger.debug("v_3d.shape=%s",str(v_3d.shape))
    logger.debug("w_3d.shape=%s",str(w_3d.shape))
    logger.debug("um.shape=%s",str(um.shape))
    logger.debug("vm.shape=%s",str(vm.shape))
    logger.debug("wm.shape=%s",str(wm.shape))
    logger.debug("up.shape=%s",str(up.shape))
    logger.debug("vp.shape=%s",str(vp.shape))
    logger.debug("wp.shape=%s",str(wp.shape))
    logger.debug("uw.shape=%s",str(uw.shape))
    logger.debug("vw.shape=%s",str(vw.shape))
    logger.debug("wind_strength.shape=%s",str(wind_strength.shape))
    logger.debug("mean_dir.shape=%s",str(mean_dir.shape))
    logger.debug("uw.shape=%s",str(uw.shape))
    ## Define height level limits (DEPRECATED: moved to plotgen_3d and introduced as parameter, as limits not reconstructable from sliced h)
    #h_limits = np.array([0])
    #logger.debug(h)
    #for i in range(len(h)):
        #h_limits = np.append(h_limits, 2*h[i]-h_limits[i])
    #logger.debug(h_limits)
    #h_extent = np.diff(h_limits)
    #logger.debug(h_extent)
    # Create arrays to save conditional u'w'
    uw_cld = np.zeros(um.size)
    vw_cld = np.zeros(vm.size)
    uw_dilated = np.zeros(um.size)
    vw_dilated = np.zeros(vm.size)
    
    ### Get colormaps
    # TODO: what if uw.min*uw.max>0?
    # TODO: for uw and vw, what to do with outliers???
    logger.info('Setup color maps')
    ## Get cloud fraction color map
    # Calculate cloud fraction for shading of height indicator
    cloud_frac_cmap = mpl.cm.get_cmap(plots.cloud_frac_cmap)
    cloud_frac = cld_mask.mean(axis=(1,2))
    cloud_frac = cloud_frac/cloud_frac.max()
    # Generate colors from color map
    cloud_colors = cloud_frac_cmap(cloud_frac)
    ## Generate quiver colormap
    # Define segment colors:
    #red = mpl.colors.to_rgb('xkcd:red')
    red = (.9,0,0)
    #green = mpl.colors.to_rgb('tab:green')
    green = (0,.75,0)
    cld_col = 'skyblue'
    dilated_col = 'mediumblue'
    grey_tone = 1       # grey value at 0
    ### Get colormaps and normalizations for projection onto colormap
    ## UW
    # Get limits for uw data in order to generate colormap
    uwlim = uw.min(), uw.max()
    logger.debug("uwlim=(%f,%f)",uwlim[0], uwlim[1])
    logger.debug('uwlim[0]=(%f,%f)',uw[0].min(), uw[0].max())
    norm_uw = mpl.colors.Normalize(vmin=uwlim[0], vmax=uwlim[1])
    # Calculate position of 0 in normalized data (projected onto [0;1])
    zero_uw = (0-uwlim[0])/(uwlim[1]-uwlim[0])
    # Define color map segments
    cdict_uw = {"red":((0,red[0],red[0]),(zero_uw,grey_tone,grey_tone),(1,green[0],green[0])),
                "green":((0,red[1],red[1]),(zero_uw,grey_tone,grey_tone),(1,green[1],green[1])),
                "blue":((0,red[2],red[2]),(zero_uw,grey_tone,grey_tone),(1,green[2],green[2])),
                #"alpha":((0,1,1),(zero,0,0),(1,1,1))
                }
    # Create color map
    quiver_cmap_uw = mpl.colors.LinearSegmentedColormap('quiver_cmap_uw', cdict_uw, N=1000)
    
    
    ## VW
    # Get limits for vw data in order to generate colormap
    vwlim = vw.min(), vw.max()
    logger.debug("vwlim=(%f,%f)",vwlim[0], vwlim[1])
    logger.debug('vwlim[0]=(%f,%f)',vw[0].min(), vw[0].max())
    norm_vw = mpl.colors.Normalize(vmin=vwlim[0], vmax=vwlim[1])
    # Calculate position of 0 in normalized data (projected onto [0;1])
    zero_vw = (0-vwlim[0])/(vwlim[1]-vwlim[0])
    # Define color map segments
    cdict_vw = {"red":((0,red[0],red[0]),(zero_vw,grey_tone,grey_tone),(1,green[0],green[0])),
                "green":((0,red[1],red[1]),(zero_vw,grey_tone,grey_tone),(1,green[1],green[1])),
                "blue":((0,red[2],red[2]),(zero_vw,grey_tone,grey_tone),(1,green[2],green[2])),
                #"alpha":((0,1,1),(zero,0,0),(1,1,1))
                }
    # Create color map
    quiver_cmap_vw = mpl.colors.LinearSegmentedColormap('quiver_cmap_vw', cdict_vw, N=1000)
    
    ## W
    # Get limits of w data in order to generate colormap
    wlim = w_3d.min(),w_3d.max()
    norm_w = mpl.colors.Normalize(vmin=wlim[0], vmax=wlim[1])
    # Calculate position of 0 in normalized data (projected onto [0;1])
    zero_w = (0-wlim[0])/(wlim[1]-wlim[0])
    # Define color map segments
    cdict_w = {"red":((0,red[0],red[0]),(zero_w,grey_tone,grey_tone),(1,green[0],green[0])),
             "green":((0,red[1],red[1]),(zero_w,grey_tone,grey_tone),(1,green[1],green[1])),
             "blue":((0,red[2],red[2]),(zero_w,grey_tone,grey_tone),(1,green[2],green[2])),
             #"alpha":((0,1,1),(zero,0,0),(1,1,1))
             }
    # Create color map
    quiver_cmap_w = mpl.colors.LinearSegmentedColormap('quiver_cmap_w',cdict_w, N=1000)
    
    ## Plotting preparations
    # Get shape of arrays, here (z,x,y)
    dims = qn_3d.shape
    # Create meshgrid for 2d plot
    logger.info('Create meshgrid')
    grid_1d = np.arange(dims[1])
    xgrid, ygrid = np.meshgrid(grid_1d, grid_1d)
    # Define base figsize (DEPRECATED: moved to plot_defs.py)
    #figsize = np.array(mpp.rcParams['figure.figsize'])*4
    #figsize = (32.5,20)
    logger.info('Setup mp4 output')
    # Initialize ffmpeg
    Writer = manim.writers['ffmpeg']
    writer = Writer(fps=fps, metadata={'artist':'Steffen Domke'})
    # Generate title (TODO: denote time interval in case file?)
    title = plots.title_template.format(
        case = cf.case,
        x = dims[1],
        y = dims[2],
        z = dims[0],
        dx = prm_vars['dx'],
        dz = h_extent.mean(),
        t = (cf.time_3d*cf.dt)/60.)
    # Define update function for mp4/gif creation
    def update_plot(framenumber, wt):
        logger.debug('Generating frame %d', framenumber)
        #logger.debug('wt: %s', wt)
        
        # Fill conditional u'w' arrays
        if cld_mask[framenumber].sum()>0:
            uw_cld[framenumber] = uw[framenumber, cld_mask[framenumber]].mean()
            vw_cld[framenumber] = vw[framenumber, cld_mask[framenumber]].mean()
        if cld_mask_dilated[framenumber].sum()>0:
            uw_dilated[framenumber] = uw[framenumber, cld_mask_dilated[framenumber]].mean()
            vw_dilated[framenumber] = vw[framenumber, cld_mask_dilated[framenumber]].mean()
        
        # Clear figures
        fig.clear()
        # Prepare canvas
        ax1, ax2 = fig.subplots(1, 2, gridspec_kw={'width_ratios':width_ratio})
        # Get data at height level framenumber
        #qn = qn_3d[framenumber]
        u = u_3d[framenumber]
        v = v_3d[framenumber]
        #w = w_3d[framenumber]
        # Print vertical wind component array as background image
        if 'uw' in wt:
            pc = ax1.imshow(uw[framenumber], cmap=quiver_cmap_uw, norm=norm_uw, interpolation=cloud_interpolation, origin='lower')
        elif 'vw' in wt:
            pc = ax1.imshow(vw[framenumber], cmap=quiver_cmap_vw, norm=norm_vw, interpolation=cloud_interpolation, origin='lower')
        else:
            pc = ax1.imshow(w_3d[framenumber], cmap=quiver_cmap_w, norm=norm_w, interpolation=cloud_interpolation, origin='lower')
        # Plot cloud contour fields
        cs = ax1.contourf(xgrid, ygrid, cld_mask[framenumber], levels=[.5,1], colors=[cld_col], extend='neither', alpha=.5)
        # Plot contour of halo
        ax1.contour(xgrid, ygrid, cld_mask_dilated[framenumber], levels=[.5,1], colors=[dilated_col], extend='neither', alpha=.5)
        # Plot wind fields
        if 'total' in wt:
            #q = ax1.quiver(xgrid, ygrid, u, v, w, angles='xy', minshaft=2, minlength=0, scale=quiver_scale_factor[cf.case], cmap=quiver_cmap, norm=norm)
            q = ax1.quiver(xgrid, ygrid, u, v, angles='xy', minshaft=2, minlength=0, scale=quiver_scale_factor[cf.case], pivot='mid')
        else:
            # Calculate deviation from horizontal mean
            #u_dev = u - um[framenumber]
            #v_dev = v - v[framenumber]
            #q = ax1.quiver(xgrid, ygrid, u_dev, v_dev, w, angles='xy', minshaft=2, minlength=0, scale=quiver_scale_factor[cf.case]/5, cmap=quiver_cmap, norm=norm)
            q = ax1.quiver(xgrid, ygrid, u-um[framenumber], v-vm[framenumber], angles='xy', minshaft=2, minlength=0, scale=quiver_scale_factor[cf.case]/5, pivot='mid')
        
        # Add colorbar to ax1
        divider = make_axes_locatable(ax1)
        cax = divider.append_axes('right', size='5%',pad=.1)
        cb = fig.colorbar(pc, cax=cax)
        if 'uw' in wt:
            cb.ax.set_ylabel(r"u'w' $\mathrm{\left[\frac{m^2}{s^2}\right]}$", fontsize=fontsizes['labels'])
        elif 'vw' in wt:
            cb.ax.set_ylabel(r"v'w' $\mathrm{\left[\frac{m^2}{s^2}\right]}$", fontsize=fontsizes['labels'])
        else:
            cb.ax.set_ylabel(r'w $\mathrm{\left[\frac{m}{s}\right]}$', fontsize=fontsizes['labels'])
        cb.ax.tick_params(labelsize=fontsizes['labels'])
        # Generate quiver key
        ax1.quiverkey(q, .8, 1.01, wind_strength[framenumber], 'Mean wind: {:.1f}'.format(wind_strength[framenumber])+r'$\frac{m}{s}$', angle=mean_dir[framenumber], labelsep=.2*np.abs(um).max(), labelpos='E', fontproperties={'size':20})
        # Fill second subplot
        for line in plots.sortPlots_std:
            ax2.plot(data[line][1], h, label=data[line][0], color=col_3d[line], ls='-', lw=5)
        # Generate Rectangles
        xlims = ax2.get_xlim()
        for idx in range(len(h)):
            ax2.add_patch(mpl.patches.Rectangle((xlims[0], h_limits[idx]), xlims[1]-xlims[0], h_extent[idx], color=cloud_colors[idx]))
        # Format plot
        wt_title = wt.replace('_',' ')
        fig.suptitle(title.format(h=h[framenumber], wt=wt_title), fontsize=fontsizes['title'])
        #logger.debug("Set ylim to (%f, %f)", -prm_vars['dx'], dims[1]*prm_vars['dx'])
        ax1.set_xlim(-.6, dims[1]-.4)
        ax1.set_ylim(-.6, dims[2]-.4)
        ax2.set_ylim(0, h_limits[-1])
        ax2.yaxis.tick_right()
        ax2.yaxis.set_label_position('right')
        ax1.set_xlabel(r'Eastward grid dimension $\mathrm{{\left[{:.0f} m\right]}}$'.format(prm_vars['dx']), fontsize=fontsizes['labels'])
        ax1.set_ylabel(r'Northward grid dimension $\mathrm{{\left[{:.0f} m\right]}}$'.format(prm_vars['dx']), fontsize=fontsizes['labels'])
        ax1.set_title('Cloud cover and {}'.format(wt_title.split('+')[0]), loc ='left', fontsize=fontsizes['title'])
        ax2.set_title('Profile of momentum fluxes', fontsize=fontsizes['title'])
        ax2.set_ylabel(r'Cloud fraction / Height $\mathrm{\left[m\right]}$', fontsize=fontsizes['labels'])
        ax2.set_xlabel(r"Vertical momentum fluxes $\mathrm{\left[\frac{m^2}{s^2}\right]}$", fontsize=fontsizes['labels'])
        # Generate output names
        #logger.debug(plot_case_name)
        plot_case_name_frame = plot_case_name.format(wt=wt, plot=int(h[framenumber]))
        #logger.debug(os.path.join(out_dir, plot_case_name.format(plot=int(h[framenumber]))))
        ## Modify ticks
        # Set fontsize
        ax1.tick_params(axis='both', labelsize=fontsizes['ticks'])
        ax2.tick_params(axis='both', labelsize=fontsizes['ticks'])
        # Change labels
        #labels = []
        #logger.debug('xticks old: %s',str(list(ax1.get_xticklabels())))
        #for el in ax1.get_xticks():
            #try:
                #labels.append(int(el*prm_vars['dx']))
            #except ValueError as ve:
                #labels.append('')
        #ax1.set_xticklabels(labels)
        ##logger.debug('xticks new: %s',str(list(ax1.get_xticklabels())))
        #labels = []
        ##logger.debug('yticks old: %s',str(list(ax1.get_yticklabels())))
        #for el in ax1.get_yticks():
            #try:
                #labels.append(int(el*prm_vars['dx']))
            #except ValueError as ve:
                #labels.append('')
        #ax1.set_yticklabels(labels)
        ax2.legend(loc=legend_pos_3d, prop={'size': fontsizes['legend']})
        ax2.axvline(x=0, color='k', ls='--', alpha=.5)
        # Add height indicator
        ax2.axhline(h[framenumber], 0, 1, color='black')
        ax2.arrow(xlims[0]-(xlims[1]-xlims[0])/10, h[framenumber], (xlims[1]-xlims[0])/10, 0, length_includes_head=True, color='black', head_width=h_limits[-1]/50., head_length=(xlims[1]-xlims[0])/10).set_clip_on(False)
        # Set tick label format to scalar formatter with length sensitive format (switch to scientific format when a set order of magnitude is reached)
        ticks = stick(useMathText=True)
        ticks.set_powerlimits((-pow_lim,pow_lim))
        ax2.xaxis.get_offset_text().set_size(20)
        #logger.debug('yticks new: %s',str(list(ax1.get_yticklabels())))
        # Reduce blank space between subplots (using fig.subplots_adjust, because tight_layout isnt working properly)
        #fig.tight_layout()
        fig.subplots_adjust(wspace=0, hspace=0)
        # Save as png
        #logger.debug(os.path.join(jpg_dir, plot_case_name_frame))
        fig.savefig(os.path.join(jpg_dir, plot_case_name_frame))
        # Save to pdf
        pdf.savefig(fig)
        #logger.debug("frame %d saved",framenumber)
        fig.canvas.draw_idle()
    logger.info('Generate ouput')
    for wt in plots.wind_types:
        logger.info(wt)
        # Create figures for output
        fig = mpp.figure(figsize=np.array(figsize)*figure_scale[cf.case])
        logger.info("Creating output pdf: %s", out_pdf)
        pdf = PdfPages(out_pdf.format(wt=wt))
        logger.info('Generate output')
        ani = manim.FuncAnimation(fig, update_plot, frames=len(h), fargs=(wt,), interval=1000./fps, repeat=False)
        logger.info('Save files')
        logger.debug('%d frames',len(h))
        if gif:
            ani.save(os.path.join(out_dir, plot_case_name.format(wt=wt, plot='mov'))+'.gif', writer='imagemagick', fps=fps)
        else:
            ani.save(os.path.join(out_dir, plot_case_name.format(wt=wt, plot='mov'))+'.mp4', fps=fps)
        pdf.close()
    # Create plots for conditional u'w':
    logger.debug(plot_case_name)
    uw_pdf = PdfPages(os.path.join(out_dir, out_pdf.format(wt='conditional_uw_profiles')))
    pb.plot_profiles([["In-cloud mean of u'w'",True,0,uw_cld]], h, r"Cloud Conditional $\mathrm{\overline{u'w'}\ \left[\frac{m^2}{s^2}\right]}$", cf.yLabel, "In-cloud mean of u'w'", os.path.join(jpg_dir,'uw_cld'), startLevel=0, lw=cf.lw, grid=False, centering=False, pdf=uw_pdf)
    
    pb.plot_profiles([["Extended in-cloud mean of u'w'",True,0,uw_dilated]], h, r"Dilated Cloud Conditional $\mathrm{\overline{u'w'}\ \left[\frac{m^2}{s^2}\right]}$", cf.yLabel, "Extended in-cloud mean of u'w'", os.path.join(jpg_dir,'uw_dilated'), startLevel=0, lw=cf.lw, grid=False, centering=False, pdf=uw_pdf)
    
    pb.plot_profiles([["In-cloud mean of v'w'",True,0,vw_cld]], h, r"Cloud Conditional $\mathrm{\overline{v'w'}\ \left[\frac{m^2}{s^2}\right]}$", cf.yLabel, "In-cloud mean of v'w'", os.path.join(jpg_dir,'vw_cld'), startLevel=0, lw=cf.lw, grid=False, centering=False, pdf=uw_pdf)
    
    pb.plot_profiles([["Extended in-cloud mean of v'w'",True,0,vw_dilated]], h, r"Dilated Cloud Conditional $\mathrm{\overline{v'w'}\ \left[\frac{m^2}{s^2}\right]}$", cf.yLabel, "Extended in-cloud mean of v'w'", os.path.join(jpg_dir,'vw_dilated'), startLevel=0, lw=cf.lw, grid=False, centering=False, pdf=uw_pdf)
    uw_pdf.close()
    
    

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
        if plot_old_clubb:
            pb.plot_comparison(data_clubb[plot_label], data_sam[plot_label], h_clubb, h_sam, units, cf.yLabel, title, os.path.join(jpg_dir, plot_case_name.format(plot=plot_label)), startLevel=0, grid=False, pdf=pdf, plot_old_clubb=plot_old_clubb, data_old=data_old[plot_label], level_old=h_old)
        else:
            pb.plot_comparison(data_clubb[plot_label], data_sam[plot_label], h_clubb, h_sam, units, cf.yLabel, title, os.path.join(jpg_dir, plot_case_name.format(plot=plot_label)), startLevel=plots.startLevel, grid=False, pdf=pdf)
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
    TODO: No data from std sam nc file needed
    """
    logger.info("plotgen_3d")
    gif = None
    for n in range(ntrials):
        user_input = raw_input('Movie output as .mp4 (1) or .gif(2)? -> ')
        try:
            if int(user_input)==2:
                gif = True
                break
            elif int(user_input)==1:
                gif = False
                break
            else:
                logger.warning('Invalid input. Choice of valid input options: 1 (.mp4), 2 (.gif)')
                continue
        except ValueError as err:
            logger.warning('Invalid input. Choice of valid input options: 1 (.mp4), 2 (.gif)')
            continue
    if gif is None:
        logger.error('Too many invalid trials. Exiting...')
        sys.exit()
    nc_std,nc_3d = load_nc(plots, cf)
    logger.info('Fetching information from prm file')
    prm_vars = get_values_from_prm(cf)
    logger.info('Fetching dimension variables')
    # t not necessarily needed for momentary plot TODO
    t = get_t_dim(nc_std, create=False)
    idx_t0 = (np.abs(t - cf.startTime)).argmin()
    idx_t1 = (np.abs(t - cf.endHeight)).argmin()+1
    h = get_h_dim(nc_std, model='sam', create=False, cf=cf)
    # Calculate height level limits
    h_limits = np.array([0])
    logger.debug(h)
    for i in range(len(h)):
        h_limits = np.append(h_limits, 2*h[i]-h_limits[i])
    # Calculate height indices based on parameters in case file
    idx_h0 = (np.abs(h - cf.startHeight)).argmin()
    idx_h1 = (np.abs(h - cf.endHeight)).argmin()+1
    logger.debug('Height indices: [%d,%d]',idx_h0, idx_h1)
    # Slice arrays
    h = h[idx_h0:idx_h1]
    h_limits = h_limits[idx_h0:idx_h1+1]
    logger.debug(h_limits)
    # Calculate height level extents
    h_extent = np.diff(h_limits)
    logger.debug(h_extent)
    logger.info("Fetching sam data")
    data_std = get_all_variables(nc_std, plots.lines_std, plots.sortPlots_std, len(h), len(t), idx_t0, idx_t1, idx_h0, idx_h1, filler=plots.filler)
    logger.info("Fetching sam_3d data")
    data_3d = get_all_variables(nc_3d, plots.lines_3d, plots.sortPlots_3d, len(h), 1, 0, 0, idx_h0, idx_h1, filler=plots.filler)
    data = {key:value[0][-1]  for (key,value) in data_3d.iteritems()}
    data.update({key:[value[0][0], value[0][-1]] for (key,value) in data_std.iteritems()})
    logger.debug("3d data: %s", str(data))
    logger.info("Create plots")
    plot_3d(plots, cf, data, h, h_limits, h_extent, prm_vars, gif=gif)
    
    
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
    h_clubb_old = None
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
    data_clubb_old = None
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