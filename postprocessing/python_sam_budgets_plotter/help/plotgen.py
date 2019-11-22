# TODO: Hide differences between SAM and CLUBB in functions
#       -> Generate list of nc files from CLUBB and make variable getter process a list of nc files
#       -> After that everything else should behave in the sa,e way for CLUBB and SAM

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
    Initialise paths, create subfolders, fill out naming templates
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
    Load netcdf files based on information from case files
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
            if 'clubb_zm' in plots.nc_files:
                nc_list.append(Dataset(cf.old_clubb_zm_file, 'r'))
            if 'clubb_zt' in plots.nc_files:
                nc_list.append(Dataset(cf.old_clubb_zt_file, 'r'))
        except IOError as e:
            logger.error('The file {}, specified in the {} case file, could not be opened: {}'.format(e.filename, cf.case, e.message))
            sys.exit()
    return nc_list

def get_t_dim(nc, create=False, dt=None):
    """
    TODO:   - implement nsave3Dstart/end
            - match time scaling to consistent units of minutes/seconds
    Read from netcdf or generate array for time dimension (In some netcdf files the dimension arrays are invalid)
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
    Read from netcdf or generate array for height (z) dimension (In some netcdf files the dimension arrays are invalid)
    TODO:   - Use grd if reading from nc fails
            - If both fails, abort
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
    Read parameters from prm file. Path is stored in case file
    """
    logger.info('get_values_from_prm')
    try:
        with open(cf.sam_prm,'r') as f:
            prm = f.read()
    except IOError as ioe:
        logger.error('The file {}, specified as sam_prm in the {} case file, could not be opened: {}'.format(ioe.filename, cf.case, ioe.message))
        return None
    prm_vars = {}
    logger.debug(prm)
    for var in prm_patterns:
        #logger.debug(var)
        try:
            logger.debug(re.search(prm_patterns[var], prm).group(1))
            prm_vars[var] = float(re.search(prm_patterns[var], prm).group(1))
        except AttributeError as ve:
            logger.error('In file {}, no entry for {} could be found.',cf.sam_prm, var)
            sys.exit()
    return prm_vars


def get_all_variables(nc, lines, plotLabels, nh, nt, t0, t1, h0, h1, filler=0):
    """
    TODO:   - Include pb.get_units in output structure?
            - Adapt variable to clubb's multiple nc files: Pass a list of nc files and try each one for the given variable name
                This would hide the differences between SAM and CLUBB
            - Use constants to access line elements
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
       2   | Profile identifier from netcdf
       3   | plot axis (0: left, 1: right)
       4   | data of variable (numpy array)
    
    (*) List structure in plots.lines:
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
      nc            -- Dataset created from netcdf file
      lines         -- List of all individual plots, which in turn contain all lines, see (*) for structure
      plotLabels    -- Identifier string for plot, used as keys for dictionary returned by this function
      nh            -- Size of height dimension in netcdf data
      nt            -- Size of time dimension in netcdf data
      t0, t1        -- Boundaries of sampling time interval
      h0, h1        -- Height bounds for output plots
    Output:
      dictionary with keys from plotLabels and entries of plot_lines as values
    """
    logger.info('get_all_variables')
    logger.debug(plotLabels)
    # Initialize list of all plots
    plot_data = []
    # Rename to make some of the expressions in SAM standalone work
    n = nh
    for i,l in enumerate(lines):
        logger.info("Getting variables for plot %s", plotLabels[i])
        # If CLD (cloud fraction) data is used in this plot,
        # we will have to evaluate expressions before averaging over sampling time interval,
        # because cloud fraction changes over time, so conditional averaging changes with each time step
        # Check if CLD needed:
        condPlot = any([x[2]=='CLD' for x in l])
        logger.debug('Cloud conditional plot? %s',condPlot)
        # Initialize list of all lines in plot l
        plot_lines = []
        # Initialize list of function lines in plot l
        functions = []
        logger.debug('Iterate over lines')
        for j,var in enumerate(l):
            logger.info("Getting variable %s", var[0])
            if var[2] is None:
                # TODO: conversion factor is not applied to function entries!!!
                logger.debug("Adding dummy line")
                plot_lines.append([var[0], var[1], var[2], var[4], None])
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
                    data = np.where(np.logical_or(np.isnan(data), data<=-1000), filler, data)
                # Average over given time indices
                logger.debug("%d>%d?",t1,t0)
                logger.debug(t1>t0)
                # TODO: Fix error when key not found -> data already has reduced size
                # t1>t0 should only apply for profile data!
                if t1>t0:
                    if not condPlot:
                        if data.shape[1]>nh:
                            data = pb.mean_profiles(data, t0, t1, h0, h1)
                        else:
                            data = pb.mean_profiles(data, t0, t1, 0, nh)
                    else:
                        # Slice here in order to have right dimensions for function calculations
                        if data.shape[1]>nh:
                            data = data[t0:t1,h0:h1]
                        logger.debug(data.shape)
                else:
                    # No time averaging
                    if data.ndim==2 and data.shape[1]>nh:
                        # 2d data here has the dimensions 0:time, 1:height(z)
                        # usually, 2d data should be averaged over time, this might not work! (TODO: test)
                        data=data[:,h0:h1]
                    elif data.ndim==3 and data.shape[0]>nh:
                        # 3d data here has the dimenions 0:z, 1:x, 2:y
                        data = data[h0:h1]
                # Save to plots, use var[2] as identifier string for replacement
                #plot_lines.append([var[0], var[1], var[2], var[4], data])
                plot_lines.append([var[0], var[1], var[2], var[4], data])
            else:
                logger.debug("Is function")
                plot_lines.append(None)
                functions.append([var[0], var[1], var[2], var[4], j])
        logger.debug('Iterating over functions')
        for func in functions:
            expression = func[2]
            logger.debug('Calculate %s', func[0])
            #logger.debug('Expression before replacement: %s',expression)
            # Insert plot_lines entries into function string
            all_nan = True # Check if all values in expression are NAN
            for k in range(len(plot_lines)):
                #logger.debug('Replace plot_line[%d]',k)
                if plot_lines[k] is not None and plot_lines[k][4] is not None:
                    if not np.all(np.isnan(plot_lines[k][4])):
                        all_nan = False
                    # Replace variables in expression with their respective entries from plot_lines and replace NAN values with zero
                    expression = expression.replace(plot_lines[k][2], 'np.where(np.isnan(plot_lines['+str(k)+'][4]),0,plot_lines['+str(k)+'][4])')
            # if all entries are NAN, return array containing only NANs
            if all_nan:
                data = np.full(nh, np.nan)
            else:
                #logger.debug('Expression after replacement: %s',expression)
                try:
                    data = eval(expression)
                    if np.any(np.isnan(data)) or np.any(np.isinf(data)):
                        # if there are invalid data points in the variable replace by filler value
                        logger.warning("Invalid data in variable %s of plot %s. The entries will be replaced by the fill value %s", var[0], plotLabels[i], str(filler))
                        data = np.where(np.isinf(data), filler, data)
                except Exception as e:
                    logger.error('Expression %s of line %s in plot %s could not be evaluated: %s', expression, plotLabels[i], func[0], e.message)
                    data = np.full(nh, np.nan)
            #logger.debug("Line insertion index: %d", func[4])
            plot_lines[func[4]] = [func[0], func[1], func[2], func[3], data]
        # Apply late averaging. Only for profile data!:
        if condPlot and t1>t0:
            logger.debug('Applying late averaging')
            for i in range(len(plot_lines)):
                data = plot_lines[i][4]
                if data is not None:
                    plot_lines[i][4] = pb.mean_profiles(data, 0, nt, 0, nh)
                    #if data.shape[1]>nh:
                        #plot_lines[i][4] = pb.mean_profiles(data, t0, t1, h0, h1)
                    #else:
                        #plot_lines[i][4] = pb.mean_profiles(data, t0, t1, 0, nh)
        plot_data.append(plot_lines)
    
    return dict(zip(plotLabels, plot_data))


def plot_default(plots, cf, data, h, centering):
    """
    Default plotting routine. Plots a series of figures containing height profiles from data.
    TODO: pass startLevel etc. to plot routine (style file?)
    """
    logger.info('plot_default')
    # Initialise names
    out_dir, jpg_dir, plot_case_name, out_pdf = init_plotgen(plots, cf)
    # Create  output pdf
    pdf = PdfPages(out_pdf)
    # Loop over plots
    imageNames = []
    
    #for i,plot_label in enumerate(plots.sortPlots):
        #logger.info('Generating plot %s', plot_label)
        #title, units = plots.plotNames[i]
        ##title = plots.plots[i][PLOT_TITLE]
        ##units = plots.plots[i][PLOT_XLABEL]
        ## pb.get_units() needs nc
        #name = plot_case_name.format(plot=plot_label)
        #imageNames.append(name)
        #if 'budget' in plots.name:
            #prefix = plots.prefix+', '
            #else:
                #prefix = ''
                #if cf.endTime>cf.startTime:
                    #title2 = ', {case}, t={startTime:.0f}-{endTime:.0f} min'.format(case=cf.full_name, startTime=cf.startTime, endTime=cf.endTime)
                    #else:
                        #title2 = ', {case}, t={startTime:.0f} min'.format(case=cf.full_name, startTime=cf.startTime)
                        #pb.plot_profiles(data[plot_label], h, units, cf.yLabel, prefix+title+title2, os.path.join(jpg_dir,name), textEntry=plots.plotText, textPos=plots.textPos,startLevel=plots.startLevel, lw=cf.lw, grid=False, centering=centering, pdf=pdf)
                        ##pb.plot_profiles(data[plot_label], h, units, cf.yLabel, title+title2, os.path.join(jpg_dir,name), startLevel=plots.startLevel, lw=cf.lw, grid=False, centering=np.any([el[3]==1 for el in data[plot_label]]), pdf=pdf)
    for plot_info in plots.plots:
        logger.info('Generating plot %s', plot_info[PLOT_ID])
        title = plot_info[PLOT_TITLE]
        units = plot_info[PLOT_XLABEL]
        # pb.get_units() needs nc
        name = plot_case_name.format(plot=plot_info[PLOT_ID])
        imageNames.append(name)
        if cf.endTime>cf.startTime:
            title2 = ', {case}, t={startTime:.0f}-{endTime:.0f} min'.format(case=cf.full_name, startTime=cf.startTime, endTime=cf.endTime)
        else:
            title2 = ', {case}, t={startTime:.0f} min'.format(case=cf.full_name, startTime=cf.startTime)
        #pb.plot_profiles(data[plot_info[PLOT_ID]], h, units, cf.yLabel, plots.prefix+title+title2, os.path.join(jpg_dir,name), textEntry=plot_info[PLOT_TEXT], textPos=plot_info[PLOT_TEXTPOS], startLevel=plots.startLevel, lw=cf.lw, grid=False, centering=centering, pdf=pdf)
        pb.plot_profiles(data[plot_info[PLOT_ID]], h, units, cf.yLabel, plots.prefix+title+title2, os.path.join(jpg_dir,name), textEntry=plot_info[PLOT_TEXT], textPos=plot_info[PLOT_TEXTPOS], lw=cf.lw, grid=False, centering=centering, pdf=pdf)
    pdf.close()
    
    # write html page
    logger.info("Writing HTML page")
    index = out_dir + 'index.html'
    mode = 'Splotgen'
    logger.debug("HTML header template: "+cf.headerText)
    headerText = cf.headerText.format(type=plots.header)
    logger.debug("HTML header: "+headerText)
    ow.writeNavPage(out_dir, headerText)
    ow.writePlotsPage(out_dir, headerText, mode, imageNames)
    ow.writeIndex(index, mode)
    
    ## Add bookmarks to pdf output:
    logger.info("Adding bookmarks to pdf output")
    writer = pdfw() # create pdf output class
    with open(out_pdf, 'rb') as inpdf: # open input file
        reader = pdfr(inpdf) # open pdf input class
        writer.appendPagesFromReader(reader) # copy all pages to output
        for i in range(writer.getNumPages()):
            #writer.addBookmark(plots.sortPlots[i], i, parent=None) # add bookmark
            writer.addBookmark(plots.plots[i][PLOT_ID], i, parent=None) # add bookmark
        writer.setPageMode("/UseOutlines") # make pdf open bookmarks overview
        with open(os.path.join(out_dir,'tmp.pdf'),'wb') as tmp_pdf: # open output pdf file
            writer.write(tmp_pdf) # write pdf content to output file
    # file streams closed
    os.rename(os.path.join(out_dir,'tmp.pdf'), out_pdf) # rename newly created pdf
    
    # Write plot parameters to params file
    param_file = open(os.path.join(out_dir, 'plot_params_{pn}_{date}'.format(pn=plots.name, date=plot_case_name.split('_')[-2])), "w")
    # dilation length dil_len
    # percentile value quant
    # cloud_interpolation cloud_interpolation
    # grid skip length skip
    # cloud recognition limit cld_lim
    # (boundary times)
    # (boundary heights)
    # (grid size)
    # (grid spacing)
    prms = [
        '# simulation parameters',
        '(nx, ny, nz) = ({}, {}, {})'.format(cf.nx, cf.ny, cf.nz),
        'sampling time = {} - {} min'.format(cf.startTime, cf.endTime),
        'dx = {} m'.format(cf.dxy),
        'dy = {} m'.format(cf.dxy),
        'dz = {} m'.format(cf.dz),
        'height = {} - {} m'.format(cf.startHeight, cf.endHeight)
        ]
    if 'clubb' in plots.name.lower():
        prms.extend([
            'CLUBB zm input file name = {}'.format(cf.clubb_zm_file),
            'CLUBB zt input file name = {}'.format(cf.clubb_zt_file)
            ])
    else:
        prms.append('SAM input file name = {}'.format(cf.sam_file))
    param_file.write('\n'.join(prms))
    param_file.close()


def plot_3d(plots, cf, data, h, h_limits, h_extent, prm_vars, fps=2, dil_len=1, gif=False):
    """
    Generates horizontal output showing cloud outlines and wind speeds
    Additionally, creates cloud conditional profiles
    TODO:   Special handling of RICO, because of huge grid -> split up grid in multiple subgrids?
            Split up creation of cloud conditionals and horizontal plots into different routines
    """
    logger.info('plot_3d')
    # Initialise names
    out_dir, jpg_dir, plot_case_name, out_pdf = init_plotgen(plots, cf)
    
    ### Get data
    logger.info('Unpacking data and perform calculations')
    ## Unpack arrays
    qn_3d = data['qn_3d']
    qv_3d = data['qv_3d']
    qt_3d = qn_3d+qv_3d
    #thv_3d = data['thetav']
    u_3d = data['u_3d']
    v_3d = data['v_3d']
    w_3d = data['w_3d']
    u_2d = data['u']
    v_2d = data['v']
    w_2d = data['w']
    ucld_2d = data['ucld']
    vcld_2d = data['vcld']
    wcld_2d = data['wcld']
    uw_2d = data['uw']
    vw_2d = data['vw']
    #logger.debug('u_2d = %s',u_2d)
    #logger.debug('v_2d = %s',v_2d)
    #logger.debug('w_2d = %s',w_2d)
    #logger.debug('ucld_2d = %s',ucld_2d)
    #logger.debug('vcld_2d = %s',vcld_2d)
    #logger.debug('wcld_2d = %s',wcld_2d)
    #logger.debug('uw_2d = %s',uw_2d)
    #logger.debug('vw_2d = %s',vw_2d)
    # interolate on the stupid s-grid
    u_3d = np.cumsum(u_3d,axis=1)
    v_3d = np.cumsum(v_3d,axis=2)
    w_3d = np.cumsum(w_3d,axis=0)
    u_3d[:,2:] = u_3d[:,2:] - u_3d[:,:-2]
    v_3d[:,:,2:] = v_3d[:,:,2:] - v_3d[:,:,:-2]
    w_3d[2:] = w_3d[2:] - w_3d[:-2]
    u_3d[:,:-1] = u_3d[:,1:]/2.
    v_3d[:,:,:-1] = v_3d[:,:,1:]/2.
    w_3d[:-1] = w_3d[1:]/2.
    # Set last entries to 0 or NAN, as they cannot b running averaged
    u_3d[:,-1] = np.nan
    v_3d[:,:,-1] = np.nan
    w_3d[-1] = np.nan
    # Due to the interpolation of w, we might have to exclude the last level of h here (TODO)
    #logger.debug('u=%s', str(u_3d.shape))
    #logger.debug('v=%s', str(v_3d.shape))
    #logger.debug('w=%s', str(w_3d.shape))
    #logger.debug('qn=%s', str(qn_3d.shape))
    ## Calculate needed values
    # Create cloud mask
    cld_mask = qn_3d>cld_lim
    cld_cnt = cld_mask.sum(axis=(1,2))
    # Create dilated cloud mask with cross kernel
    # Create dilation kernel
    # Restrict input
    #logger.debug("dil_len=%i",dil_len)
    dil_len = min(5,max(1,dil_len))
    #logger.debug("dil_len=%i",dil_len)
    # Calculate kernel size
    s = 1 + 2*dil_len
    # dot is outward most kernel layer
    dot = np.zeros((s,s))
    dot[dil_len,dil_len]=1
    # Dilate dot with cross to generate additional layers
    cross = np.zeros((3,3))
    cross[1,:]=1
    cross = np.logical_or(cross, cross.T)
    # Generate kernel
    kern = np.empty((s,s,s))
    kern[0] = dot
    kern[-1] = dot
    for i in range(1,dil_len+1):
        tmp = bin_dilation(kern[i-1],cross)
        kern[i] = tmp
        kern[-i-1] = tmp
    logger.debug("Dilation kernel:\n%s", str(kern))
    # Apply dilation
    dilated_mask = bin_dilation(cld_mask, kern)
    dilated_cnt = dilated_mask.sum(axis=(1,2))
    # Difference is halo
    halo_mask = np.logical_xor(dilated_mask, cld_mask)
    halo_cnt = halo_mask.sum(axis=(1,2))
    # Complement of dilated is no cloud mask
    nocld_mask = np.logical_not(dilated_mask)
    nocld_cnt = nocld_mask.sum(axis=(1,2))
    logger.debug("Masks cover all grid points: %s", str(np.all((dilated_mask+nocld_mask)==1)))
    # Profile of mean wind component
    um = np.nanmean(u_3d, axis=(1,2))
    vm = np.nanmean(v_3d, axis=(1,2))
    wm = np.nanmean(w_3d, axis=(1,2))
    qm = np.nanmean(qt_3d, axis=(1,2))
    #thvm = np.nanmean(thv_3d, axis=(1,2))
    # 3d field of u'w' and v'w'
    up = u_3d - um[:,None,None]
    wp = w_3d - wm[:,None,None]
    vp = v_3d - vm[:,None,None]
    qp = qt_3d - qm[:,None,None]
    #thvp = thv_3d - thvm[:,None,None]
    uw = up*wp
    vw = vp*wp
    wq = wp*qp
    #upthvp = up*thvp
    #vpthvp = vp*thvp
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
    uw_cld = np.full(wm.size, np.nan)
    vw_cld = np.full(wm.size, np.nan)
    wq_cld = np.full(wm.size, np.nan)
    uw_halo = np.full(wm.size, np.nan)
    vw_halo = np.full(wm.size, np.nan)
    wq_halo = np.full(wm.size, np.nan)
    uw_nocld = np.full(wm.size, np.nan)
    vw_nocld = np.full(wm.size, np.nan)
    wq_nocld = np.full(wm.size, np.nan)
    u_cld_mean = np.full(um.size, np.nan)
    v_cld_mean = np.full(vm.size, np.nan)
    w_cld_mean = np.full(vm.size, np.nan)
    u_halo_mean = np.full(um.size, np.nan)
    v_halo_mean = np.full(vm.size, np.nan)
    w_halo_mean = np.full(vm.size, np.nan)
    u_nocld_mean = np.full(um.size, np.nan)
    v_nocld_mean = np.full(vm.size, np.nan)
    w_nocld_mean = np.full(vm.size, np.nan)
    #upthvp_cld = np.full(um.size, np.nan)
    #upthvp_halo = np.full(um.size, np.nan)
    #upthvp_nocld = np.full(um.size, np.nan)
    #vpthvp_cld = np.full(vm.size, np.nan)
    #vpthvp_halo = np.full(vm.size, np.nan)
    #vpthvp_nocld = np.full(vm.size, np.nan)
    #thvp_cld = np.full(thvm.size, np.nan)
    #thvp_halo = np.full(thvm.size, np.nan)
    #thvp_nocld = np.full(thvm.size, np.nan)
    
    
    ### Get colormaps
    # TODO: what if uw.min*uw.max>0?
    logger.info('Setup color maps')
    ## Get cloud fraction color map
    # Calculate cloud fraction for shading of height indicator
    cloud_frac_cmap = mpl.cm.get_cmap(plots.cloud_frac_cmap)
    cloud_frac = cld_mask.mean(axis=(1,2))
    cloud_frac = cloud_frac/cloud_frac.max()
    # Generate colors from color map
    cloud_colors = cloud_frac_cmap(cloud_frac)
    ## Generate quiver colormap
    ### Get colormaps and normalizations for projection onto colormap
    ## UW
    # Get limits for uw data in order to generate colormap
    #uwlim = uw.min(), uw.max() # Do not use this, because of outliers
    # Amount of outliers seeems to be constant with data array size, but to achieve good distribution over color scale, percentiles are better
    uwlim = np.nanpercentile(uw,quant,interpolation='lower'), np.nanpercentile(uw,100-quant,interpolation='higher')
    # Get norm boundaries for color scale by discarding outliers only
    #uwsort = np.sort(uw.flatten())
    #uwlim = uwsort[trash], uwsort[-trash]
    logger.debug('uw.min=%f, uw.max=%f',np.nanmin(uw),np.nanmax(uw))
    logger.debug("uwlim=(%f,%f)",uwlim[0], uwlim[1])
    logger.debug('uw[0].min=%f,uw[0].max=%f',np.nanmin(uw[0]), np.nanmax(uw[0]))
    norm_uw = mpl.colors.Normalize(vmin=uwlim[0], vmax=uwlim[1])
    # Calculate position of 0 in normalized data (projected onto [0;1])
    zero_uw = (0-uwlim[0])/(uwlim[1]-uwlim[0])
    # Define color map segments
    cdict_uw = {"red":((0,red[0],red[0]),(zero_uw,zero_col,zero_col),(1,green[0],green[0])),
                "green":((0,red[1],red[1]),(zero_uw,zero_col,zero_col),(1,green[1],green[1])),
                "blue":((0,red[2],red[2]),(zero_uw,zero_col,zero_col),(1,green[2],green[2])),
                #"alpha":((0,1,1),(zero,0,0),(1,1,1))
                }
    # Create color map
    quiver_cmap_uw = mpl.colors.LinearSegmentedColormap('quiver_cmap_uw', cdict_uw, N=1000)
    
    
    ## VW
    # Get limits for vw data in order to generate colormap
    #vwlim = vw.min(), vw.max() # Do not use this, because of outliers
    # Amount of outliers seeems to be constant with data array size, but to achieve good distribution over color scale, percentiles are better
    vwlim = np.nanpercentile(vw,quant,interpolation='lower'), np.nanpercentile(vw,100-quant,interpolation='higher')
    # Get norm boundaries for color scale by discarding outliers only
    #vwsort = np.sort(vw.flatten())
    #vwlim = vwsort[trash], vwsort[-trash]
    logger.debug('uw.min=%f, uw.max=%f',np.nanmin(uw),np.nanmax(uw))
    logger.debug("vwlim=(%f,%f)",vwlim[0], vwlim[1])
    logger.debug('vw[0].min=%f,vw[0].max=%f',np.nanmin(vw[0]), np.nanmax(vw[0]))
    norm_vw = mpl.colors.Normalize(vmin=vwlim[0], vmax=vwlim[1])
    # Calculate position of 0 in normalized data (projected onto [0;1])
    zero_vw = (0-vwlim[0])/(vwlim[1]-vwlim[0])
    # Define color map segments
    cdict_vw = {"red":((0,red[0],red[0]),(zero_vw,zero_col,zero_col),(1,green[0],green[0])),
                "green":((0,red[1],red[1]),(zero_vw,zero_col,zero_col),(1,green[1],green[1])),
                "blue":((0,red[2],red[2]),(zero_vw,zero_col,zero_col),(1,green[2],green[2])),
                #"alpha":((0,1,1),(zero,0,0),(1,1,1))
                }
    # Create color map
    quiver_cmap_vw = mpl.colors.LinearSegmentedColormap('quiver_cmap_vw', cdict_vw, N=1000)
    
    ## UP
    # Get limits of w data in order to generate colormap, min/max still work for w, as there are no outliers
    uplim = np.nanmin(up), np.nanmax(up)
    # Use in order to exclude outliers
    #wlim = np.percentile(w,quant), np.percentile(w,100-quant)
    norm_up = mpl.colors.Normalize(vmin=uplim[0], vmax=uplim[1])
    # Calculate position of 0 in normalized data (projected onto [0;1])
    zero_up = (0-uplim[0])/(uplim[1]-uplim[0])
    # Define color map segments
    cdict_w = {"red":((0,red[0],red[0]),(zero_up,zero_col,zero_col),(1,green[0],green[0])),
               "green":((0,red[1],red[1]),(zero_up,zero_col,zero_col),(1,green[1],green[1])),
               "blue":((0,red[2],red[2]),(zero_up,zero_col,zero_col),(1,green[2],green[2])),
               #"alpha":((0,1,1),(zero,0,0),(1,1,1))
    }
    # Create color map
    quiver_cmap_up = mpl.colors.LinearSegmentedColormap('quiver_cmap_up',cdict_w, N=1000)
    
    ## VP
    # Get limits of w data in order to generate colormap, min/max still work for w, as there are no outliers
    vplim = np.nanmin(vp), np.nanmax(vp)
    # Use in order to exclude outliers
    #wlim = np.percentile(w,quant), np.percentile(w,100-quant)
    norm_vp = mpl.colors.Normalize(vmin=vplim[0], vmax=vplim[1])
    # Calculate position of 0 in normalized data (projected onto [0;1])
    zero_vp = (0-vplim[0])/(vplim[1]-vplim[0])
    # Define color map segments
    cdict_w = {"red":((0,red[0],red[0]),(zero_vp,zero_col,zero_col),(1,green[0],green[0])),
               "green":((0,red[1],red[1]),(zero_vp,zero_col,zero_col),(1,green[1],green[1])),
               "blue":((0,red[2],red[2]),(zero_vp,zero_col,zero_col),(1,green[2],green[2])),
               #"alpha":((0,1,1),(zero,0,0),(1,1,1))
               }
    # Create color map
    quiver_cmap_vp = mpl.colors.LinearSegmentedColormap('quiver_cmap_vp',cdict_w, N=1000)
    
    ## W
    # Get limits of w data in order to generate colormap, min/max still work for w, as there are no outliers
    wlim = np.nanmin(w_3d), np.nanmax(w_3d)
    # Use in order to exclude outliers
    #wlim = np.percentile(w,quant), np.percentile(w,100-quant)
    norm_w = mpl.colors.Normalize(vmin=wlim[0], vmax=wlim[1])
    # Calculate position of 0 in normalized data (projected onto [0;1])
    zero_w = (0-wlim[0])/(wlim[1]-wlim[0])
    # Define color map segments
    cdict_w = {"red":((0,red[0],red[0]),(zero_w,zero_col,zero_col),(1,green[0],green[0])),
             "green":((0,red[1],red[1]),(zero_w,zero_col,zero_col),(1,green[1],green[1])),
             "blue":((0,red[2],red[2]),(zero_w,zero_col,zero_col),(1,green[2],green[2])),
             #"alpha":((0,1,1),(zero,0,0),(1,1,1))
             }
    # Create color map
    quiver_cmap_w = mpl.colors.LinearSegmentedColormap('quiver_cmap_w',cdict_w, N=1000)
    
    ## Plotting preparations
    # Get shape of arrays, here (z,x,y)
    dims = qn_3d.shape
    N = dims[1]*dims[2]
    # Create meshgrid for 2d plot
    logger.info('Create meshgrids')
    # contour prep
    # Set hatching options:
    mpp.rcParams['hatch.color'] = hatch_col
    mpp.rcParams['hatch.linewidth'] = hw
    contour_1d = np.arange(dims[1])
    contour_xgrid, contour_ygrid= np.meshgrid(contour_1d, contour_1d)
    #logger.debug('contour_1d=%s', str(contour_1d.shape))
    #logger.debug('contour_x=%s, contour_y=%s', str(contour_xgrid.shape), str(contour_ygrid.shape))
    # quiver
    quiver_1d = np.arange(0,dims[1],skip)
    quiver_xgrid, quiver_ygrid= np.meshgrid(quiver_1d, quiver_1d)
    #logger.debug('quiver_1d=%s', str(quiver_1d.shape))
    #logger.debug('quiver_x=%s, quiver_y=%s', str(quiver_xgrid.shape), str(quiver_ygrid.shape))
    #logger.debug('um=%s, u_skip=%s', str(um.shape), str(u_3d[0,::skip,::skip].shape))

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
        #dx = prm_vars['dx'],
        dx = cf.dxy,
        dz = h_extent.mean(),
        t = (cf.time_3d*cf.dt)/60.
        )
    # Define update function for mp4/gif creation
    def update_plot(framenumber, wt):
        logger.info('Generating frame %d', framenumber)
        #logger.debug('wt: %s', wt)
        
        # Get data at height level framenumber
        #qn = qn_3d[framenumber]
        u = u_3d[framenumber]
        v = v_3d[framenumber]
        
        # Fill conditional u'w' arrays (TEST: divide by dimx*dimy or reduced because of NANs?)
        if cld_mask[framenumber].sum()>0:
            uw_cld[framenumber] = np.nanmean(uw[framenumber, cld_mask[framenumber]])
            vw_cld[framenumber] = np.nanmean(vw[framenumber, cld_mask[framenumber]])
            wq_cld[framenumber] = np.nanmean(wq[framenumber, cld_mask[framenumber]])
            u_cld_mean[framenumber] = np.nanmean(u[cld_mask[framenumber]])
            v_cld_mean[framenumber] = np.nanmean(v[cld_mask[framenumber]])
            w_cld_mean[framenumber] = np.nanmean(w_3d[framenumber, cld_mask[framenumber]])
            #uw_cld[framenumber] = np.nansum(uw[framenumber, cld_mask[framenumber]])/N
            #vw_cld[framenumber] = np.nansum(vw[framenumber, cld_mask[framenumber]])/N
            #u_cld_mean[framenumber] = np.nansum(u[cld_mask[framenumber]])/N
            #v_cld_mean[framenumber] = np.nansum(v[cld_mask[framenumber]])/N
            #w_cld_mean[framenumber] = np.nansum(w_3d[framenumber, cld_mask[framenumber]])/N
            #upthvp_cld[framenumber] = np.nansum(upthvp[framenumber, cld_mask[framenumber]])/N
            #vpthvp_cld[framenumber] = np.nansum(vpthvp[framenumber, cld_mask[framenumber]])/N
            #thvp_cld[framenumber] = np.nansum(thvp[framenumber, cld_mask[framenumber]])/N
        if halo_mask[framenumber].sum()>0:
            uw_halo[framenumber] = np.nanmean(uw[framenumber, halo_mask[framenumber]])
            vw_halo[framenumber] = np.nanmean(vw[framenumber, halo_mask[framenumber]])
            wq_halo[framenumber] = np.nanmean(wq[framenumber, halo_mask[framenumber]])
            u_halo_mean[framenumber] = np.nanmean(u[halo_mask[framenumber]])
            v_halo_mean[framenumber] = np.nanmean(v[halo_mask[framenumber]])
            w_halo_mean[framenumber] = np.nanmean(w_3d[framenumber, halo_mask[framenumber]])
            #uw_halo[framenumber] = np.nansum(uw[framenumber, halo_mask[framenumber]])/N
            #vw_halo[framenumber] = np.nansum(vw[framenumber, halo_mask[framenumber]])/N
            #u_halo_mean[framenumber] = np.nansum(u[halo_mask[framenumber]])/N
            #v_halo_mean[framenumber] = np.nansum(v[halo_mask[framenumber]])/N
            #w_halo_mean[framenumber] = np.nansum(w_3d[framenumber, halo_mask[framenumber]])/N
            #upthvp_halo[framenumber] = np.nansum(upthvp[framenumber, halo_mask[framenumber]])/N
            #vpthvp_halo[framenumber] = np.nansum(vpthvp[framenumber, halo_mask[framenumber]])/N
            #thvp_halo[framenumber] = np.nansum(thvp[framenumber, halo_mask[framenumber]])/N
        if nocld_mask[framenumber].sum()>0:
            uw_nocld[framenumber] = np.nanmean(uw[framenumber, nocld_mask[framenumber]])
            vw_nocld[framenumber] = np.nanmean(vw[framenumber, nocld_mask[framenumber]])
            wq_nocld[framenumber] = np.nanmean(wq[framenumber, nocld_mask[framenumber]])
            u_nocld_mean[framenumber] = np.nanmean(u[nocld_mask[framenumber]])
            v_nocld_mean[framenumber] = np.nanmean(v[nocld_mask[framenumber]])
            w_nocld_mean[framenumber] = np.nanmean(w_3d[framenumber, nocld_mask[framenumber]])
            #uw_nocld[framenumber] = np.nansum(uw[framenumber, nocld_mask[framenumber]])/N
            #vw_nocld[framenumber] = np.nansum(vw[framenumber, nocld_mask[framenumber]])/N
            #u_nocld_mean[framenumber] = np.nansum(u[nocld_mask[framenumber]])/N
            #v_nocld_mean[framenumber] = np.nansum(v[nocld_mask[framenumber]])/N
            #w_nocld_mean[framenumber] = np.nansum(w_3d[framenumber, nocld_mask[framenumber]])/N
            #upthvp_nocld[framenumber] = np.nansum(upthvp[framenumber, nocld_mask[framenumber]])/N
            #vpthvp_nocld[framenumber] = np.nansum(vpthvp[framenumber, nocld_mask[framenumber]])/N
            #thvp_nocld[framenumber] = np.nansum(thvp[framenumber, nocld_mask[framenumber]])/N
        
        # Clear figures
        fig.clear()
        # Prepare canvas
        ax1, ax2 = fig.subplots(1, 2, gridspec_kw={'width_ratios':width_ratio})
        #w = w_3d[framenumber]
        # Print vertical wind component array as background image
        if 'uw' in wt:
            pc = ax1.imshow(uw[framenumber], cmap=quiver_cmap_uw, norm=norm_uw, interpolation=cloud_interpolation, origin='lower')
        elif 'vw' in wt:
            pc = ax1.imshow(vw[framenumber], cmap=quiver_cmap_vw, norm=norm_vw, interpolation=cloud_interpolation, origin='lower')
        elif 'up' in wt:
            pc = ax1.imshow(up[framenumber], cmap=quiver_cmap_up, norm=norm_up, interpolation=cloud_interpolation, origin='lower')
        elif 'vp' in wt:
            pc = ax1.imshow(vp[framenumber], cmap=quiver_cmap_vp, norm=norm_vp, interpolation=cloud_interpolation, origin='lower')
        else:
            pc = ax1.imshow(w_3d[framenumber], cmap=quiver_cmap_w, norm=norm_w, interpolation=cloud_interpolation, origin='lower')
        # Plot cloud contour fields (Changed to line instead of patch) TODO: Add legend
        # TEST: increase linewidths? -> 2 or 3
        clds = ax1.contour(contour_xgrid, contour_ygrid, cld_mask[framenumber], levels=[.5,1], colors=[cld_col], extend='neither')
        ax1.contourf(contour_xgrid, contour_ygrid, cld_mask[framenumber], levels=[.5,1], colors=[(0,0,0,0)], extend='neither', hatches=['x'])
        # Plot contour of halo
        halos = ax1.contour(contour_xgrid, contour_ygrid, dilated_mask[framenumber], levels=[.5,1], colors=[halo_col], extend='neither')
        # Plot wind fields
        if 'total' in wt:
            #q = ax1.quiver(xgrid, ygrid, u, v, w, angles='xy', minshaft=2, minlength=0, scale=quiver_scale_factor[cf.case], cmap=quiver_cmap, norm=norm)
            q = ax1.quiver(quiver_xgrid, quiver_ygrid, u[::skip,::skip], v[::skip,::skip], angles='xy', minshaft=2, minlength=0, scale=quiver_scale_factor[cf.case], pivot='mid')
        else:
            # Calculate deviation from horizontal mean
            #u_dev = u - um[framenumber]
            #v_dev = v - v[framenumber]
            #q = ax1.quiver(xgrid, ygrid, u_dev, v_dev, w, angles='xy', minshaft=2, minlength=0, scale=quiver_scale_factor[cf.case]/5, cmap=quiver_cmap, norm=norm)
            q = ax1.quiver(quiver_xgrid, quiver_ygrid, u[::skip,::skip]-um[framenumber], v[::skip,::skip]-vm[framenumber], angles='xy', minshaft=2, minlength=0, scale=quiver_scale_factor[cf.case]/5, pivot='mid')
        
        # Add colorbar to ax1
        divider = make_axes_locatable(ax1)
        cax = divider.append_axes('right', size='5%',pad=.1)
        cb = fig.colorbar(pc, cax=cax)
        if 'uw' in wt:
            cb.ax.set_ylabel(r"u'w' $\mathrm{\left[\frac{m^2}{s^2}\right]}$", fontsize=fontsizes['labels'])
        elif 'vw' in wt:
            cb.ax.set_ylabel(r"v'w' $\mathrm{\left[\frac{m^2}{s^2}\right]}$", fontsize=fontsizes['labels'])
        elif 'up' in wt:
            cb.ax.set_ylabel(r"u' $\mathrm{\left[\frac{m}{s}\right]}$", fontsize=fontsizes['labels'])
        elif 'vp' in wt:
            cb.ax.set_ylabel(r"v' $\mathrm{\left[\frac{m}{s}\right]}$", fontsize=fontsizes['labels'])
        else:
            cb.ax.set_ylabel(r'w $\mathrm{\left[\frac{m}{s}\right]}$', fontsize=fontsizes['labels'])
        cb.ax.tick_params(labelsize=fontsizes['labels'])
        # Generate quiver key
        ax1.quiverkey(q, .8, 1.01, wind_strength[framenumber], 'Mean wind: {:.1f}'.format(wind_strength[framenumber])+r'$\frac{m}{s}$', angle=mean_dir[framenumber], labelsep=.2*np.nanmax(np.abs(um)), labelpos='E', fontproperties={'size':20})
        # Fill second subplot
        for line in plots.profiles:
            ax2.plot(data[line][1], h, label=data[line][0], color=col_3d[line], ls='-', lw=5)
        # Generate Rectangles
        xlims = ax2.get_xlim()
        for idx in range(len(h)):
            ax2.add_patch(mpl.patches.Rectangle((xlims[0], h_limits[idx]), xlims[1]-xlims[0], h_extent[idx], color=cloud_colors[idx]))
        # Format plot
        wt_title = wt.replace('_',' ')
        fig.suptitle(title.format(h=h[framenumber], wt=wt_title), fontsize=fontsizes['title'])
        #logger.debug("Set ylim to (%f, %f)", -cf.dxy, dims[1]*cf.dxy)
        ax1.set_xlim(-.6, dims[1]-.4)
        ax1.set_ylim(-.6, dims[2]-.4)
        ax2.set_ylim(0, h_limits[-1])
        ax2.yaxis.tick_right()
        ax2.yaxis.set_label_position('right')
        ax1.set_xlabel(r'Eastward grid dimension $\mathrm{{\left[{:.0f} m\right]}}$'.format(cf.dxy), fontsize=fontsizes['labels'])
        ax1.set_ylabel(r'Northward grid dimension $\mathrm{{\left[{:.0f} m\right]}}$'.format(cf.dxy), fontsize=fontsizes['labels'])
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
                #labels.append(int(el*cf.dxy))
            #except ValueError as ve:
                #labels.append('')
        #ax1.set_xticklabels(labels)
        ##logger.debug('xticks new: %s',str(list(ax1.get_xticklabels())))
        #labels = []
        ##logger.debug('yticks old: %s',str(list(ax1.get_yticklabels())))
        #for el in ax1.get_yticks():
            #try:
                #labels.append(int(el*cf.dxy))
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
        ticks.set_powerlimits(pow_lims)
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
    cutoff=-10
    uw_pdf = PdfPages(os.path.join(out_dir, out_pdf.format(wt='conditional_uw_profiles')))
    #uw_plots = [[r"Cloud averaged, $\overline{u'w'}^\mathrm{cld}$", True, 'dummy', 0, uw_cld[:cutoff]],
                ##[r"In-halo mean of $\mathrm{\overline{u'w'}}$", True, 'dummy', 0, uw_halo],
                ##[r"Out-cloud mean of $\mathrm{\overline{u'w'}}$", True, 'dummy', 0, uw_nocld],
                #[r"Cloud+Halo averaged, $\overline{u'w'}^\mathrm{(c\!+\!h)}$", True, 'dummy', 1, (np.nansum(np.stack((uw_cld*cld_cnt,uw_halo*halo_cnt)),0)/dilated_cnt)[:cutoff]],
                ##[r"Extended in-cloud mean of $\mathrm{\overline{u'w'}}$", True, 'dummy', 0, np.nansum(np.stack((uw_cld,uw_halo)),0)],
                ##[r"Added means of $\mathrm{\overline{u'w'}}$", True, 'dummy', 0, np.nansum(np.stack((uw_cld,uw_halo,uw_nocld)),0)],
                #[r"Layer averaged, $\overline{u'w'}$", True, 'dummy', 1, np.nansum(np.stack((uw_cld*cld_cnt,uw_halo*halo_cnt,uw_nocld*nocld_cnt))/N,0)[:cutoff]],
                ##[r"$\mathrm{\overline{u'w'}}$ from std output", True, 'dummy', 0, data['uw'][1]]
                #]
                
    uw_plots = [[r"Cloud avg., $\overline{u'w'}^\mathrm{cld}$", True, 'dummy', 0, uw_cld[:cutoff]],
                [r"Halo avg., $\overline{u'w'}^\mathrm{halo}$", True, 'dummy', 0, uw_halo[:cutoff]],
                #[r"Halo avg., $\overline{u'w'}^\mathrm{halo}$", True, 'dummy', 0, None],
                [r"Env. avg., $\overline{u'w'}^\mathrm{env}$", True, 'dummy', 0, uw_nocld[:cutoff]],
                #[r"Environment avg., $\overline{u'w'}\mathrm{env}$", True, 'dummy', 0, None],
                #[r"Cloud+Halo averaged, $\overline{u'w'}^\mathrm{(c\!+\!h)}$", True, 'dummy', 1, (np.nansum(np.stack((uw_cld*cld_cnt,uw_halo*halo_cnt)),0)/dilated_cnt)[:cutoff]],
                [r"Cloud+Halo avg., $\overline{u'w'}^\mathrm{(c\!+\!h)}$", True, 'dummy', 0, None],
                #[r"Extended in-cloud mean of $\mathrm{\overline{u'w'}}$", True, 'dummy', 0, np.nansum(np.stack((uw_cld,uw_halo)),0)],
                #[r"Added means of $\mathrm{\overline{u'w'}}$", True, 'dummy', 0, np.nansum(np.stack((uw_cld,uw_halo,uw_nocld)),0)],
                [r"Layer avg., $\overline{u'w'}$", True, 'dummy', 0, np.nansum(np.stack((uw_cld*cld_cnt,uw_halo*halo_cnt,uw_nocld*nocld_cnt))/N,0)[:cutoff]],
                #[r"$\mathrm{\overline{u'w'}}$ from std output", True, 'dummy', 0, data['uw'][1]]
                ]
    vw_plots = [[r"Cloud averaged, $\overline{v'w'}^\mathrm{cld}$", True, 'dummy', 0, vw_cld[:cutoff]],
                [r"Halo averaged, $\overline{v'w'}^\mathrm{halo}$", True, 'dummy', 0, vw_halo[:cutoff]],
                [r"Environment averaged, $\overline{v'w'}^\mathrm{env}$", True, 'dummy', 0, vw_nocld[:cutoff]],
                [r"Cloud+Halo averaged, $\overline{v'w'}^\mathrm{(c\!+\!h)}$", True, 'dummy', 0, (np.nansum(np.stack((vw_cld*cld_cnt,vw_halo*halo_cnt)),0)/dilated_cnt)[:cutoff]],
                #[r"Extended in-cloud mean of $\mathrm{\overline{v'w'}}$", True, 'dummy', 0, np.nansum(np.stack((vw_cld,vw_halo)),0)],
                #[r"Added means of $\mathrm{\overline{v'w'}}$", True, 'dummy', 0, np.nansum(np.stack((vw_cld,vw_halo,vw_nocld)),0)],
                [r"Layer averaged, $\overline{v'w'}$", True, 'dummy', 1, np.nansum(np.stack((vw_cld*cld_cnt,vw_halo*halo_cnt,vw_nocld*nocld_cnt))/N,0)[:cutoff]],
                #[r"$\mathrm{\overline{v'w'}}$ from std output", True, 'dummy', 0, data['vw'][1]]
                ]
    
    wq_plots = [[r"Cloud weighted, $\overline{w'q_t'}^\mathrm{cld}$", True, 'dummy', 0, wq_cld[:cutoff]*cld_cnt[:cutoff]/N],
                #[r"Halo averaged, $\overline{u'w'}^\mathrm{halo}$", True, 'dummy', 0, uw_halo],
                [r"Halo weighted, $\overline{w'q_t'}^\mathrm{halo}$", True, 'dummy', 0, None],
                #[r"Environment averaged, $\overline{u'w'}\mathrm{env}$", True, 'dummy', 0, uw_nocld],
                [r"Environment weighted, $\overline{w'q_t'}^\mathrm{env}$", True, 'dummy', 0, (wq_halo*halo_cnt+wq_nocld*nocld_cnt)[:cutoff]/N],
                #[r"Cloud+Halo averaged, $\overline{u'w'}^\mathrm{(c\!+\!h)}$", True, 'dummy', 1, (np.nansum(np.stack((uw_cld*cld_cnt,uw_halo*halo_cnt)),0)/dilated_cnt)[:cutoff]],
                [r"Cloud+Halo avg., $\overline{w'q_t'}^\mathrm{(c\!+\!h)}$", True, 'dummy', 0, None],
                #[r"Extended in-cloud mean of $\mathrm{\overline{u'w'}}$", True, 'dummy', 0, np.nansum(np.stack((uw_cld,uw_halo)),0)],
                #[r"Added means of $\mathrm{\overline{u'w'}}$", True, 'dummy', 0, np.nansum(np.stack((uw_cld,uw_halo,uw_nocld)),0)],
                [r"Layer avg., $\overline{w'q_t'}$", True, 'dummy', 0, np.nansum(np.stack((wq_cld*cld_cnt,wq_halo*halo_cnt,wq_nocld*nocld_cnt)),0)[:cutoff]/N],
                #[r"$\mathrm{\overline{u'w'}}$ from std output", True, 'dummy', 0, data['uw'][1]]
                ]
    
    um_plots = [[r"Cloud avg., $\bar{u}^\mathrm{cld}$", True, 'dummy', 0, u_cld_mean[:cutoff]],
                [r"Halo avg., $\bar{u}^\mathrm{halo}$", True, 'dummy', 0, u_halo_mean[:cutoff]],
                [r"Env. avg.,  $\bar{u}^\mathrm{env}$", True, 'dummy', 0, u_nocld_mean[:cutoff]],
                #[r"Cloud+Halo averaged, $\overline{u}^\mathrm{(c\!+\!h)}$", True, 'dummy', 0, (np.nansum(np.stack((u_cld_mean*cld_cnt,u_halo_mean*halo_cnt)),0)/dilated_cnt)[:cutoff]],
                [r"Cloud+Halo avg., $\bar{u}^\mathrm{(c\!+\!h)}$", True, 'dummy', 0, None],
                #["Extended in-cloud mean of u", True, 'dummy', 0, np.nansum(np.stack((u_cld_mean,u_halo_mean)),0)],
                #["Extended in-cloud mean of u", True, 'dummy', 0, None],
                #["Added means of u", True, 'dummy', 0, np.nansum(np.stack((u_cld_mean,u_halo_mean,u_nocld_mean)),0)],
                [r"Layer avg., $\bar{u}$", True, 'dummy', 0, np.nansum(np.stack((u_cld_mean*cld_cnt,u_halo_mean*halo_cnt,u_nocld_mean*nocld_cnt))/N,0)[:cutoff]],
                # Add means from std ouput
                #[r"$\mathrm{\bar{u}}$ from std output", True, 'dummy', 0, data['u'][1]],
                #[r"In-cloud mean of u from std output", True, 'dummy', 0, data['ucld'][1]]
                ]
    vm_plots = [["Cloud averaged, v", True, 'dummy', 0, v_cld_mean[:cutoff]],
                ["Halo averaged, v", True, 'dummy', 0, v_halo_mean[:cutoff]],
                ["Environment averaged, v", True, 'dummy', 0, v_nocld_mean[:cutoff]],
                ["Cloud+Halo averaged, v", True, 'dummy', 0, (np.nansum(np.stack((v_cld_mean*cld_cnt,v_halo_mean*halo_cnt)),0)/dilated_cnt)[:cutoff]],
                #["Extended in-cloud mean of v", True, 'dummy', 0, np.nansum(np.stack((v_cld_mean,v_halo_mean)),0)],
                #["Added means of v", True, 'dummy', 0, np.nansum(np.stack((v_cld_mean,v_halo_mean,v_nocld_mean)),0)],
                ["Layer averaged, v", True, 'dummy', 0, np.nansum(np.stack((v_cld_mean*cld_cnt,v_halo_mean*halo_cnt,v_nocld_mean*nocld_cnt))/N,0)[:cutoff]],
                # Add means from std ouput
                #[r"$\mathrm{\bar{v}}$ from std output", True, 'dummy', 0, data['v'][1]],
                #[r"In-cloud mean of v from std output", True, 'dummy', 0, data['vcld'][1]]
                ]
    wm_plots = [["Cloud averaged, w", True, 'dummy', 0, w_cld_mean[:cutoff]],
                ["Halo averaged, w", True, 'dummy', 0, w_halo_mean[:cutoff]],
                ["Environment averaged, w", True, 'dummy', 0, w_nocld_mean[:cutoff]],
                ["Extended in-cloud mean of w", True, 'dummy', 0, (np.nansum(np.stack((w_cld_mean*cld_cnt,w_halo_mean*halo_cnt)),0)/dilated_cnt)[:cutoff]],
                #["Extended in-cloud mean of w", True, 'dummy', 0, np.nansum(np.stack((w_cld_mean,w_halo_mean)),0)],
                #["Added means of w", True, 'dummy', 0, np.nansum(np.stack((w_cld_mean,w_halo_mean,w_nocld_mean)),0)],
                ["Added means of w", True, 'dummy', 0, np.nansum(np.stack((w_cld_mean*cld_cnt,w_halo_mean*halo_cnt,w_nocld_mean*nocld_cnt))/N,0)[:cutoff]],
                # Add means from std ouput
                #[r"$\mathrm{\bar{w}}$ from std output", True, 'dummy', 0, data['w'][1]],
                #[r"In-cloud mean of w from std output", True, 'dummy', 0, data['wcld'][1]]
                ]
    # theta_v not included in 3d data
    #upthvp_plots = [[r"In-cloud mean of $\mathrm{\overline{u'\theta_v'}}$", True, 'dummy', 0, upthvp_cld],
                    #[r"In-halo mean of $\mathrm{\overline{u'\theta_v'}}$", True, 'dummy', 0, upthvp_halo],
                    #[r"Extended in-cloud mean of $\mathrm{\overline{u'\theta_v'}}$", True, 'dummy', 0, np.nansum(np.stack((upthvp_cld,upthvp_halo)), 0)],
                    #[r"Out-cloud mean of $\mathrm{\overline{u'\theta_v'}}$", True, 'dummy', 0, upthvp_nocld],
                    #[r"Total mean of $\mathrm{\overline{u'\theta_v'}}$", True, 'dummy', 0, upthvp_total],
                    #[r"Added means of $\mathrm{\overline{u'\theta_v'}}$", True, 'dummy', 0, np.nansum(np.stack((upthvp_cld,upthvp_halo,upthvp_nocld)),0)]
                    #]
    #vpthvp_plots = [["In-cloud mean of $\mathrm{\overline{v'\theta_v'}}$", True, 'dummy', 0, vpthvp_cld],
                    #["In-halo mean of $\mathrm{\overline{v'\theta_v'}}$", True, 'dummy', 0, vpthvp_halo],
                    #["Extended in-cloud mean of $\mathrm{\overline{v'\theta_v'}}$", True, 'dummy', 0, np.nansum(np.stack((vpthvp_cld,vpthvp_halo)),0)],
                    #["Out-cloud mean of $\mathrm{\overline{v'\theta_v'}}$", True, 'dummy', 0, vpthvp_nocld],
                    #["Total mean of $\mathrm{\overline{v'\theta_v'}}$", True, 'dummy', 0, vpthvp_total],
                    #["Added means of $\mathrm{\overline{v'\theta_v'}}$", True, 'dummy', 0, np.nansum(np.stack((vpthvp_cld,vpthvp_halo,vpthvp_nocld)),0)]
                    #]
    #thvp_plots = [["In-cloud mean of $\mathrm{\overline{\theta_v'}}$", True, 'dummy', 0, thvp_cld],
                  #["In-halo mean of $\mathrm{\overline{\theta_v'}}$", True, 'dummy', 0, thvp_halo],
                  #["Extended in-cloud mean of $\mathrm{\overline{\theta_v'}}$", True, 'dummy', 0, np.nansum(np.stack((thvp_cld,thvp_halo)),0)],
                  #["Out-cloud mean of $\mathrm{\overline{\theta_v'}}$", True, 'dummy', 0, thvp_nocld],
                  #["Total mean of $\mathrm{\overline{\theta_v'}}$", True, 'dummy', 0, thvp_total],
                  #["Added means of $\mathrm{\overline{\theta_v'}}$", True, 'dummy', 0, np.nansum(np.stack((thvp_cld,thvp_halo,thvp_nocld)),0)]
                  #]
    # UW plots
    for k in range(len(uw_plots)):
        #pb.plot_profiles([uw_plots[k]], h[:cutoff], r"Cloud Conditional $\mathrm{\overline{u'w'}\ \left[\frac{m^2}{s^2}\right]}$", cf.yLabel, uw_plots[k][0], os.path.join(jpg_dir,'uw_cld{}'.format(k)), startLevel=0, lw=cf.lw, grid=False, centering=False, pdf=uw_pdf)
        pb.plot_profiles([uw_plots[k]], h[:cutoff], r"Cloud Conditional $\mathrm{\overline{u'w'}\ \left[\frac{m^2}{s^2}\right]}$", cf.yLabel, uw_plots[k][0], os.path.join(jpg_dir,'uw_cld{}'.format(k)), lw=cf.lw, grid=False, centering=False, pdf=uw_pdf)

    #pb.plot_profiles([["u'^2",True,'dummy',0,np.nanmean(up*up, axis=(1,2))]], h, "u'^2", cf.yLabel, "u'^2", os.path.join(jpg_dir,'up2'), startLevel=0, lw=cf.lw, grid=False, centering=False, pdf=uw_pdf)
    pb.plot_profiles([["u'^2",True,'dummy',0,np.nanmean(up*up, axis=(1,2))]], h, "u'^2", cf.yLabel, "u'^2", os.path.join(jpg_dir,'up2'), lw=cf.lw, grid=False, centering=False, pdf=uw_pdf)
    
    #pb.plot_profiles(uw_plots, h[:cutoff], r"Cloud conditional $\mathrm{\overline{u'w'}\ \left[\frac{m^2}{s^2}\right]}$", cf.yLabel, r"Conditional means of $\mathrm{\overline{u'w'}}$", os.path.join(jpg_dir,'uw_cld_all'), startLevel=0, lw=cf.lw, grid=False, centering=True, pdf=uw_pdf)
    pb.plot_profiles(uw_plots, h[:cutoff], r"Cloud conditional $\mathrm{\overline{u'w'}\ \left[\frac{m^2}{s^2}\right]}$", cf.yLabel, r"Conditional means of $\mathrm{\overline{u'w'}}$", os.path.join(jpg_dir,'uw_cld_all'), lw=cf.lw, grid=False, centering=True, pdf=uw_pdf)
    
    # TODO: Add standard mean u
    ##pb.plot_profiles([["Extended in-cloud mean of u",True,0,u_cld_mean]], h, r"Extended Cloud Conditional $\mathrm{\bar{u}\ \left[\frac{m}{s}\right]}$", cf.yLabel, "Extended Conditional mean of u", os.path.join(jpg_dir,'u_cld_mean'), startLevel=0, lw=cf.lw, grid=False, centering=False, pdf=uw_pdf)
    #pb.plot_profiles([["Extended in-cloud mean of u",True,0,u_cld_mean]], h, r"Extended Cloud Conditional $\mathrm{\bar{u}\ \left[\frac{m}{s}\right]}$", cf.yLabel, "Extended Conditional mean of u", os.path.join(jpg_dir,'u_cld_mean'), lw=cf.lw, grid=False, centering=False, pdf=uw_pdf)
    
    # U plots
    #pb.plot_profiles(um_plots, h[:cutoff], r"Cloud conditional $\mathrm{\bar{u}\ \left[\frac{m}{s}\right]}$", cf.yLabel, "Comparison of conditional means of u", os.path.join(jpg_dir, 'u_cld_all'), startLevel=0, lw=cf.lw, grid=False, centering=False, pdf=uw_pdf)
    pb.plot_profiles(um_plots, h[:cutoff], r"Cloud conditional $\mathrm{\bar{u}\ \left[\frac{m}{s}\right]}$", cf.yLabel, "Comparison of conditional means of u", os.path.join(jpg_dir, 'u_cld_all'), lw=cf.lw, grid=False, centering=False, pdf=uw_pdf)

    # VW plots
    for k in range(len(vw_plots)):
        #pb.plot_profiles([vw_plots[k]], h[:cutoff], r"Cloud conditional $\mathrm{\overline{v'w'}\ \left[\frac{m^2}{s^2}\right]}$", cf.yLabel, vw_plots[k][0], os.path.join(jpg_dir,'vw_cld{}'.format(k)), startLevel=0, lw=cf.lw, grid=False, centering=False, pdf=uw_pdf)
        pb.plot_profiles([vw_plots[k]], h[:cutoff], r"Cloud conditional $\mathrm{\overline{v'w'}\ \left[\frac{m^2}{s^2}\right]}$", cf.yLabel, vw_plots[k][0], os.path.join(jpg_dir,'vw_cld{}'.format(k)), lw=cf.lw, grid=False, centering=False, pdf=uw_pdf)
    
    #pb.plot_profiles([["v'^2",True,'dummy', 0,np.nanmean(vp*vp, axis=(1,2))]], h, "v'^2", cf.yLabel, "v'^2", os.path.join(jpg_dir,'vp2'), startLevel=0, lw=cf.lw, grid=False, centering=False, pdf=uw_pdf)
    pb.plot_profiles([["v'^2",True,'dummy', 0,np.nanmean(vp*vp, axis=(1,2))]], h, "v'^2", cf.yLabel, "v'^2", os.path.join(jpg_dir,'vp2'), lw=cf.lw, grid=False, centering=False, pdf=uw_pdf)
    
    #pb.plot_profiles(vw_plots, h[:cutoff], r"Cloud conditional $\mathrm{\overline{v'w'}\ \left[\frac{m^2}{s^2}\right]}$", cf.yLabel, r"Conditional means of $\mathrm{\overline{v'w'}}$", os.path.join(jpg_dir,'vw_cld_all'), startLevel=0, lw=cf.lw, grid=False, centering=True, pdf=uw_pdf)
    pb.plot_profiles(vw_plots, h[:cutoff], r"Cloud conditional $\mathrm{\overline{v'w'}\ \left[\frac{m^2}{s^2}\right]}$", cf.yLabel, r"Conditional means of $\mathrm{\overline{v'w'}}$", os.path.join(jpg_dir,'vw_cld_all'), lw=cf.lw, grid=False, centering=True, pdf=uw_pdf)
    
    #pb.plot_profiles(wq_plots, h[:cutoff], r"Cloud conditional $\mathrm{\overline{w'q_n'}\ \left[g\,m\,kg^{-1}\,s^{-1}\right]}$", cf.yLabel, r"Conditional means of $\overline{w'q_n'}}$", os.path.join(jpg_dir,'wq_cld_all'), startLevel=0, lw=cf.lw, grid=False, centering=True, pdf=uw_pdf)
    pb.plot_profiles(wq_plots, h[:cutoff], r"Cloud conditional $\mathrm{\overline{w'q_n'}\ \left[g\,m\,kg^{-1}\,s^{-1}\right]}$", cf.yLabel, r"Conditional means of $\overline{w'q_n'}}$", os.path.join(jpg_dir,'wq_cld_all'), lw=cf.lw, grid=False, centering=True, pdf=uw_pdf)
    
    # TODO: Add standard mean v
    ##pb.plot_profiles([["Extended in-cloud mean of v",True,0,v_cld_mean]], h, r"Extended Cloud Conditional $\mathrm{\bar{v}\ \left[\frac{m}{s}\right]}$", cf.yLabel, "Extended Conditional mean of v", os.path.join(jpg_dir,'v_cld_mean'), startLevel=0, lw=cf.lw, grid=False, centering=False, pdf=uw_pdf)
    #pb.plot_profiles([["Extended in-cloud mean of v",True,0,v_cld_mean]], h, r"Extended Cloud Conditional $\mathrm{\bar{v}\ \left[\frac{m}{s}\right]}$", cf.yLabel, "Extended Conditional mean of v", os.path.join(jpg_dir,'v_cld_mean'), lw=cf.lw, grid=False, centering=False, pdf=uw_pdf)
    
    # V plots
    #pb.plot_profiles(vm_plots, h[:cutoff], r"Cloud conditional $\mathrm{\bar{v}\ \left[\frac{m}{s}\right]}$", cf.yLabel, "Comparison of conditional means of v", os.path.join(jpg_dir, 'v_cld_all'), startLevel=0, lw=cf.lw, grid=False, centering=False, pdf=uw_pdf)
    pb.plot_profiles(vm_plots, h[:cutoff], r"Cloud conditional $\mathrm{\bar{v}\ \left[\frac{m}{s}\right]}$", cf.yLabel, "Comparison of conditional means of v", os.path.join(jpg_dir, 'v_cld_all'), lw=cf.lw, grid=False, centering=False, pdf=uw_pdf)
    
    # W plots
    #pb.plot_profiles(wm_plots, h[:cutoff], r"Cloud conditional $\mathrm{\bar{w}\ \left[\frac{m}{s}\right]}$", cf.yLabel, "Comparison of conditional means of w", os.path.join(jpg_dir, 'w_cld_all'), startLevel=0, lw=cf.lw, grid=False, centering=False, pdf=uw_pdf)
    pb.plot_profiles(wm_plots, h[:cutoff], r"Cloud conditional $\mathrm{\bar{w}\ \left[\frac{m}{s}\right]}$", cf.yLabel, "Comparison of conditional means of w", os.path.join(jpg_dir, 'w_cld_all'), lw=cf.lw, grid=False, centering=False, pdf=uw_pdf)
    
    ## UTHV plot
    #for k in range(len(upthvp_plots)):
        #pb.plot_profiles([upthvp_plots[k]], h, r"Cloud Conditional $\mathrm{\overline{u'\theta_v'}\ \left[\frac{m^2K^2}{s^2}\right]}$", cf.yLabel, upthvp_plots[k][0], os.path.join(jpg_dir,'upthvp_cld{}'.format(k)), startLevel=0, lw=cf.lw, grid=False, centering=False, pdf=uw_pdf)
    
    #pb.plot_profiles(upthvp_plots, h, r"Cloud Conditional $\mathrm{\overline{u'\theta_v'}\ \left[\frac{m^2K^2}{s^2}\right]}$", cf.yLabel, r"Conditional $\mathrm{\overline{u'\theta_v'}}$", os.path.join(jpg_dir,'upthvp_cld_all'), startLevel=0, lw=cf.lw, grid=False, centering=False, pdf=uw_pdf)
    
    ## VTHV plot
    #for k in range(len(vpthvp_plots)):
        #pb.plot_profiles([vpthvp_plots[k]], h, r"Cloud Conditional $\mathrm{\overline{v'\theta_v'}\ \left[\frac{m^2K^2}{s^2}\right]}$", cf.yLabel, vpthvp_plots[k][0], os.path.join(jpg_dir,'vpthvp_cld{}'.format(k)), startLevel=0, lw=cf.lw, grid=False, centering=False, pdf=uw_pdf)
        
    #pb.plot_profiles(vpthvp_plots, h, r"Cloud Conditional $\mathrm{\overline{v'\theta_v'}\ \left[\frac{m^2K^2}{s^2}\right]}$", cf.yLabel, r"Conditional $\mathrm{\overline{v'\theta_v'}}$", os.path.join(jpg_dir,'vpthvp_cld_all'), startLevel=0, lw=cf.lw, grid=False, centering=False, pdf=uw_pdf)
        
    ## THV plot
    #for k in range(len(thvp_plots)):
        #pb.plot_profiles([thvp_plots[k]], h, r"Cloud Conditional $\mathrm{\overline{\theta_v'}\ \left[K^2}\right]}$", cf.yLabel, thvp_plots[k][0], os.path.join(jpg_dir,'thvp_cld{}'.format(k)), startLevel=0, lw=cf.lw, grid=False, centering=False, pdf=uw_pdf)

    #pb.plot_profiles(thvp_plots, h, r"Cloud Conditional $\mathrm{\overline{\theta_v'}\ \left[K^2\right]}$", cf.yLabel, r"Conditional $\mathrm{\overline{\theta_v'}}$", os.path.join(jpg_dir,'thvp_cld_all'), startLevel=0, lw=cf.lw, grid=False, centering=False, pdf=uw_pdf)

    # UM plots
    
    uw_pdf.close()
    
    # Write plot parameters to params file
    param_file = open(os.path.join(out_dir, 'plot_3d_params_{pn}_{date}'.format(pn=plots.name, date=plot_case_name.split('_')[-2])), "w")
    # dilation length dil_len
    # percentile value quant
    # cloud_interpolation cloud_interpolation
    # grid skip length skip
    # cloud recognition limit cld_lim
    # (boundary times)
    # (boundary heights)
    # (grid size)
    # (grid spacing)
    prms = [
        '# simulation parameters',
        '(nx, ny, nz) = ({}, {}, {})'.format(dims[1], dims[2], dims[0]),
        't = {}'.format(float(cf.time_3d)*float(cf.dt)),
        'dx = {}'.format(cf.dxy),
        'dy = {}'.format(cf.dxy),
        'dz = {}'.format(cf.dz),
        'height = {} - {}'.format(cf.startHeight, cf.endHeight),
        '# plotting parameters:',
        'number of arrows skipped = {}'.format(skip),
        'cloud recognition limit for water vapor = {}'.format(cld_lim),
        'background interpolation method = {}'.format(cloud_interpolation),
        'cloud halo dilation length'.format(dil_len),
        'convariance percentile cutoff = {}'.format(quant),
        '3D input file name = {}'.format(cf.sam_3d_file),
        '2D input file name = {}'.format(cf.sam_file),
        ]
    param_file.write('\n'.join(prms))
    param_file.close()
    
    

def plot_comparison(plots, cf, data_clubb, data_sam, h_clubb, h_sam, plot_old_clubb=False, data_old=None, h_old=None):
    """
    Create a series of figures containing profile plots from SAM and CLUBB data.
    Added functionality to include old CLUBB data in plots.
    TODO: get units from clubb long name
    """
    logger.info('plot_comparison')
    # Initialise names
    out_dir, jpg_dir, plot_case_name, out_pdf = init_plotgen(plots, cf)
    # Create  output pdf
    pdf = PdfPages(out_pdf)
    # Create image name list for html page
    imageNames = []
    
    # Loop over plots
    #for i, plot_label in enumerate(plots.sortPlots):
        #logger.info('Generating plot %s', plot_label)
        #title, units = plots.plotNames[i]
        ##title = plots.plots[i][PLOT_TITLE]
        ##units = plots.plots[i][PLOT_XLABEL]
        #if 'budget' in plots.name.lower():
            #prefix = plots.prefix+', '
        #else:
            #prefix = ''
        #if cf.endTime>cf.startTime:
            #title2 = ', {case}, t={startTime:.0f}-{endTime:.0f} min'.format(case=cf.full_name, startTime=cf.startTime, endTime=cf.endTime)
        #else:
            #title2 = ', {case}, t={startTime:.0f} min'.format(case=cf.full_name, startTime=cf.startTime)
        #if plot_old_clubb:
            #pb.plot_comparison(data_clubb[plot_label], data_sam[plot_label], h_clubb, h_sam, units, cf.yLabel, prefix+title+title2, os.path.join(jpg_dir, plot_case_name.format(plot=plot_label)), textEntry=plots.plotText, textPos=plots.textPos, startLevel=plots.startLevel, grid=False, pdf=pdf, plot_old_clubb=plot_old_clubb, data_old=data_old[plot_label], level_old=h_old)
        #else:
            #pb.plot_comparison(data_clubb[plot_label], data_sam[plot_label], h_clubb, h_sam, units, cf.yLabel, prefix+title+title2, os.path.join(jpg_dir, plot_case_name.format(plot=plot_label)), textEntry=plots.plotText, textPos=plots.textPos, startLevel=plots.startLevel, grid=False, pdf=pdf)
    
    for plot_info in plots.plots:
        #logger.debug(plot_info)
        logger.info('Generating plot %s', plot_info[PLOT_ID])
        #title, units = plots.plotNames[i]
        title = plot_info[PLOT_TITLE]
        units = plot_info[PLOT_XLABEL]
        name = plot_case_name.format(plot=plot_info[PLOT_ID])
        imageNames.append(name)
        if cf.endTime>cf.startTime:
            title2 = ', {case}, t={startTime:.0f}-{endTime:.0f} min'.format(case=cf.full_name, startTime=cf.startTime, endTime=cf.endTime)
        else:
            title2 = ', {case}, t={startTime:.0f} min'.format(case=cf.full_name, startTime=cf.startTime)
        if plot_old_clubb:
            #pb.plot_comparison(data_clubb[plot_info[PLOT_ID]], data_sam[plot_info[PLOT_ID]], h_clubb, h_sam, units, cf.yLabel, plots.prefix+title+title2, os.path.join(jpg_dir, plot_case_name.format(plot=plot_info[PLOT_ID])), textEntry=plot_info[PLOT_TEXT], textPos=plot_info[PLOT_TEXTPOS], startLevel=plots.startLevel, grid=False, pdf=pdf, plot_old_clubb=plot_old_clubb, data_old=data_old[plot_info[PLOT_ID]], level_old=h_old)
            pb.plot_comparison(data_clubb[plot_info[PLOT_ID]], data_sam[plot_info[PLOT_ID]], h_clubb, h_sam, units, cf.yLabel, plots.prefix+title+title2, os.path.join(jpg_dir, name), textEntry=plot_info[PLOT_TEXT], textPos=plot_info[PLOT_TEXTPOS], grid=False, pdf=pdf, plot_old_clubb=plot_old_clubb, data_old=data_old[plot_info[PLOT_ID]], level_old=h_old)
        else:
            #pb.plot_comparison(data_clubb[plot_info[PLOT_ID]], data_sam[plot_info[PLOT_ID]], h_clubb, h_sam, units, cf.yLabel, plots.prefix+title+title2, os.path.join(jpg_dir, plot_case_name.format(plot=plot_info[PLOT_ID])), textEntry=plot_info[PLOT_TEXT], textPos=plot_info[PLOT_TEXTPOS], startLevel=plots.startLevel, grid=False, pdf=pdf)
            pb.plot_comparison(data_clubb[plot_info[PLOT_ID]], data_sam[plot_info[PLOT_ID]], h_clubb, h_sam, units, cf.yLabel, plots.prefix+title+title2, os.path.join(jpg_dir, name), textEntry=plot_info[PLOT_TEXT], textPos=plot_info[PLOT_TEXTPOS], grid=False, pdf=pdf)
    pdf.close()
    
    # write html page
    logger.info("Writing HTML page")
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
            #writer.addBookmark(plots.sortPlots[i], i, parent=None) # add bookmark
            writer.addBookmark(plots.plots[i][PLOT_ID], i, parent=None) # add bookmark
        writer.setPageMode("/UseOutlines") # make pdf open bookmarks overview
        with open(os.path.join(out_dir,'tmp.pdf'),'wb') as tmp_pdf: # open output pdf file
            writer.write(tmp_pdf) # write pdf content to output file
    # file streams closed
    os.rename(os.path.join(out_dir,'tmp.pdf'), out_pdf) # rename newly created pdf
    
    # Write plot parameters to params file TODO
    param_file = open(os.path.join(out_dir, 'plot_params_{pn}_{date}'.format(pn=plots.name, date=plot_case_name.split('_')[-2])), "w")
    # dilation length dil_len
    # percentile value quant
    # cloud_interpolation cloud_interpolation
    # grid skip length skip
    # cloud recognition limit cld_lim
    # (boundary times)
    # (boundary heights)
    # (grid size)
    # (grid spacing)
    prms = [
        '# simulation parameters',
        '(nx, ny, nz) = ({}, {}, {})'.format(cf.nx, cf.ny, cf.nz),
        'sampling time = {}-{}min'.format(cf.startTime, cf.endTime),
        'dx = {}'.format(cf.dxy),
        'dy = {}'.format(cf.dxy),
        'dz = {}'.format(cf.dz),
        'height = {} - {}'.format(cf.startHeight, cf.endHeight),
        'SAM input file name = {}'.format(cf.sam_file),
        'New CLUBB zm input file name = {}'.format(cf.clubb_zm_file),
        'New CLUBB zt input file name = {}'.format(cf.clubb_zt_file),
        'Old CLUBB zm input file name = {}'.format(cf.old_clubb_zm_file),
        'Old CLUBB zt input file name = {}'.format(cf.old_clubb_zt_file),
        ]
    param_file.write('\n'.join(prms))
    param_file.close()



#-------------------------------------------------------------------------------
#   MAIN PLOT FUNCTIONS
#-------------------------------------------------------------------------------

def plotgen_default(plots, cf):
    """
    Routine setting up default height profile plotting for either SAM or CLUBB data.
    Data is processed here and afterwards, the plot_default routine is called to create plots.
    """
    logger.info("plotgen_default")
    ncs = load_nc(plots, cf)
    if 'clubb' in plots.name.lower():
        # Separate CLUBB nc files
        #if plots.sortPlots_zm and plots.sortPlots_zt:
        if all([ncfile in plots.nc_files for ncfile in ['clubb_zm', 'clubb_zt']]):
            nc_zm, nc_zt = ncs
            nc = nc_zm
        elif 'clubb_zm' in plots.nc_files:
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
    if 'clubb' in plots.name.lower():
        h = get_h_dim(nc, model='clubb')
    else:
        h = get_h_dim(nc_sam, model='sam')
    logger.info('Find nearest levels to case setup')
    idx_h0 = (np.abs(h - cf.startHeight)).argmin()
    idx_h1 = min((np.abs(h - cf.endHeight)).argmin()+1, h.size-1)
    logger.debug("Height indices: %d, %d", idx_h0, idx_h1)
    logger.debug("Height limits: %f, %f", h[idx_h0], h[idx_h1])
    # Times for CLUBB are given in seconds, for SAM in minutes
    if 'clubb' in plots.name.lower():
        idx_t0 = (np.abs(t - cf.startTime*60)).argmin()
        idx_t1 = min((np.abs(t - cf.endTime*60)).argmin()+1, t.size-1)
    else:
        idx_t0 = (np.abs(t - cf.startTime)).argmin()
        idx_t1 = min((np.abs(t - cf.endTime)).argmin()+1, t.size-1)
    h = h[idx_h0:idx_h1]
    nt = idx_t1 - idx_t0
    logger.debug('t0 = %f, t1 = %f, nt = %d', idx_t0, idx_t1, nt)
    
    # TODO: figure out how to adapt this to the new structure in the setup file
    if 'clubb' in plots.name.lower():
        # Initialize data variables as empty dicts
        data_zm = {}
        data_zt = {}
        if nc_zm:
            logger.info("Fetching clubb_zm data")
            #data_zm = get_all_variables(nc_zm, plots.lines_zm, plots.sortPlots_zm, len(h), len(t), idx_t0, idx_t1, idx_h0, idx_h1, filler=plots.filler)
            data_zm = get_all_variables(nc_zm, [l[PLOT_LINES] for l in plots.plots_zm], [l[PLOT_ID] for l in plots.plots_zm], len(h), nt, idx_t0, idx_t1, idx_h0, idx_h1, filler=plots.filler)
        if nc_zt:
            logger.info("Fetching clubb_zt data")
            #data_zt = get_all_variables(nc_zt, plots.lines_zt, plots.sortPlots_zt, len(h), len(t), idx_t0, idx_t1, idx_h0, idx_h1, filler=plots.filler)
            data_zt = get_all_variables(nc_zt, [l[PLOT_LINES] for l in plots.plots_zt], [l[PLOT_ID] for l in plots.plots_zt], len(h), nt, idx_t0, idx_t1, idx_h0, idx_h1, filler=plots.filler)
        # Combine both dictionaries
        data = data_zm.copy()
        data.update(data_zt)
    else:
        logger.info("Fetching sam data")
        #data = get_all_variables(nc_sam, plots.lines, plots.sortPlots, len(h), len(t), idx_t0, idx_t1, idx_h0, idx_h1, filler=plots.filler)
        data = get_all_variables(nc_sam, [l[PLOT_LINES] for l in plots.plots], [l[PLOT_ID] for l in plots.plots], len(h), nt, idx_t0, idx_t1, idx_h0, idx_h1, filler=plots.filler)
    logger.info("Create plots")
    # Center plots if plotting budgets
    #logger.debug("Centering=%s", str('budget' in plots.name))
    plot_default(plots, cf, data, h, centering='budget' in plots.name.lower())
    #logger.debug(mpp.rcParams['text.usetex'])


def plotgen_3d(plots, cf):
    """
    Routine setting up horizontal and cloud conditional plotting from SAM 3D snapshots.
    Data is processed here and afterwards, the plot_3d routine is called to create plots and movies.
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
    idx_t1 = min((np.abs(t - cf.endTime)).argmin()+1, t.size-1)
    nt = idx_t1 - idx_t0
    h = get_h_dim(nc_std, model='sam', create=False, cf=cf)
    # Calculate height level limits
    h_limits = np.array([0])
    logger.debug(h)
    for i in range(len(h)):
        h_limits = np.append(h_limits, 2*h[i]-h_limits[i])
    # Calculate height indices based on parameters in case file
    idx_h0 = (np.abs(h - cf.startHeight)).argmin()
    idx_h1 = min((np.abs(h - cf.endHeight)).argmin()+1, h.size-1)
    logger.debug('Height indices: [%d,%d]',idx_h0, idx_h1)
    # Slice arrays
    h = h[idx_h0:idx_h1]
    h_limits = h_limits[idx_h0:idx_h1+1]
    logger.debug(h_limits)
    # Calculate height level extents
    h_extent = np.diff(h_limits)
    logger.debug(h_extent)
    logger.info("Fetching sam data")
    # TODO: Change structure of 3d setup file as well?
    data_std = get_all_variables(nc_std, plots.lines_std, plots.sortPlots_std, len(h), nt, idx_t0, idx_t1, idx_h0, idx_h1, filler=plots.filler)
    logger.info("Fetching sam_3d data")
    data_3d = get_all_variables(nc_3d, plots.lines_3d, plots.sortPlots_3d, len(h), 1, 0, 0, idx_h0, idx_h1, filler=plots.filler)
    data = {key:value[0][-1]  for (key,value) in data_3d.iteritems()}
    data.update({key:[value[0][0], value[0][-1]] for (key,value) in data_std.iteritems()})
    #logger.debug("3d data: %s", str(data))
    logger.info("Create plots")
    #plot_3d(plots, cf, data, h, h_limits, h_extent, prm_vars, dil_len=dil_len, gif=gif)
    plot_3d(plots, cf, data, h, h_limits, h_extent, dil_len=dil_len, gif=gif)
    
    
def plotgen_comparison(plots, cf):
    """
    Routine setting up height profile comparison plotting for both SAM and CLUBB data.
    Data is processed here and afterwards, the plot_comparison routine is called to create plots
    Added code to include old CLUBB data in comparison plots.
    TODO: Find better way to implement this.
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
        if 'clubb_zt' in plots.nc_files:
            nc_zt_old = ncs.pop()
        if 'clubb_zm' in plots.nc_files:
            nc_zm_old = ncs.pop()
            nc_clubb_old = nc_zm_old
        else:
            nc_clubb_old = nc_zt_old
    # Set variables depending on which nc files are used
    #if plots.sortPlots_zm and plots.sortPlots_zt:
    if all([ncfile in plots.nc_files for ncfile in ['clubb_zm', 'clubb_zt']]):
        nc_zm, nc_zt, nc_sam = ncs
        nc_clubb = nc_zm
    #elif plots.sortPlots_zm:
    elif 'clubb_zm' in plots.nc_files:
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
    idx_h1_clubb = min((np.abs(h_clubb - cf.endHeight)).argmin()+1, h_clubb.size-1)
    # Times for CLUBB are given in seconds, for SAM in minutes
    idx_t0_clubb = (np.abs(t_clubb - cf.startTime*60)).argmin()
    idx_t1_clubb = min((np.abs(t_clubb - cf.endTime*60)).argmin()+1, t_clubb.size-1)
    h_clubb = h_clubb[idx_h0_clubb:idx_h1_clubb]
    logger.info('SAM')
    t_sam = get_t_dim(nc_sam)
    h_sam = get_h_dim(nc_sam, model='sam')
    idx_h0_sam = (np.abs(h_sam - cf.startHeight)).argmin()
    idx_h1_sam = min((np.abs(h_sam - cf.endHeight)).argmin()+1, h_sam.size-1)
    idx_t0_sam = (np.abs(t_sam - cf.startTime)).argmin()
    idx_t1_sam = min((np.abs(t_sam - cf.endTime)).argmin()+1, t_sam.size-1)
    h_sam = h_sam[idx_h0_sam:idx_h1_sam]
    h_clubb_old = None
    if old_clubb:
        t_clubb_old = get_t_dim(nc_clubb_old)
        h_clubb_old = get_h_dim(nc_clubb_old, model='clubb')
        idx_h0_clubb_old = (np.abs(h_clubb_old - cf.startHeight)).argmin()
        idx_h1_clubb_old = min((np.abs(h_clubb_old - cf.endHeight)).argmin()+1, h_clubb_old.size-1)
        # Times for CLUBB are given in seconds, for SAM in minutes
        idx_t0_clubb_old = (np.abs(t_clubb_old - cf.startTime*60)).argmin()
        idx_t1_clubb_old = min((np.abs(t_clubb_old - cf.endTime*60)).argmin()+1, t_clubb_old.size-1)
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
    # TODO: Change to new data structure in variables setup file
    logger.info("Fetching SAM data")
    #data_sam = get_all_variables(nc_sam, plots.lines_sam, plots.sortPlots, len(h_sam), len(t_sam), idx_t0_sam, idx_t1_sam, idx_h0_sam, idx_h1_sam, filler=plots.filler)
    data_sam = get_all_variables(nc_sam, [l[PLOT_LINES_SAM] for l in plots.plots], [l[PLOT_ID] for l in plots.plots], len(h_sam), len(t_sam), idx_t0_sam, idx_t1_sam, idx_h0_sam, idx_h1_sam, filler=plots.filler)
    #logger.debug(data_sam)
    logger.info("Fetching CLUBB data")
    data_zm = None
    data_zt = None
    #if plots.sortPlots_zm:
    if 'clubb_zm' in plots.nc_files:
        logger.info("Fetching clubb_zm data")
        #data_zm = get_all_variables(nc_zm, plots.lines_zm, plots.sortPlots_zm, len(h_clubb), len(t_clubb), idx_t0_clubb, idx_t1_clubb, idx_h0_clubb, idx_h1_clubb, filler=plots.filler)
        data_zm = get_all_variables(nc_zm, [l[PLOT_LINES_CLUBB] for l in plots.plots_zm], [l[PLOT_ID] for l in plots.plots_zm], len(h_clubb), len(t_clubb), idx_t0_clubb, idx_t1_clubb, idx_h0_clubb, idx_h1_clubb, filler=plots.filler)
        #logger.debug(data_clubb_zm)
    #if plots.sortPlots_zt:
    if 'clubb_zt' in plots.nc_files:
        logger.info("Fetching clubb_zt data")
        #data_zt = get_all_variables(nc_zt, plots.lines_zt, plots.sortPlots_zt, len(h_clubb), len(t_clubb), idx_t0_clubb, idx_t1_clubb, idx_h0_clubb, idx_h1_clubb, filler=plots.filler)
        data_zt = get_all_variables(nc_zt, [l[PLOT_LINES_CLUBB] for l in plots.plots_zt], [l[PLOT_ID] for l in plots.plots_zt], len(h_clubb), len(t_clubb), idx_t0_clubb, idx_t1_clubb, idx_h0_clubb, idx_h1_clubb, filler=plots.filler)
        #logger.debug(data_clubb_zt)
    # combine both dictionaries
    data_clubb = data_zm.copy()
    data_clubb.update(data_zt)
    logger.info("Fetching old CLUBB data")
    data_clubb_old = None
    if old_clubb:
        data_zm_old = None
        data_zt_old = None
        #if plots.sortPlots_zm:
        if 'clubb_zm' in plots.nc_files:
            logger.info("Fetching old clubb_zm data")
            #data_zm_old = get_all_variables(nc_zm_old, plots.lines_zm, plots.sortPlots_zm, len(h_clubb_old), len(t_clubb_old), idx_t0_clubb_old, idx_t1_clubb_old, idx_h0_clubb_old, idx_h1_clubb_old, filler=plots.filler)
            data_zm_old = get_all_variables(nc_zm_old, [l[PLOT_LINES_CLUBB] for l in plots.plots_zm], [l[PLOT_ID] for l in plots.plots_zm], len(h_clubb_old), len(t_clubb_old), idx_t0_clubb_old, idx_t1_clubb_old, idx_h0_clubb_old, idx_h1_clubb_old, filler=plots.filler)
            #logger.debug(data_zm_old)
        #if plots.sortPlots_zt:
        if 'clubb_zt' in plots.nc_files:
            logger.info("Fetching old clubb_zt data")
            #data_zt_old = get_all_variables(nc_zt_old, plots.lines_zt, plots.sortPlots_zt, len(h_clubb_old), len(t_clubb_old), idx_t0_clubb_old, idx_t1_clubb_old, idx_h0_clubb_old, idx_h1_clubb_old, filler=plots.filler)
            data_zt_old = get_all_variables(nc_zt_old, [l[PLOT_LINES_CLUBB] for l in plots.plots_zt], [l[PLOT_ID] for l in plots.plots_zt], len(h_clubb_old), len(t_clubb_old), idx_t0_clubb_old, idx_t1_clubb_old, idx_h0_clubb_old, idx_h1_clubb_old, filler=plots.filler)
            #logger.debug(data_zt_old)
        # combine both dictionaries
        data_clubb_old = data_zm_old.copy()
        data_clubb_old.update(data_zt_old)
        #logger.debug(data_clubb_old)
    logger.info("Create plots")
    plot_comparison(plots, cf, data_clubb, data_sam, h_clubb, h_sam, plot_old_clubb=old_clubb, data_old=data_clubb_old, h_old=h_clubb_old)


# Compile plot functions in dict
plot_dict = {
    "budget" : plotgen_default,
    "standalone" : plotgen_default,
    "3d" : plotgen_3d,
    "comparison" : plotgen_comparison,
    }