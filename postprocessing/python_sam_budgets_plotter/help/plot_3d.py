"""
This module defines functions to generate plane and cloud conditional plots from 3d SAM data
"""

#--------------------
# IMPORTS
#--------------------
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

def generate_cmap(data, cmapName, maxCol='g', minCol='r', zeroCol='w', N=1000, quant=0, interpol=None):
    """
    Function generating segmented colormap from data array
    mapping the minimum value to minCol and the maximum value to maxCol
    Input:
    data - Array of values to be used as basis for colormap
    cmapName - String used to register colormap
    maxCol - Color onto which the maximum value in data is mapped
    minCol - Color onto which the minimum value in data is mapped
    zeroCol - Color onto which zero in data is mapped
    N - Parameter for matplotlib.colors.LinearSegmentedColormap, number of rgb quantization levels
    interpol - Parameter for numpy.nanpercentile specifying the methods of interpolation to find percentile value
    Output:
    matplotlib.colors.LinearSegmentedColormap(cmapName,cdict,N=N)
    """
    # Convert colors to rgba
    maxRGBA = mpl.colors.to_rgba(maxCol)
    minRGBA = mpl.colors.to_rgba(minCol)
    zeroRGBA = mpl.colors.to_rgba(zeroCol)
    
    # Get limits from data in order to generate colormap
    #lims = data.min(), data.max() # Do not use this, because of outliers
    # Get norm boundaries for color scale by discarding outliers only
    #dataSort = np.sort(data.flatten())
    #lims = dataSort[trash], dataSort[-trash]
    # Amount of outliers seeems to be constant with data array size, but to achieve good distribution over color scale, percentiles are better
    try:
        interpol[1]
    except Exception as e:
        interpol = (interpol, interpol)
    lims = np.nanpercentile(data,quant,interpolation=interpol[0]), np.nanpercentile(data,100-quant,interpolation=interpol[1])
    #logger.debug('uw.min=%f, uw.max=%f',np.nanmin(uw),np.nanmax(uw))
    #logger.debug("uwlim=(%f,%f)",uwlim[0], uwlim[1])
    #logger.debug('uw[0].min=%f,uw[0].max=%f',np.nanmin(uw[0]), np.nanmax(uw[0]))
    
    # Define norm mapping values t [0;1]
    norm = mpl.colors.Normalize(vmin=lims[0], vmax=lims[1])
    
    # Calculate mapping of 0 in normalized data
    zero = (0-lims[0])/(lims[1]-lims[0])
    
    # Define color map segments
    cdict = {"red":((0,minRGBA[0],minRGBA[0]),(zero,zeroRGBA[0],zeroRGBA[0]),(1,maxRGBA[0],maxRGBA[0])),
             "green":((0,minRGBA[1],minRGBA[1]),(zero,zeroRGBA[1],zeroRGBA[1]),(1,maxRGBA[1],maxRGBA[1])),
             "blue":((0,minRGBA[2],minRGBA[2]),(zero,zeroRGBA[2],zeroRGBA[2]),(1,maxRGBA[2],maxRGBA[2])),
             "alpha":((0,minRGBA[3],minRGBA[3]),(zero,zeroRGBA[3],zeroRGBA[3]),(1,maxRGBA[3],maxRGBA[3]))
             }
    # Create and return color map
    return mpl.colors.LinearSegmentedColormap(cmapName, cdict, N=N)

def generate_meshgrids(d, skip):
    """
    Function generating 2 meshgrids.
    The first one is a full integer meshgrid for square grid of edge length d
    The second one is an integer meshgrid generated from a 1d range array starting at 0 and going up to d, but skipping skip elements
    """
    grid1d = np.arange(d)
    fullGrid= np.meshgrid(grid1d, grid1d)
    skipGrid = np.meshgrid(grid1d[::skip], grid1d[::skip])
    return (fullGrid, skipGrid)


#def plot_3d(plots, cf, data, h, h_limits, h_extent, prm_vars, fps=2, dil_len=1, gif=False):
def plot_3d(plots, cf, data, h, h_limits, h_extent, fps=2, dil_len=1, gif=False):
    """
    Generate horizontal output showing cloud outlines and wind speeds
    TODO: Special handling of RICO, because of huge grid -> split up grid in multiple subgrids?
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
    quiver_cmap_uw = generate_cmap(uw, 'quiver_cmap_uw', green, red, zero_col, 1000, quant, ('lower','higher'))

    ## VW
    quiver_cmap_vw = generate_cmap(vw, 'quiver_cmap_vw', green, red, zero_col, 1000, quant, ('lower','higher'))

    ## UP
    quiver_cmap_up = generate_cmap(up, 'quiver_cmap_up', green, red, zero_col, 1000, 0, ('lower','higher'))

    ## VP
    quiver_cmap_vp = generate_cmap(vp, 'quiver_cmap_vp', green, red, zero_col, 1000, 0, ('lower','higher'))

    ## W
    quiver_cmap_w = generate_cmap(w, 'quiver_cmap_w', green, red, zero_col, 1000, 0, ('lower','higher'))

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
        #ax1.set_xlabel(r'Eastward grid dimension $\mathrm{{\left[{:.0f} m\right]}}$'.format(prm_vars['dx']), fontsize=fontsizes['labels'])
        #ax1.set_ylabel(r'Northward grid dimension $\mathrm{{\left[{:.0f} m\right]}}$'.format(prm_vars['dx']), fontsize=fontsizes['labels'])
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
        pb.plot_profiles([uw_plots[k]], h[:cutoff], r"Cloud Conditional $\mathrm{\overline{u'w'}\ \left[\frac{m^2}{s^2}\right]}$", cf.yLabel, uw_plots[k][0], os.path.join(jpg_dir,'uw_cld{}'.format(k)), startLevel=0, lw=cf.lw, grid=False, centering=False, pdf=uw_pdf)

    pb.plot_profiles([["u'^2",True,'dummy',0,np.nanmean(up*up, axis=(1,2))]], h, "u'^2", cf.yLabel, "u'^2", os.path.join(jpg_dir,'up2'), startLevel=0, lw=cf.lw, grid=False, centering=False, pdf=uw_pdf)

    pb.plot_profiles(uw_plots, h[:cutoff], r"Cloud conditional $\mathrm{\overline{u'w'}\ \left[\frac{m^2}{s^2}\right]}$", cf.yLabel, r"Conditional means of $\mathrm{\overline{u'w'}}$", os.path.join(jpg_dir,'uw_cld_all'), startLevel=0, lw=cf.lw, grid=False, centering=True, pdf=uw_pdf)

    # TODO: Add standard mean u
    #pb.plot_profiles([["Extended in-cloud mean of u",True,0,u_cld_mean]], h, r"Extended Cloud Conditional $\mathrm{\bar{u}\ \left[\frac{m}{s}\right]}$", cf.yLabel, "Extended Conditional mean of u", os.path.join(jpg_dir,'u_cld_mean'), startLevel=0, lw=cf.lw, grid=False, centering=False, pdf=uw_pdf)

    # U plots
    pb.plot_profiles(um_plots, h[:cutoff], r"Cloud conditional $\mathrm{\bar{u}\ \left[\frac{m}{s}\right]}$", cf.yLabel, "Comparison of conditional means of u", os.path.join(jpg_dir, 'u_cld_all'), startLevel=0, lw=cf.lw, grid=False, centering=False, pdf=uw_pdf)
                                                            
    # VW plots
    for k in range(len(vw_plots)):
        pb.plot_profiles([vw_plots[k]], h[:cutoff], r"Cloud conditional $\mathrm{\overline{v'w'}\ \left[\frac{m^2}{s^2}\right]}$", cf.yLabel, vw_plots[k][0], os.path.join(jpg_dir,'vw_cld{}'.format(k)), startLevel=0, lw=cf.lw, grid=False, centering=False, pdf=uw_pdf)
        
    pb.plot_profiles([["v'^2",True,'dummy', 0,np.nanmean(vp*vp, axis=(1,2))]], h, "v'^2", cf.yLabel, "v'^2", os.path.join(jpg_dir,'vp2'), startLevel=0, lw=cf.lw, grid=False, centering=False, pdf=uw_pdf)

    pb.plot_profiles(vw_plots, h[:cutoff], r"Cloud conditional $\mathrm{\overline{v'w'}\ \left[\frac{m^2}{s^2}\right]}$", cf.yLabel, r"Conditional means of $\mathrm{\overline{v'w'}}$", os.path.join(jpg_dir,'vw_cld_all'), startLevel=0, lw=cf.lw, grid=False, centering=True, pdf=uw_pdf)

    pb.plot_profiles(wq_plots, h[:cutoff], r"Cloud conditional $\mathrm{\overline{w'q_n'}\ \left[g\,m\,kg^{-1}\,s^{-1}\right]}$", cf.yLabel, r"Conditional means of $\overline{w'q_n'}}$", os.path.join(jpg_dir,'wq_cld_all'), startLevel=0, lw=cf.lw, grid=False, centering=True, pdf=uw_pdf)
                                                                                    
    # TODO: Add standard mean v
    #pb.plot_profiles([["Extended in-cloud mean of v",True,0,v_cld_mean]], h, r"Extended Cloud Conditional $\mathrm{\bar{v}\ \left[\frac{m}{s}\right]}$", cf.yLabel, "Extended Conditional mean of v", os.path.join(jpg_dir,'v_cld_mean'), startLevel=0, lw=cf.lw, grid=False, centering=False, pdf=uw_pdf)

    # V plots
    pb.plot_profiles(vm_plots, h[:cutoff], r"Cloud conditional $\mathrm{\bar{v}\ \left[\frac{m}{s}\right]}$", cf.yLabel, "Comparison of conditional means of v", os.path.join(jpg_dir, 'v_cld_all'), startLevel=0, lw=cf.lw, grid=False, centering=False, pdf=uw_pdf)
                                                            
    # W plots
    pb.plot_profiles(wm_plots, h[:cutoff], r"Cloud conditional $\mathrm{\bar{w}\ \left[\frac{m}{s}\right]}$", cf.yLabel, "Comparison of conditional means of w", os.path.join(jpg_dir, 'w_cld_all'), startLevel=0, lw=cf.lw, grid=False, centering=False, pdf=uw_pdf)

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
        #'dx = {}'.format(prm_vars['dx']),
        #'dy = {}'.format(prm_vars['dx']),
        'dx = {}'.format(cf.dxy),
        'dy = {}'.format(cf.dxy),
        'dz = {}'.format(cf.dz),
        'height = {} - {}'.format(cf.startHeight, cf.endHeight),
        '# plotting parameters:',
        'number of arrows skipped = {}'.format(skip),
        'cloud recognition limit for water vapor = {}'.format(cld_lim),
        'background interpolation method = {}'.format(cloud_interpolation),
        'cloud halo dilation length'.format(dil_len),
        'convariance percentile cutoff = {}'.format(quant)
        ]
    param_file.write('\n'.join(prms))
    param_file.close()