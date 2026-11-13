# -*- coding: utf-8 -*-
'''
Child class derived from Panel.
Takes the same variable definitions as a Panel for profile plots.
But needs the unaveraged netcdf data.
The changed plot routine outputs animations of profile plots instead of simple profile plots.

:author: Benjamin A. Stephens
:date: October 2020
'''

import glob
import os
import shutil
import warnings
from datetime import datetime
from textwrap import fill

import matplotlib.colors as mpc
import matplotlib.pyplot as plt
import numpy as np
from cycler import cycler

from config import Style_definitions
from src.interoperability import clean_path, clean_title

try:
    import cv2  #opencv-python for writing the movies
except ImportError:
    pass

from src.Panel import Panel

class AnimationPanel(Panel):
    """
    AnimationPanel class derived from Panel
    The difference is that the plot routine will generate an animated plot instead of a standard profile plot.

    For information on the input parameters of this class, please see the documentation for the
    ``__init__()`` method.
    """
    def __init__(self, plots, bkgrnd_rcm=None, altitude_bkgrnd_rcm=None, time_bkgrnd_rcm=None, start_alt_idx=None, end_alt_idx=None,
                 panel_type="profile", title="Unnamed panel", dependent_title="dependent variable",
                 sci_scale = None, centered = False, bkgrnd_rcm_flag = False):
        """
        Creates a new panel

        :param plots: list of Line objects to plot onto the panel
        :param panel_type: Type of panel being plotted (i.e. budget, profile, timeseries)
        :param title: The title of this plot (e.g. 'Liquid water potential tempature')
        :param dependent_title: Label of the dependent axis (labels the x-axis for all panel types except timeseries).
        :param sci_scale: The scale at which to display the x axis (e.g. to scale to 1e-3 set sci_scale=-3).
            If not specified, the matplotlib default sci scaling will be used.
        :param centered: If True, the Panel will be centered around 0.
            Profile plots are usually centered, while budget plots are not.
        Warning! Arguments `sci_scale` and `centered` are unused here!
        """

        self.time_bkgrnd_rcm = time_bkgrnd_rcm

        super().__init__(plots, bkgrnd_rcm=bkgrnd_rcm, altitude_bkgrnd_rcm=altitude_bkgrnd_rcm,
                         start_alt_idx=start_alt_idx, end_alt_idx=end_alt_idx, panel_type=panel_type, title=title,
                         dependent_title=dependent_title, sci_scale=None, centered=False, bkgrnd_rcm_flag = bkgrnd_rcm_flag)

    def plot(self, output_folder, casename, replace_images = False, no_legends = True, thin_lines = False,
             alphabetic_id="", paired_plots = True, image_extension=".png", movie_extension=".mp4"):
        """
        New version of plot routine to generate movies of profiles.

        :param output_folder: String containing path to folder in which the image files should be created
        :param casename: The name of the case that is plotted in this panel
        :param replace_images: Switch to tell pyplotgen if existing files should be overwritten
        :param no_legends: If False, a legend will be generated for this Panel
        :param thin_lines: If True, the line_width for this Panel is specified in Style_definitions.THIN_LINE_THICKNESS
        :param alphabetic_id: A string printed into the Panel at coordinates (.9,.9) as an identifier.
        :param paired_plots: If no format is specified and paired_plots is True,
            use the color/style rotation specified in Style_definitions.py
        :param image_extension: Present in case movies can be made from different image types (only .png for now)
        :param movie_extension: Passed so the movies are output to the user's desired format (mp4, avi, etc.)
        :return: None
        Warning! Argument `replace_images` is unused here!
        """
        tmax = np.inf ; idx = 0
        sim_lengths=[]
        for i,var in enumerate(self.all_plots):
            sim_lengths.append(len(var.x))
            if len(var.x) < tmax:
                x_dataset=var.x
                tmax=len(x_dataset)
                idx_max = idx
            idx+=1

        #if discrepanies in number of time steps, filter out the
        #extraneous time steps from those simulations that have extra steps
        if np.all(sim_lengths==tmax):
            filteringFlag=False
        else:
            filteringFlag=True
            idx = 0
            for var in self.all_plots:
                if idx == idx_max:
                    idx+=1
                else:
                    temp_x_data=var.x
                    temp_y_data=var.y
                    temp_data=var.data
                    filtered_data = np.zeros((len(x_dataset),len(temp_y_data)))
                    for i in range(len(x_dataset)):
                        for j in range(len(temp_x_data)):
                            if x_dataset[i] == temp_x_data[j] or abs(x_dataset[i]-temp_x_data[j]) < 0.5:
                                filtered_data[i,:]=temp_data[j,:]
                    var.data = filtered_data
                    idx+=1

        if self.bkgrnd_rcm_flag:
            ## Set up rcm contours
            # Find matching rcm times
            rcm_time_idcs = []
            for x in x_dataset*60:
                rcm_time_idcs.append(np.argmin(np.abs(self.time_bkgrnd_rcm-x)))
            # Set up colorbar that is constant over all time frames
            min_rcm_value = self.bkgrnd_rcm.min()
            max_rcm_value = self.bkgrnd_rcm.max()
            cmap = plt.get_cmap("gist_yarg")
            # Modify color map to only use grays
            cmap = mpc.LinearSegmentedColormap.from_list('Only grays', cmap(np.linspace(0.0, 0.5, 256)))
            norm = mpc.Normalize(vmin=min_rcm_value, vmax=max_rcm_value)
            cont_map = plt.cm.ScalarMappable(norm=norm, cmap=cmap)

        min_x_value = np.inf ; max_x_value = -1*np.inf  #set large to be overwritten during first pass below
        for t in range(0,tmax):

            # Suppress deprecation warnings
            with warnings.catch_warnings():
                # Create new figure and axis
                warnings.simplefilter("ignore")
                plt.subplot(111)

            # Set line color/style. This will cycle through all colors,
            # then once colors run out use a new style and cycle through colors again
            default_cycler = ( cycler(linestyle=Style_definitions.STYLE_ROTATION)
                               * cycler(color=Style_definitions.COLOR_ROTATION))
            plt.rc('axes', prop_cycle=default_cycler)

            # Set graph size
            plt.figure(figsize=Style_definitions.FIGSIZE)

            # Set font sizes
            plt.rc('font', size=Style_definitions.DEFAULT_TEXT_SIZE)          # controls default text sizes
            plt.rc('axes', titlesize=Style_definitions.AXES_TITLE_FONT_SIZE)     # fontsize of the axes title
            plt.rc('axes', labelsize=Style_definitions.AXES_LABEL_FONT_SIZE)    # fontsize of the x and y labels
            plt.rc('xtick', labelsize=Style_definitions.X_TICKMARK_FONT_SIZE)    # fontsize of the tick labels
            plt.rc('ytick', labelsize=Style_definitions.Y_TICKMARK_FONT_SIZE)    # fontsize of the tick labels
            plt.rc('legend', fontsize=Style_definitions.LEGEND_FONT_SIZE)    # legend fontsize
            plt.rc('figure', titlesize=Style_definitions.TITLE_TEXT_SIZE)  # fontsize of the figure title

            label_scale_factor = ""
            # Use custom sci scaling
            if self.sci_scale is not None:
                scalepower = -1 * self.sci_scale
                if self.sci_scale != 0:
                    label_scale_factor = "x 1e" + str(self.sci_scale)
                math_scale_factor =  10 ** (scalepower)
                plt.ticklabel_format(style='plain', axis='x')
            # Use pyplot's default sci scaling
            else:
                plt.ticklabel_format(style='sci', axis='x', scilimits=Style_definitions.POW_LIMS)

            # Prevent x-axis label from getting cut off
            plt.gcf().subplots_adjust(bottom=0.15)

            # Plot dashed line. This var will oscillate between true and false
            plot_dashed = True

            for var in self.all_plots:
                legend_char_wrap_length = 17
                var.label = var.label.replace('_', ' ') # replace _'s in foldernames with spaces for the legend label
                var.label = fill(var.label, width=legend_char_wrap_length)
                c_data = var.data
                if self.sci_scale is not None:
                    c_data[t,:] = c_data[t,:] * math_scale_factor
                y_data = var.y

                # Find min/max values for fixed x-axis
                # Suppress "All-NaN slice encountered" warning
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    # Exclude first frame from finding min/max since it
                    # may have some very extreme values
                    current_min=np.nanmin(np.ndarray.flatten(c_data[1:,:]))
                    current_max=np.nanmax(np.ndarray.flatten(c_data[1:,:]))
                    if current_min < min_x_value:
                        min_x_value = current_min
                    elif np.isnan(current_min):
                        min_x_value = -1.0
                    if current_max > max_x_value:
                        max_x_value = current_max
                    elif np.isnan(current_max):
                        max_x_value = 1.0

                #shape sanity check
                if c_data[t,:].shape[0] != y_data.shape[0]:
                    raise ValueError("X and Y dependent_data have different shapes X: "+str(c_data[t,:].shape)
                                     + "  Y:" + str(y_data.shape) + ". Attempted to plot " + self.title +
                                     " using X: " + self.x_title + "  Y: " + self.y_title)

                # Set correct line formatting and plot data
                if var.line_format == Style_definitions.BENCHMARK_LINE_STYLES['coamps']:
                    line_width = Style_definitions.LES_LINE_THICKNESS
                elif var.line_format == Style_definitions.BENCHMARK_LINE_STYLES['sam']:
                    line_width = Style_definitions.LES_LINE_THICKNESS
                elif var.line_format == Style_definitions.BENCHMARK_LINE_STYLES['r408']:
                    line_width = Style_definitions.ARCHIVED_CLUBB_LINE_THICKNESS
                elif var.line_format == Style_definitions.BENCHMARK_LINE_STYLES['e3sm']:
                    line_width = Style_definitions.E3SM_LINE_THICKNESS
                else:
                    line_width = Style_definitions.CLUBB_LINE_THICKNESS
                if thin_lines:
                    line_width = Style_definitions.THIN_LINE_THICKNESS
                plotting_benchmark = var.line_format != ""
                if plotting_benchmark:
                    plt.plot(c_data[t,:], y_data, var.line_format, label=var.label, linewidth=line_width)
                    # If a benchmark defines a custom color (e.g. "gray" or "#404040) this messes up the color rotation.
                    # Setting the prop cycle to None and then redefining it fixes the color rotation.
                    # This fix may be dependent on benchmarks being plotted first. If this stops being the case, colors may
                    # repeat themselves sooner than expected.
                    plt.gca().set_prop_cycle(None)
                    plt.rc('axes', prop_cycle=default_cycler)
                # If format is not specified and paired_plots are enabled,
                # use the color/style rotation specified in Style_definitions.py
                elif paired_plots:
                    if plot_dashed:
                        line_width = Style_definitions.DASHED_LINE_THICKNESS
                        line_style = '--'
                        plot_dashed = False
                    else:
                        line_width = Style_definitions.FLAT_LINE_THICKNESS
                        line_style = '-'
                        plot_dashed = True

                    plt.plot(c_data[t,:], y_data, linestyle=line_style, label=var.label, linewidth=line_width)
                else:
                    plt.plot(c_data[t,:], y_data, label=var.label, linewidth=line_width)

            # Show grid if enabled
            ax = plt.gca()
            ax.grid(Style_definitions.SHOW_GRID)

            ax.set_prop_cycle(default_cycler)

            # Set titles---top title includes minute counter for reference
            plt.title(self.title +'\nMinute = {}'.format(int(x_dataset[t])))
            plt.ylabel(self.y_title)
            plt.text(1, -0.15, label_scale_factor, transform=ax.transAxes, fontsize=Style_definitions.MEDIUM_FONT_SIZE)
            plt.xlabel(self.x_title)

            # Add alphabetic ID
            if alphabetic_id != "":
                ax.text(0.9, 0.9, '('+alphabetic_id+')', ha='center', va='center', transform=ax.transAxes,
                        fontsize=Style_definitions.LARGE_FONT_SIZE) # Add letter label to panels

            # Plot legend
            if no_legends is False:
                # Shrink current axis by 20%
                box = ax.get_position()
                ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
                # Put a legend to the right of the current axis
                if not self.bkgrnd_rcm_flag:
                    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
                else:
                    ax.legend(loc='upper right', bbox_to_anchor=(1, 1))

            # Fix x-axis
            if min_x_value != 0 and max_x_value != 0:
                plt.xlim( min_x_value - abs(min_x_value) * Style_definitions.MOVIE_XAXIS_SCALE_FACTOR,
                          max_x_value + abs(max_x_value) * Style_definitions.MOVIE_XAXIS_SCALE_FACTOR)
            if min_x_value == 0 and max_x_value == 0:
                plt.xlim(-1,1)
            elif min_x_value == 0:
                plt.xlim( min_x_value - abs(max_x_value) * Style_definitions.MOVIE_XAXIS_SCALE_FACTOR,
                          max_x_value + abs(max_x_value) * Style_definitions.MOVIE_XAXIS_SCALE_FACTOR)
            elif max_x_value == 0:
                plt.xlim( min_x_value - abs(min_x_value) * Style_definitions.MOVIE_XAXIS_SCALE_FACTOR,
                          max_x_value + abs(min_x_value) * Style_definitions.MOVIE_XAXIS_SCALE_FACTOR)

            # Emphasize 0 line in profile plots if 0 is in x-axis range
            xlim = plt.xlim()
            if xlim[0] <= 0 <= xlim[1]:
                plt.axvline(x=0, color='grey', ls='-')

            # Background rcm contour plot
            if self.bkgrnd_rcm_flag:
                num_points = self.end_alt_idx - self.start_alt_idx + 1
                bkgrnd_rcm_contours = np.eye( num_points )
                for i in range(num_points):
                    bkgrnd_rcm_contours[i,:] = self.bkgrnd_rcm[rcm_time_idcs[t], self.start_alt_idx+i]
                x_vector_contour = np.zeros( num_points )
                x_diff = xlim[1] - xlim[0]
                x_interval = x_diff / ( num_points - 1 )
                for k in range(num_points):
                    x_vector_contour[k] = xlim[0] + float(k) * x_interval
                plt.contourf( x_vector_contour, self.altitude_bkgrnd_rcm[self.start_alt_idx:self.end_alt_idx+1],
                              bkgrnd_rcm_contours, norm=norm, cmap=cmap )
                plt.gcf().colorbar( cont_map, ax=ax, cmap=cmap, norm=norm, label="rcm [kg/kg]", orientation="vertical" )

            # Create folders
            # Because os.mkdir("output") can fail and prevent os.mkdir("output/" + casename) from being called we must
            # use two separate try blocks
            try:
                os.mkdir(output_folder)
            except FileExistsError:
                pass # do nothing
            try:
                os.mkdir(output_folder + "/" + casename)
            except FileExistsError:
                pass # do nothing

            #create tmp folder to house original images for movie
            temp_dir='tmp'
            try:
                os.mkdir(output_folder + "/" + casename + "/" + temp_dir)
            except FileExistsError:
                pass #do nothing

            # Generate image filename
            filename = self.panel_type + "_"+ str(datetime.now())
            # Force subcolumn plots to show up on top
            if self.panel_type == Panel.TYPE_SUBCOLUMN:
                filename = 'aaa' + filename
            if self.panel_type == Panel.TYPE_BUDGET:
                filename = filename + "_"+ self.title
            else:
                filename = filename + '_' + self.y_title + "_VS_" + self.x_title
            filename = self.__removeInvalidFilenameChars__(filename)
            # Concatenate with output foldername
            rel_filename = output_folder + "/" + casename + '/' + temp_dir + '/' + filename
            rel_filename = clean_path(rel_filename)
            # Save image file
            plt.savefig(rel_filename+image_extension, dpi=Style_definitions.IMG_OUTPUT_DPI, bbox_inches='tight')
            plt.close()

        # Lights, camera, action!
        img_array=[]
        for filename2 in sorted(glob.glob(output_folder + '/' + casename + '/' + temp_dir + '/*'+image_extension)):
            img = cv2.imread(filename2)
            height, width, layers = img.shape
            size = (width,height)
            img_array.append(img)

        # Set proper codec
        if movie_extension == ".mp4":
            fourcc = cv2.VideoWriter_fourcc(*'mp4v')
        elif movie_extension == ".avi":
            fourcc = cv2.VideoWriter_fourcc(*'XVID')

        # check to see if FFMPEG is available--determines if we need a temporary name
        if shutil.which('ffmpeg') is not None and movie_extension == '.mp4':
            moviename = output_folder + '/' + casename + '/' + "movie" + movie_extension
        else:
            moviename = output_folder + '/' + casename + '/' + filename + movie_extension

        # We're rolling...
        out = cv2.VideoWriter(moviename, fourcc, Style_definitions.FRAMES_PER_SECOND , size)
        for i in range(len(img_array)):
            out.write(img_array[i])

        # Cut!
        out.release()

        # rename if FFMPEG is present
        if shutil.which('ffmpeg') is not None and movie_extension == '.mp4':
            final_name = output_folder + '/' + casename + '/' + filename + movie_extension
            command="ffmpeg -hide_banner -loglevel panic -i " + moviename + " -vcodec libx264 " + final_name
            os.system(command)
            os.remove(moviename)

        # Delete temp folder
        shutil.rmtree(output_folder + "/" + casename + "/" + temp_dir)

        return filteringFlag