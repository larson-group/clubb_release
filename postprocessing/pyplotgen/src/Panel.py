'''
:author: Nicolas Strike
:date: Mid 2019
'''
import os
import warnings
from datetime import datetime
from textwrap import fill

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from cycler import cycler

from config import Style_definitions
from src.interoperability import clean_path, clean_title

class Panel:
    """
    Represents an individual panel/graph. Each panel contains a number of details
    specific to it, such as a title, axis labels, and lines. Each panel can be plotted/saved to a file.

    For information on the input parameters of this class, please see the documentation for the
    ``__init__()`` method.
    """
    TYPE_PROFILE = 'profile'
    TYPE_BUDGET = 'budget'
    TYPE_TIMESERIES = 'timeseries'
    TYPE_TIMEHEIGHT = 'timeheight'
    TYPE_ANIMATION = 'animation'
    TYPE_SUBCOLUMN = 'subcolumn'

    VALID_PANEL_TYPES = [TYPE_PROFILE, TYPE_BUDGET, TYPE_TIMESERIES, TYPE_TIMEHEIGHT, TYPE_ANIMATION, TYPE_SUBCOLUMN]

    def __init__(self, plots, bkgrnd_rcm, altitude_bkgrnd_rcm, start_alt_idx, end_alt_idx,
                 panel_type="profile", title="Unnamed panel", dependent_title="dependent variable", sci_scale = None,
                 centered = False, bkgrnd_rcm_flag = False, file_identifier=''):
        """
        Creates a new panel

        :param plots: list of Line objects to plot onto the panel
        :param bkgrnd_rcm: Vertical profile of rcm that can be displayed in the background of plots
                           Time-averaged values for regular profile plots, not time-averaged for animations
        :param altitude_bkgrnd_rcm: Heights corresponding with the background rcm profile.
        :param start_alt_idx: Index of bkgrnd_rcm_tavg that corresponds with bottom of the plot.
        :param end_alt_idx: Index of bkgrnd_rcm_tavg that corresponds with top of the plot.
        :param panel_type: Type of panel being plotted (i.e. budget, profile, timeseries)
        :param title: The title of this plot (e.g. 'Liquid water potential tempature')
        :param dependent_title: Label of the dependent axis (labels the x-axis for all panel types except timeseries).
        :param sci_scale: The scale at which to display the x axis (e.g. to scale to 1e-3 set sci_scale=-3).
            If not specified, the matplotlib default sci scaling will be used.
        :param centered: If True, the Panel will be centered around 0.
            Profile plots are usually centered, while budget plots are not.
        :param bkgrnd_rcm_flag: Show a height-based "contour" plot of time-averaged rcm behind CLUBB profiles.
        """
        self.panel_type = panel_type
        self.all_plots = plots
        self.bkgrnd_rcm = bkgrnd_rcm
        self.altitude_bkgrnd_rcm = altitude_bkgrnd_rcm
        self.start_alt_idx = start_alt_idx
        self.end_alt_idx = end_alt_idx
        self.title = title
        self.dependent_title = dependent_title
        self.x_title = "x title unassigned"
        self.y_title = "y title unassigned"
        self.__init_axis_titles__()
        self.sci_scale = sci_scale
        self.centered = centered
        self.bkgrnd_rcm_flag = bkgrnd_rcm_flag
        self.file_identifier = file_identifier

    def __init_axis_titles__(self):
        """
        Sets the axis titles for this Panel depending on self.panel_type

        :return: None
        """
        if self.panel_type is Panel.TYPE_PROFILE:
            self.x_title = self.dependent_title
            self.y_title = "Height [m]"
        elif self.panel_type is Panel.TYPE_BUDGET:
            self.y_title = "Height [m]"
            self.x_title = self.dependent_title
        elif self.panel_type is Panel.TYPE_SUBCOLUMN:
            self.x_title = self.dependent_title
            self.y_title = "Height [m]"
        elif self.panel_type is Panel.TYPE_TIMESERIES:
            self.x_title = "Time [min]"
            self.y_title = self.dependent_title
        elif self.panel_type == Panel.TYPE_TIMEHEIGHT:
            self.x_title = "Time [min]"
            self.y_title = "Height [m]"
        else:
            raise ValueError('Invalid panel type ' + self.panel_type +
                             '. Valid options are: ' + str(Panel.VALID_PANEL_TYPES))

    def plot(self, output_folder, casename, replace_images = False, no_legends = True, thin_lines = False,
             alphabetic_id="", paired_plots = True, image_extension=".png", 
             generate_grid_adapt_plot=False, grid_comparison_plot=False, read_file_paths=[]):
        """
        Saves a single panel/graph as image to the output directory specified by the pyplotgen launch parameters

        :param output_folder: String containing path to folder in which the image files should be created
        :param casename: The name of the case that is plotted in this panel
        :param replace_images: Switch to tell pyplotgen if existing files should be overwritten
        :param no_legends: If False, a legend will be generated for this Panel
        :param thin_lines: If True, the line_width for this Panel is specified in Style_definitions.THIN_LINE_THICKNESS
        :param alphabetic_id: A string printed into the Panel at coordinates (.9,.9) as an identifier.
        :paired_plots: If no format is specified and paired_plots is True,
            use the color/style rotation specified in Style_definitions.py
        :return: None
        Warning! Argument `replace_images` is unused here!
        """
        # Suppress deprecation warnings
        with warnings.catch_warnings():
            # Create new figure and axis
            warnings.simplefilter("ignore")
            plt.subplot(111)

        # Set line color/style. This will cycle through all colors,
        # then once colors run out use a new style and cycle through
        # colors again
        default_cycler = (
                cycler(linestyle=Style_definitions.STYLE_ROTATION) * cycler(color=Style_definitions.COLOR_ROTATION))
        if not generate_grid_adapt_plot:
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
        plot_dashed = False


        max_panel_value = 0
        for var in self.all_plots:
            legend_char_wrap_length = 17
            legend_char_wrap_length = 100
            #var.label = var.label.replace('_', ' ') # replace _'s in foldernames with spaces for the legend label
            var.label = fill(var.label, width=legend_char_wrap_length)
            x_data = var.x
            if self.sci_scale is not None:
                x_data = x_data * math_scale_factor
            y_data = var.y

            # Find absolutely greatest value in x_data in Panel for centering
            # Suppress "All-NaN slice encountered" warning
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                max_variable_value = max(abs(np.nanmin(x_data)),np.nanmax(x_data))
            max_panel_value = max(max_panel_value,max_variable_value)

            if x_data.shape[0] != y_data.shape[0]:
                raise ValueError("X and Y dependent_data have different shapes X: "+str(x_data.shape)
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
            if generate_grid_adapt_plot:
                read_dir = read_file_paths[0]
                output_dir = output_folder
                file_ending = '_grid_adapt.txt'
                for filename in os.listdir(read_dir):
                    if filename.endswith(file_ending):
                        read_file = read_dir + '/' + filename
                        write_file_plot = output_dir + '/' + filename.split('.')[0] + '.png'

                        with open(read_file, 'r') as file:
                            lines = file.readlines()
                        matrix_grid = []
                        times_grid = []
                        restrict_time_frame = False
                        restrict_shown_grid_levs = False
                        max_time = 800
                        max_grid_levs = 15
                        for line in lines:
                            splitted_line = line.split()
                            is_grid = (splitted_line[0]).strip() == 'g'
                            itime = line.split()[1]
                            itime = float(itime)
                            if itime <= max_time or not restrict_time_frame:
                                if is_grid:
                                    times_grid.append(itime)
                                n = len(splitted_line)
                                grid = splitted_line[2:n]
                                if is_grid:
                                    matrix_grid.append([float(level) for level in grid])
                        times = np.array(times_grid)
                        matrix = np.array(matrix_grid)
                        if restrict_shown_grid_levs:
                            plt.plot(times, matrix[:,0:max_grid_levs])
                        else:
                            plt.plot(times, matrix)
                        #plt.xlabel('time [min]')
                        title = ''
                        if ('arm' in read_file):
                            title = 'ARM'
                        elif ('astex' in read_file):
                            title = 'ASTEX'
                        elif ('gabls2' in read_file):
                            title = 'GABLS2'

                        ax = plt.gca()
                        box = ax.get_position()
                        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
                        plt.title(title)
                        plt.xlabel('iteration')
                        plt.ylabel('z [m]')
                        plt.savefig(write_file_plot, dpi=300, bbox_inches='tight')
                        plt.clf()
                raise Exception('only calculated grid adapt plot, set generate_grid_adapt_plot to False')
            elif grid_comparison_plot:
                hi_res = [i for i in range(0,9001,20)]
                hi_res_grid_spacings = []
                hi_res_grid_spacings_heights = []
                for i in range(len(hi_res)-1):
                    hi_res_grid_spacings.append(hi_res[i+1]-hi_res[i])
                    hi_res_grid_spacings_heights.append((hi_res[i+1]+hi_res[i])/2)

                # This is the path to the file where the dycore grid is stored in, so dycore.grd
                dycore_file_path = '../../input/grid/dycore.grd'

                # Read in dycore grid and make it an array
                dycore = []
                with open(dycore_file_path) as f:
                    for line in f:
                        if (float(line) <= 9000.0):
                            dycore.append(float(line))
                dycore_grid_spacings = []
                dycore_grid_spacings_heights = []
                for i in range(len(dycore)-1):
                    dycore_grid_spacings.append(dycore[i+1]-dycore[i])
                    dycore_grid_spacings_heights.append((dycore[i+1]+dycore[i])/2)

                plt.plot(1, 1)

                ax = plt.gca()
                box = ax.get_position()
                ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

                plt.scatter(hi_res_grid_spacings, hi_res_grid_spacings_heights, s=0)
                for i, (x, y) in enumerate(zip(hi_res_grid_spacings, hi_res_grid_spacings_heights)):
                    plt.text(x, y, str(i+1), fontsize=12, ha='center', va='center', color='black')

                plt.scatter(dycore_grid_spacings, dycore_grid_spacings_heights, s=0)
                for i, (x, y) in enumerate(zip(dycore_grid_spacings, dycore_grid_spacings_heights)):
                    plt.text(x, y, str(i+1), fontsize=12, ha='center', va='center', color='red')

                plt.title('Grid comparison')
                plt.xlabel('$\Delta z \quad [\mathrm{m}]$')
                plt.ylabel('$z \quad [\mathrm{m}]$')

                # Create a custom legend entry
                legend_patch_hi_res = mpatches.Patch(color='black', label="hi-res")
                legend_patch_dycore = mpatches.Patch(color='red', label="dycore")

                # Add the legend
                plt.legend(handles=[legend_patch_hi_res, legend_patch_dycore])

                plt.savefig(output_folder + '/grid_comp.png', dpi=300, bbox_inches='tight')
                plt.close()
                raise Exception('only calculated grid adapt plot, set grid_comparison_plot to False')
                #plt.savefig('grid_comp.png', dpi=300, bbox_inches='tight')
                #print('after savefig')
                #plt.clf()
                #raise Exception('only calculated grid adapt plot, set grid_comparison_plot to False')
            elif plotting_benchmark:
                plt.plot(x_data, y_data, var.line_format, label=var.label, linewidth=line_width)
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

                plt.plot(x_data, y_data, linestyle=line_style, label=var.label, linewidth=line_width)
            else:
                plt.plot(x_data, y_data, label=var.label, linewidth=line_width)

        # Show grid if enabled
        ax = plt.gca()
        ax.grid(Style_definitions.SHOW_GRID)

        ax.set_prop_cycle(default_cycler)

        # Set titles
        plt.title(self.title)
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
            ax.legend(loc='upper right', bbox_to_anchor=(1, 1))
            #if not self.bkgrnd_rcm_flag:
            #    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
            #else:
            #    ax.legend(loc='upper right', bbox_to_anchor=(1, 1))

        # Center plots
        if max_panel_value != 0:
            if self.centered:
                plt.xlim(-1 * max_panel_value * Style_definitions.BUDGET_XAXIS_SCALE_FACTOR,
                         max_panel_value * Style_definitions.BUDGET_XAXIS_SCALE_FACTOR)

        # Emphasize 0 line in profile plots if 0 is in x-axis range
        xlim = plt.xlim()
        if xlim[0] == 0 and xlim[1] == 0:
            plt.xlim=(-1,1)
        if self.panel_type == Panel.TYPE_PROFILE and xlim[0] <= 0 <= xlim[1]:
            plt.axvline(x=0, color='grey', ls='-')

        # Background rcm contour plot
        if self.bkgrnd_rcm_flag:
            num_points = self.end_alt_idx - self.start_alt_idx + 1
            bkgrnd_rcm_contours = np.eye( num_points )
            for i in range(num_points):
                bkgrnd_rcm_contours[i,:] = self.bkgrnd_rcm[self.start_alt_idx+i]
            min_value = min( self.bkgrnd_rcm )
            max_value = max( self.bkgrnd_rcm )
            x_vector_contour = np.zeros( num_points )
            x_diff = xlim[1] - xlim[0]
            x_interval = x_diff / ( num_points - 1 )
            for k in range(num_points):
                x_vector_contour[k] = xlim[0] + float(k) * x_interval
            plt.contourf( x_vector_contour, self.altitude_bkgrnd_rcm[self.start_alt_idx:self.end_alt_idx+1], bkgrnd_rcm_contours,
                          vmin=min_value, vmax=2.0*max_value, cmap=plt.set_cmap("gist_yarg") )
            plt.colorbar( label="rcm [kg/kg]", orientation="vertical" )

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

        # Generate image filename
        #filename = self.panel_type
        if len(self.file_identifier) > 0:
            filename = self.file_identifier + '_' + self.title
        else:
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
        rel_filename = output_folder + "/" +casename+'/' + filename
        rel_filename = clean_path(rel_filename)
        # Save image file
        plt.savefig(rel_filename + image_extension, dpi=Style_definitions.IMG_OUTPUT_DPI, bbox_inches='tight')
        plt.close()

    def __removeInvalidFilenameChars__(self, filename):
        """
        Removes characters from a string that are not valid for a filename

        :param filename: Filename string to have characters removed
        :return: A character stripped version of the filename
        """
        filename = filename.replace('/', '')
        filename = clean_path(filename)
        filename = clean_title(filename)
        return filename