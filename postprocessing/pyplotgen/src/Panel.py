'''
:author: Nicolas Strike
:date: Mid 2019
'''
import os
from datetime import datetime

import matplotlib.pyplot as plt
import numpy as np
from cycler import cycler

from config import Style_definitions
from src.interoperability import clean_path, clean_title


class Panel:
    """
    Represents an individual panel/graph. Each panel contains a number of details
    specific to it, such as a title, axis labels, and lines. Each panel can be plotted/saved to a file.
    """
    #x_title: str
    TYPE_PROFILE = 'profile'
    TYPE_BUDGET = 'budget'
    TYPE_TIMESERIES = 'timeseries'
    EXTENSION = '.png'

    def __init__(self, plots, panel_type="profile", title="Unnamed panel", dependent_title="dependent variable",
                 sci_scale = None, centered = False):
        """
        Creates a new panel
        :param plots: list of Line objects to plot onto the panel
        :param panel_type: Type of panel being plotted (i.e. budget, profile, timeseries)
        :param title: The title of this plot (e.g. 'Liquid water potential tempature')
        :param dependent_title: Label of the dependent axis (labels the x-axis for all panel types except timeseries).
        :param sci_scale: The scale at which to display the x axis (e.g. to scale to 1e-3 set sci_scale=-3).
            If not specified, the matplotlib default sci scaling will be used.
        """

        self.panel_type = panel_type
        self.all_plots = plots
        self.title = title
        self.dependent_title = dependent_title
        self.x_title = "x title unassigned"
        self.y_title = "y title unassigned"
        self.__init_axis_titles__()
        self.sci_scale = sci_scale
        self.centered = centered

    def __init_axis_titles__(self):
        """

        :return:
        """
        if self.panel_type is Panel.TYPE_PROFILE:
            self.x_title = self.dependent_title
            self.y_title = "Height [m]"
        elif self.panel_type is Panel.TYPE_BUDGET:
            self.y_title = "Height [m]"
            self.x_title = self.dependent_title
        elif self.panel_type is Panel.TYPE_TIMESERIES:
            self.x_title = "Time [min]"
            self.y_title = self.dependent_title
        else:
            raise ValueError('Invalid panel type ' + self.panel_type +
                             '. Valid options are profile, budget, timeseries')

    def plot(self, output_folder, casename, replace_images = False, no_legends = True, thin_lines = False,
             alphabetic_id="", paired_plots = True):
        """
         Saves a single panel/graph to the output directory specified by the pyplotgen launch parameters

        :param casename: The name of the case that's plotted in this panel
        :param replace_images: Switch to tell pyplotgen if existing files should be overwritten
        :param no_legends: 
        :param thin_lines:
        :param alphabetic_id:
        :paired_plots:
        :return: None
        """
        print('Plotting panel {}'.format(self.title))
        plt.subplot(111)

        # Set line color/style. This will cycle through all colors,
        # then once colors run out use a new style and cycle through
        # colors again
        default_cycler = (
                cycler(linestyle=Style_definitions.STYLE_ROTATION) * cycler(color=Style_definitions.COLOR_ROTATION))
        plt.rc('axes', prop_cycle=default_cycler)

        # Set graph size
        plt.figure(figsize=(10,6))

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
                label_scale_factor = "\t1e" + str(self.sci_scale)
            math_scale_factor =  10 ** (scalepower)
            plt.ticklabel_format(style='plain', axis='x')
        # Use pyplot's default sci scaling
        else:
            plt.ticklabel_format(style='sci', axis='x', scilimits=Style_definitions.POW_LIMS)

        # prevent x-axis label from getting cut off
        plt.gcf().subplots_adjust(bottom=0.15)

        # Plot dashed line. This var will oscillate between true and false
        plot_dashed = True

        max_panel_value = 0
        for var in self.all_plots:
            x_data = var.x
            if self.sci_scale is not None:
                x_data = x_data * math_scale_factor
            y_data = var.y

            max_variable_value = max(abs(np.nanmin(x_data)),np.nanmax(x_data))

            max_panel_value = max(max_panel_value,max_variable_value)

            if x_data.shape[0] != y_data.shape[0]:
                raise ValueError("X and Y dependent_data have different shapes X: "+str(x_data.shape)
                                 + "  Y:" + str(y_data.shape) + ". Attempted to plot " + self.title + " using X: " +
                                 self.x_title + "  Y: " + self.y_title)
            if var.line_format == Style_definitions.BENCHMARK_LINE_STYLES['sam']:
                linewidth = Style_definitions.LES_LINE_THICKNESS
            elif var.line_format == Style_definitions.BENCHMARK_LINE_STYLES['r408']:
                linewidth = Style_definitions.ARCHIVED_CLUBB_LINE_THICKNESS
            elif var.line_format == Style_definitions.BENCHMARK_LINE_STYLES['e3sm']:
                linewidth = Style_definitions.E3SM_LINE_THICKNESS
            else:
                linewidth = Style_definitions.CLUBB_LINE_THICKNESS
            if thin_lines:
                linewidth = Style_definitions.THIN_LINE_THICKNESS
            if var.line_format != "":
                plt.plot(x_data, y_data, var.line_format, label=var.label, linewidth=linewidth)
            # If format is not specified and paired_plots are enabled,
            # use the color/style rotation specified in Style_definitions.py
            elif paired_plots:
                if plot_dashed:
                    linewidth = Style_definitions.DASHED_LINE_THICKNESS
                    linestyle = '--'
                    plot_dashed = False
                else:
                    linewidth = Style_definitions.FLAT_LINE_THICKNESS
                    linestyle = '-'
                    plot_dashed = True

                plt.plot(x_data, y_data, linestyle=linestyle, label=var.label, linewidth=linewidth)
            else:
                plt.plot(x_data, y_data, label=var.label, linewidth=linewidth)

        # Set titles
        plt.title(self.title)
        plt.ylabel(self.y_title)
        plt.xlabel(self.x_title + label_scale_factor)

        # Show grid if enabled
        ax = plt.gca()
        ax.grid(Style_definitions.SHOW_GRID)

        ax.set_prop_cycle(default_cycler)

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
            ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

        # Center plots
        if self.centered:
            plt.xlim(-1 * max_panel_value * Style_definitions.BUDGET_XAXIS_SCALE_FACTOR,
                     max_panel_value * Style_definitions.BUDGET_XAXIS_SCALE_FACTOR)

        # Emphasize 0 line in profile plots if 0 is in x-axis range
        xlim = plt.xlim()
        if self.panel_type == Panel.TYPE_PROFILE and 0 >= xlim[0] and 0 <= xlim[1]:
            plt.axvline(x=0, color='grey', ls='-')

        # Create folders
        # Because os.mkdir("output") can fail and prevent os.mkdir("output/" + casename) from being called we must
        # use two separate try blokcs
        try:
            os.mkdir(output_folder)
        except FileExistsError:
            pass # do nothing
        try:
            os.mkdir(output_folder + "/" + casename)
        except FileExistsError:
            pass # do nothing

        filename = self.panel_type + "_"+ str(datetime.now())
        
        if self.panel_type == Panel.TYPE_BUDGET:
            filename = filename + "_"+ self.title
        else:
            filename = filename + '_' + self.y_title + "_VS_" + self.x_title
        filename = self.__remove_invalid_filename_chars__(filename)
        rel_filename = output_folder + "/" +casename+'/' + filename
        rel_filename = clean_path(rel_filename)
        if replace_images is True or not os.path.isfile(rel_filename+Panel.EXTENSION):
            plt.savefig(rel_filename+Panel.EXTENSION)
        else: # os.path.isfile(rel_filename + Panel.EXTENSION) and replace_images is False:
            print("\n\tImage " + rel_filename+Panel.EXTENSION+
                  ' already exists. To overwrite this image during runtime pass in the --replace (-r) parameter.')
        plt.close()

    def __remove_invalid_filename_chars__(self, filename):
        """
        Removes characters from a string that are not valid for a filename

        :param filename: Filename string to have characters removed
        :return: a character stripped version of the filename
        """
        filename = filename.replace('/', '')
        filename = clean_path(filename)
        filename = clean_title(filename)
        return filename
