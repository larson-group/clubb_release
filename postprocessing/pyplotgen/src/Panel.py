'''
:author: Nicolas Strike
:date: Mid 2019
'''
import os
from datetime import datetime

import matplotlib.pyplot as plt

from config import Style_definitions


class Panel:
    """
    Represents an individual panel/graph. Each panel contains a number of details
    specific to it, such as a title, axis labels, and lines. Each panel can be plotted/saved to a file.
    """
    TYPE_PROFILE = 'profile'
    TYPE_BUDGET = 'budget'
    TYPE_TIMESERIES = 'timeseries'

    def __init__(self, plots, panel_type="profile", title="Unnamed panel", dependant_title="dependant variable"):
        """
        Creates a new panel
        :param plots: list of Line objects to plot onto the panel
        :param panel_type: Type of panel being plotted (i.e. budget, profile, timeseries)
        :param title: The title of this plot (e.g. 'Liquid water potential tempature')
        :param dependant_title: Label of the dependant axis (labels the x-axis for all panel types except timeseries).
        """

        self.panel_type = panel_type
        self.all_plots = plots
        self.title = title
        self.dependant_title = dependant_title
        self.x_title = "x title unassigned"
        self.y_title = "y title unassigned"
        self.__init_axis_titles__()

    def __init_axis_titles__(self):
        """

        :return:
        """
        if self.panel_type is Panel.TYPE_PROFILE:
            self.x_title = self.dependant_title
            self.y_title = "Height, [m]"
        elif self.panel_type is Panel.TYPE_BUDGET:
            self.y_title = "Height, [m]"
            self.x_title = self.dependant_title
        elif self.panel_type is Panel.TYPE_TIMESERIES:
            self.x_title = "Time [min]"
            self.y_title = self.dependant_title
        else:
            raise ValueError('Invalid panel type ' + self.panel_type + '. Valid options are profile, budget, timeseries')

    def __getStartEndIndex__(self, data, start_value, end_value):
        """
        Get the list floor index that contains the value to start graphing at and the
        ceiling index that contains the end value to stop graphing at

        If neither are found, returns the entire array back
        :param start_value: The first value to be graphed (may return indexes to values smaller than this)
        :param end_value: The last value that needs to be graphed (may return indexes to values larger than this)
        :return: (tuple) start_idx, end_idx   which contains the starting and ending index representing the start and end time passed into the function
        :author: Nicolas Strike
        """
        start_idx = 0
        end_idx = len(data) -1
        for i in range(0,len(data)):
            # Check for start index
            test_value = data[i]
            if test_value <= start_value and test_value > data[start_idx]:
                start_idx = i
            # Check for end index
            if test_value >= end_value and test_value < data[end_idx]:
                end_idx = i

        return start_idx, end_idx

    def plot(self, output_folder, casename, replace_images = False, no_legends = True, thin_lines = False):
        """
         Saves a single panel/graph to the output directory specified by the pyplotgen launch parameters

        :param casename: The name of the case that's plotting this panel
        :return: None
        """
        plt.figure()
        plt.subplot(111)
        for var in self.all_plots:
            x_data = var.x
            y_data = var.y
            if x_data.shape[0] != y_data.shape[0]:
                raise ValueError("X and Y data have different shapes X: "+str(x_data.shape)
                                 + "  Y:" + str(y_data.shape) + ". Attempted to plot " + self.title + " using X: " +
                                 self.x_title + "  Y: " + self.y_title)
            if var.line_format == Style_definitions.LES_LINE_STYLE:
                linewidth = Style_definitions.LES_LINE_THICKNESS
            else:
                linewidth = Style_definitions.DEFAULT_LINE_THICKNESS
            if thin_lines:
                linewidth = Style_definitions.THIN_LINE_THICKNESS
            plt.plot(x_data, y_data, var.line_format, label=var.label, linewidth=linewidth)
        # plt.rcParams.update({'font.size': 16})

        # Set font sizes
        plt.rc('font', size=Style_definitions.DEFAULT_TEXT_SIZE)          # controls default text sizes
        plt.rc('axes', titlesize=Style_definitions.AXES_TITLE_FONT_SIZE)     # fontsize of the axes title
        plt.rc('axes', labelsize=Style_definitions.AXES_LABEL_FONT_SIZE)    # fontsize of the x and y labels
        plt.rc('xtick', labelsize=Style_definitions.X_TICKMARK_FONT_SIZE)    # fontsize of the tick labels
        plt.rc('ytick', labelsize=Style_definitions.Y_TICKMARK_FONT_SIZE)    # fontsize of the tick labels
        plt.rc('legend', fontsize=Style_definitions.LEGEND_FONT_SIZE)    # legend fontsize
        plt.rc('figure', titlesize=Style_definitions.TITLE_TEXT_SIZE)  # fontsize of the figure title

        # Set titles
        plt.title(self.title)
        plt.ylabel(self.y_title)
        plt.xlabel(self.x_title)
        if no_legends is False:
            plt.figlegend()
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

        filename = self.panel_type + "_"+ str(datetime.now().strftime("%Y-%m-%d %H:%M:%S"))

        if self.panel_type == Panel.TYPE_BUDGET:
            filename = filename + "_"+ self.title
        else:
            filename = filename + '_' + self.y_title + "_VS_" + self.x_title
        # filename = filename.translate(str.maketrans('', '', string.punctuation))
        filename = self.__remove_invalid_filename_chars__(filename)
        filename = filename.replace(' ', '_')
        rel_filename = output_folder + "/" +casename+'/' + filename
        if os.path.isfile(rel_filename+'.jpg') and replace_images is True:
            plt.savefig(rel_filename+'.jpg', format='jpeg')
        if not os.path.isfile(rel_filename+'.png'):
            plt.savefig(rel_filename+'.jpg', format='jpeg')
        if os.path.isfile(rel_filename + '.jpg') and replace_images is False:
            print("\n\tImage " + rel_filename+'.jpg already exists. To overwrite this image during runtime pass in the --replace (-r) parameter.')
        plt.close()

    def __remove_invalid_filename_chars__(self, filename):
        """
        Removes characters from a string that are not valid for a filename

        :param filename: Filename string to have characters removed
        :return: a character stripped version of the filename
        """
        filename = filename.replace('.', '')
        # filename = filename.replace(':', '')
        # filename = filename.replace('[', '')
        # filename = filename.replace(']', '')
        filename = filename.replace('/', '')
        filename = filename.replace(',', '')
        return filename
