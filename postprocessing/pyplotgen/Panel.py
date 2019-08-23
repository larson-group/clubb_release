'''
:author: Nicolas Strike
:date: Mid 2019
'''
import os
import string
from time import time

import matplotlib.pyplot as plt


class Panel:
    '''
    Represents an individual panel/graph. Each graph contains a number of details
    specific to it, such as a title. Each graph can be plotted/saved to a file.
    '''
    TYPE_PROFILE = 'profile'
    TYPE_BUDGET = 'budget'
    TYPE_TIMESERIES = 'timeseries'

    def __init__(self, plots, panel_type="profile", title="Unnamed panel", dependant_title="dependant variable"):

        self.panel_type = panel_type
        self.all_plots = plots
        self.title = title
        self.dependant_title = dependant_title
        self.x_title = "x title unassigned"
        self.y_title = "y title unassigned"
        self.__init_axis_titles__()
        # self.blacklisted = blacklisted


    def __init_axis_titles__(self):
        '''

        :return:
        '''
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
        '''
        Get the list floor index that contains the value to start graphing at and the
        ceiling index that contains the end value to stop graphing at

        If neither are found, returns the entire array back
        :param start_value: The first value to be graphed (may return indexes to values smaller than this)
        :param end_value: The last value that needs to be graphed (may return indexes to values larger than this)
        :return: (tuple) start_idx, end_idx   which contains the starting and ending index representing the start and end time passed into the function
        :author: Nicolas Strike
        '''
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
        '''
         Saves a single panel/graph to the output directory specified by the pyplotgen launch paramters

        :param casename:
        :return:
        '''
        plt.figure()
        plt.subplot(111)
        for var in self.all_plots:
            x_data = var.x
            y_data = var.y
            if x_data.shape[0] != y_data.shape[0]:
                raise ValueError("X and Y data have different shapes X: "+str(x_data.shape)
                                 + "  Y:" + str(y_data.shape) + ". Attempted to plot " + self.title + " using X: " +
                                 self.x_title + "  Y: " + self.y_title)
            if var.line_format == 'k-':
                linewidth = 4
            else:
                linewidth = 2
            if thin_lines:
                linewidth = 1
            plt.plot(x_data, y_data, var.line_format, label=var.label, linewidth=linewidth)

        plt.title(self.title)
        plt.ylabel(self.y_title)
        plt.xlabel(self.x_title)
        if no_legends is False:
            plt.figlegend()
        # Create folders
        # Because os.mkdir("output") can fail and prevent os.mkdir("output/" + casename) from being called we must
        # use two separate trys
        try:
            os.mkdir(output_folder)
        except FileExistsError:
            pass # do nothing
        try:
            os.mkdir(output_folder + "/" + casename)
        except FileExistsError:
            pass # do nothing
        if self.panel_type == Panel.TYPE_BUDGET:
            filename = str(time()) + self.panel_type + self.title
        else:
            filename = str(time()) + self.panel_type+ '_' + self.y_title + "_VS_" + self.x_title
        filename = filename.translate(str.maketrans('', '', string.punctuation))
        filename = filename.replace(' ', '_')
        rel_filename = output_folder + "/" +casename+'/' + filename
        if os.path.isfile(rel_filename+'.png') and replace_images is True:
            plt.savefig(rel_filename)
        if not os.path.isfile(rel_filename+'.png'):
            plt.savefig(rel_filename)
        if os.path.isfile(rel_filename + '.png') and replace_images is False:
            print("\n\tImage " + rel_filename+'.png already exists. To overwrite this image during runtime pass in the --replace (-r) parameter.')
        plt.close()
