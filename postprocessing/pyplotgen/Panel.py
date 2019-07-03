'''
:author: Nicolas Strike
:date: Mid 2019
'''
import os
import string

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
            pass
            # self.x_title =
            # self.y_title =
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

    def plot(self, casename):
        '''
        Saves a single panel/graph to the output directory specified by the pyplotgen launch paramters

        The x_values and y_values lists must be in order, as they are associated by index

        :param x_values: A list containing the x values at which each y value shall be plotted
        :param y_values: A list containing the y to plotAll at each x value
        :param title: The title of the panel, e.g. Liquid Water Potential Temperature
        :param x_title: Label for x-axis, e.g. thlm [K]
        :param y_title: Label for y-axis, e.g. Height [m]
        :param line_format: Describes how the line shall appear. See https://matplotlib.org/2.1.2/api/_as_gen/matplotlib.pyplot.plotAll.html
        :return: n/a
        '''

        plt.figure()
        plt.subplot(111)
        for var in self.all_plots:
            x_data = var.x.data
            y_data = var.y.data
            if x_data.shape[0] != y_data.shape[0]:
                raise ValueError("X and Y data have different shapes X: "+str(x_data.shape)
                                 + "  Y:" + str(y_data.shape) + ". Attempted to plotAll " + self.title + " using X: " +
                                 self.x_title + "  Y: " + self.y_title)
            plt.plot(x_data, y_data, var.line_format, label=var.label)
        plt.title(self.title)
        plt.ylabel(self.y_title)
        plt.xlabel(self.x_title)
        plt.figlegend()

        # Create folders
        # Because os.mkdir("output") can fail and prevent os.mkdir("output/" + casename) from being called we must
        # use two separate trys
        try:
            os.mkdir("output")
        except FileExistsError:
            pass # do nothing
        try:
            os.mkdir("output/" + casename)
        except FileExistsError:
            pass # do nothing
        filename = self.title
        filename = filename.translate(str.maketrans('', '', string.punctuation))
        filename = filename.replace(' ', '_')
        rel_filename = "output/" +casename+'/' + filename
        plt.savefig(rel_filename)
