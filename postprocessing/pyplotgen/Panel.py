'''
:author: Nicolas Strike
:date: Mid 2019
'''
from pyplotgen.DataReader import NetCdfVariable
from pyplotgen.Lineplot import Lineplot
from pyplotgen.PanelType import PanelType
from pyplotgen.Plotter import Plotter
from pyplotgen.VarnameConversions import CLUBB_TO_SAM


class Panel:
    '''
    Represents an individual panel/graph. Each graph contains a number of details
    specific to it, such as a title. Each graph can be plotted/saved to a file.
    '''
    TYPE_PROFILE = 'profile'
    TYPE_BUDGET = 'budget'
    TYPE_TIMESERIES = 'timeseries'

    def __init__(self,plots, panel_type ='profile', title : str = "Unnamed panel",
                 dependant_title = "dependant variable"):

        # self.ncdf_files = ncdf_files
        # self.averaging_start_time = averaging_start_time
        # self.averaging_end_time = averaging_end_time
        # self.independent_min_value = independent_min_value
        # self.independent_max_value = independent_max_value
        self.panel_type = panel_type
        self.all_plots = plots

        # self.z = NetCdfVariable('altitude', ncdf_files, avging_start_time=averaging_start_time, avging_end_time=averaging_end_time)
        # self.z_min_idx, self.z_max_idx = self.__getStartEndIndex__(self.z.data, self.independent_min_value, self.independent_max_value)
        # self.z.data = self.z.data[self.z_min_idx:self.z_max_idx]

        # self.time = NetCdfVariable('time', ncdf_files)
        # self.time.data = self.time.data[self.independent_min_value:self.independent_max_value]

        # self.base_plot = self.__var_to_lineplot__(variable_name, self.ncdf_files, label="current clubb")
        # self.all_plots.append(self.base_plot)

        # if sam_file is not None:
        #     self.sam_plot = self.__var_to_lineplot__(CLUBB_TO_SAM[variable_name], [sam_file], label="LES output", line_format="k-",avg_axis=1)
        #     self.all_plots.append(self.sam_plot)

        self.title = title
        self.dependant_title = dependant_title
        self.x_title = "x title unassigned"
        self.y_title = "y title unassigned"
        self.__init_axis_titles__()
        # self.blacklisted = blacklisted



    def __var_to_lineplot__(self, varname, ncdf_files, label ="", line_format ="", avg_axis=0): # TODO is avg_axis appropriate here?
        '''

        :param varname:
        :return:
        '''

        lineplot = None
        if self.panel_type is Panel.TYPE_PROFILE:
            for file in ncdf_files:
                if varname in file.variables.keys():
                    variable = NetCdfVariable(varname, [file], avging_start_time=self.averaging_start_time, avging_end_time=self.averaging_end_time, avg_axis=avg_axis)
                    variable.data = variable.data[self.z_min_idx:self.z_max_idx]
                    lineplot = Lineplot(variable, self.z, label=label, line_format=line_format)

        elif self.panel_type is Panel.TYPE_BUDGET:
            pass
        elif self.panel_type is Panel.TYPE_TIMESERIES:
            variable = NetCdfVariable(varname, ncdf_files, avging_start_time=self.averaging_start_time, avging_end_time=self.averaging_end_time, avg_axis=1)
            variable.data = variable.data[self.independent_min_value:self.independent_max_value]
            lineplot = Lineplot(self.time, variable, label=label, line_format=line_format)
        else:
            raise ValueError('Invalid panel type ' + self.panel_type + '. Valid options are profile, budget, timeseries')
        return lineplot

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
        Save the panel's graphical representation to a file in the 'output' folder from launch parameters

        :param plotter:
        :return:
        '''
        plotter = Plotter()
        plotter.plot(self.all_plots, casename, title=self.title, x_title=self.x_title, y_title=self.y_title)
