'''
:author: Nicolas Strike
:date: Mid 2019
'''
from pyplotgen import Panel
from pyplotgen.DataReader import DataReader, NetCdfVariable
from pyplotgen.Lineplot import Lineplot
from pyplotgen.VarnameConversions import CLUBB_TO_SAM


class VariableGroup:
    '''
    This is the parent PanelGroup class. All other panel groups
    should be created as a subclass of this one. Properties and
    methods common to each PanelGroup can be found here.


    A PanelGroup child defines Lineplots, Panels, and is responsible for
    calculating any 'calculated' variables from netcdf
    '''
    def __init__(self, ncdf_file):
        self.ncdf_file = ncdf_file
        self.panels = []
        self.panel_type = Panel.Panel.TYPE_PROFILE

    def plot(self):
        '''
        Plots every panel in this group to the output folder specified
        in the pyplotgen launch parameters, unless the panel is blacklisted.
        
        :return: n/a
        '''
        for panel in self.panels:
            if not panel.blacklisted:
                panel.plot(self.ncdf_file)

    def get_var_from_ncdf(self, netcdf_variable:NetCdfVariable):
        '''
        Retrieve numerical data from netcdf

        :param netcdf_variable:
        :return:
        '''
        data_reader = DataReader()
        return data_reader.getVarData(netcdf_variable.ncdf_files, netcdf_variable)


        # super().getLinePlots('thlm', ncdf_files, averaging_start_time=averaging_start_time,
        #                      averaging_end_time=averaging_end_time, independent_min_value=height_min_value,
        #                      independent_max_value=height_max_value, sam_file=sam_file)
    def getLinePlots(self, varname, ncdf_files, label="", line_format="", avg_axis=0, override_panel_type=None, averaging_start_time = 0,
                     averaging_end_time=-1, sam_file=None, conversion_factor=1, sam_conv_factor=1): # TODO is avg_axis appropriate here?
        '''

        :param varname:
        :return:
        '''
        if override_panel_type is not None:
            panel_type = override_panel_type
        else:
            panel_type = self.panel_type

        all_plots = []

        if sam_file is not None:
            sam_plot = self.getLinePlots(CLUBB_TO_SAM[varname], [sam_file], label="LES output", line_format="k-",avg_axis=1, conversion_factor=sam_conv_factor)
            all_plots.extend(sam_plot)

        lineplot = None
        if panel_type is Panel.Panel.TYPE_PROFILE:
            for file in ncdf_files:
                if varname in file.variables.keys():
                    variable = NetCdfVariable(varname, [file], avging_start_time=averaging_start_time, avging_end_time=averaging_end_time, avg_axis=avg_axis, conversion_factor=conversion_factor)
                    independent_var_data = self.z
                    model_src_reader = DataReader()
                    model_src = model_src_reader.getNcdfSourceModel(file)
                    if model_src == "sam":
                        independent_var_data = self.z_sam
                    min_idx,max_idx = self.__getStartEndIndex__(independent_var_data.data, self.height_min_value,self.height_max_value)
                    variable.data = variable.data[min_idx:max_idx+1]
                    lineplot = Lineplot(variable, independent_var_data, label=label, line_format=line_format)
                    break

        elif panel_type is Panel.Panel.TYPE_BUDGET:
            pass
        elif panel_type is Panel.Panel.TYPE_TIMESERIES:
            variable = NetCdfVariable(varname, ncdf_files, avging_start_time=averaging_start_time, avging_end_time=averaging_end_time, avg_axis=1, conversion_factor=conversion_factor)
            variable.data = variable.data[self.timeseries_start_time:self.timeseries_end_time]
            lineplot = Lineplot(self.time, variable, label=label, line_format=line_format)
        else:
            raise ValueError('Invalid panel type ' + panel_type + '. Valid options are profile, budget, timeseries')

        all_plots.append(lineplot)
        return all_plots

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

