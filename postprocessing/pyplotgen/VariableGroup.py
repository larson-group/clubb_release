'''
:author: Nicolas Strike
:date: Mid 2019
'''
from netCDF4._netCDF4 import Dataset

from pyplotgen.Panel import Panel
from pyplotgen.DataReader import DataReader, NetCdfVariable
from pyplotgen.Lineplot import Lineplot


class VariableGroup:
    '''
    This is the parent PanelGroup class. All other panel groups
    should be created as a subclass of this one. Properties and
    methods common to each PanelGroup can be found here.


    A PanelGroup child defines Lineplots, Panels, and is responsible for
    calculating any 'calculated' variables from netcdf
    '''

    def __init__(self, ncdf_datasets, case, sam_file=None):
        self.panels = []
        self.panel_type = Panel.TYPE_PROFILE
        self.sam_file = sam_file
        self.ncdf_files = ncdf_datasets
        self.averaging_start_time = case.averaging_start_time
        self.averaging_end_time = case.averaging_end_time
        self.timeseries_start_time = case.timeseries_start_time
        self.timeseries_end_time = case.timeseries_end_time
        self.height_min_value = case.height_min_value
        self.height_max_value = case.height_max_value


        ### Initialize Height ###
        self.z = NetCdfVariable('altitude', ncdf_datasets['zm'], avging_start_time=self.averaging_start_time,
                                avging_end_time=self.averaging_end_time)
        self.z_min_idx, self.z_max_idx = self.__getStartEndIndex__(self.z.data, self.height_min_value,
                                                                   self.height_max_value)
        self.z.data = self.z.data[self.z_min_idx:self.z_max_idx]

        ### Initialize Time ###
        sec_per_min = 60
        sec_to_min = 1 / sec_per_min
        self.time = NetCdfVariable('time', ncdf_datasets, conversion_factor=sec_to_min)

        ### Initialize Sam Height ###
        if sam_file != None:
            self.z_sam = NetCdfVariable('z', sam_file, avging_start_time=self.averaging_start_time,
                                        avging_end_time=self.averaging_end_time)
            self.z_sam_min_idx, self.z_sam_max_idx = self.__getStartEndIndex__(self.z_sam.data, self.height_min_value,
                                                                               self.height_max_value)
            self.z_sam.data = self.z_sam.data[self.z_sam_min_idx:self.z_sam_max_idx]



        for variable in self.variables:
            self.addClubbVariable(variable)
        self.generatePanels()

    def get_var_from_ncdf(self, netcdf_variable):
        '''
        Retrieve numerical data from netcdf

        :param netcdf_variable: Of type NetCDFVariable
        :return:
        '''
        data_reader = DataReader()
        return data_reader.getVarData(netcdf_variable.ncdf_data, netcdf_variable)

    def addClubbVariable(self, variable):
        '''
        Given basic details about a variable like name,
        create plots/panels for the variable
        :return:
        '''
        data_reader = DataReader()
        clubb_name = variable['clubb_name']
        sam_name = None
        sam_file = self.sam_file
        sam_conv_factor = 1
        if 'sam_calc' in variable.keys():
            sam_file = None # don't try to autoplot sam if sam is a calculated value
        if 'sam_name' in variable.keys():
            sam_name = variable['sam_name']
            sam_file = self.sam_file # redefine sam_file incase sam_calc wiped it
        if 'sam_conv_factor' in variable.keys():
            sam_conv_factor = variable['sam_conv_factor']
        plots = self.getVarLinePlots(clubb_name, self.ncdf_files, averaging_start_time=self.averaging_start_time,
                                     averaging_end_time=self.averaging_end_time, sam_name= sam_name, sam_file=sam_file,
                                     sam_conv_factor=sam_conv_factor, label="current clubb", line_format='r--')
        variable['plots'] = plots
        if 'title' not in variable.keys():
            imported_title = data_reader.getLongName(self.ncdf_files, clubb_name)
            variable['title'] = imported_title
        if 'axis_title' not in variable.keys():
            imported_axis_title = data_reader.getAxisTitle(self.ncdf_files, clubb_name)
            variable['axis_title'] = imported_axis_title
        if 'sam_calc' in variable.keys():
            samplot = variable['sam_calc']()
            plots.append(samplot)

    def generatePanels(self):
        '''
        Generates a set of panels from the plots stored in self.
        Does not return anything, simply assigns the panels into
        self.panels
        :return:
        '''
        for variable in self.variables:
            title = variable['title']
            axis_label = variable['axis_title']
            #        liq_pot_temp = Panel(thlm_plots, title="Liquid Water Potential Temperature, Theta l", dependant_title="thlm [K]")
            plotset = variable['plots']
            panel = Panel(plotset, title=title, dependant_title=axis_label)
            self.panels.append(panel)

    def getVarLinePlots(self, varname, ncdf_datasets, label="", line_format="", avg_axis=0, override_panel_type=None,
                        averaging_start_time=0, averaging_end_time=-1, sam_name=None, sam_file=None, conversion_factor=1,
                        sam_conv_factor=1):
        '''
        Get a list of Lineplot objects for a specific clubb variable. If sam_file is specified it will also
        attempt to generate Lineplots for the SAM equivalent variables, using the name conversions found in
        VarnameConversions.py. If a SAM variable needs to be calculated (uses an equation) then it will have
        to be created within that variable group's file and not here.

        :param varname: str name of the clubb variable to be plotted, case sensitive
        :param ncdf_datasets: List of Dataset objects containing clubb or sam netcdf data
        :param label: Label to give the base-plotAll on the legend. This is normally 'current clubb', but not provided as default to help avoid debugging confusion.
        :param line_format: Line formatting string used by matplotlib's PyPlot
        :param avg_axis: Axis over which to average values. 0 - time average, 1 - height average
        :param override_panel_type: Override the VariableGroup's default panel type
        :param averaging_start_time: Beginning period of the averaging interval. Give a time VALUE, e.g. 240
        :param averaging_end_time: Ending period of the averaging interval. Give a time VALUE, e.g. 120
        :param sam_file: Dataset object containing SAM plotAll data
        :param conversion_factor: A multiplying factor used to scale clubb output. Defaults to 1.
        :param sam_conv_factor: A multiplying factor used to scale sam output. Defaults to 1.
        :return: A list of Lineplot objects containing clubb and (if requested) sam data. Returns None if requested variable is not found.
        '''

        if override_panel_type is not None:
            panel_type = override_panel_type
        else:
            panel_type = self.panel_type

        all_plots = []

        if sam_file is not None and sam_name is not None:
            sam_plot = self.getVarLinePlots(sam_name, {'sam': sam_file}, label="LES output",
                                            line_format="k-", avg_axis=1, conversion_factor=sam_conv_factor)
            all_plots.extend(sam_plot)

        if isinstance(ncdf_datasets, Dataset):
            ncdf_datasets = {'auto': ncdf_datasets}

        lineplot = None
        if panel_type is Panel.TYPE_PROFILE:

            for file in ncdf_datasets.values():
                if varname in file.variables.keys() or len(ncdf_datasets.values()) == 1:
                    variable = NetCdfVariable(varname, file, avging_start_time=averaging_start_time,
                                              avging_end_time=averaging_end_time, avg_axis=avg_axis,
                                              conversion_factor=conversion_factor)
                    independent_var_data = self.z
                    model_src_reader = DataReader()
                    model_src = model_src_reader.getNcdfSourceModel(file)
                    if model_src == "sam":
                        independent_var_data = self.z_sam
                    min_idx, max_idx = self.__getStartEndIndex__(independent_var_data.data, self.height_min_value,
                                                                 self.height_max_value)
                    variable.data = variable.data[min_idx:max_idx + 1]
                    lineplot = Lineplot(variable, independent_var_data, label=label, line_format=line_format)
                    break

        elif panel_type is Panel.TYPE_BUDGET:
            pass
        elif panel_type is Panel.TYPE_TIMESERIES:
            variable = NetCdfVariable(varname, ncdf_datasets, avging_start_time=averaging_start_time,
                                      avging_end_time=averaging_end_time, avg_axis=1,
                                      conversion_factor=conversion_factor)
            variable.data = variable.data[self.timeseries_start_time:self.timeseries_end_time]
            lineplot = Lineplot(self.time, variable, label=label, line_format=line_format)
        else:
            raise ValueError('Invalid panel type ' + panel_type + '. Valid options are profile, budget, timeseries')
        if lineplot == None:
            raise ValueError('Failed to find variable ' + varname + " in " + str(file))
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
        end_idx = len(data) - 1
        for i in range(0, len(data)):
            # Check for start index
            test_value = data[i]
            if test_value <= start_value and test_value > data[start_idx]:
                start_idx = i
            # Check for end index
            if test_value >= end_value and test_value < data[end_idx]:
                end_idx = i

        return start_idx, end_idx
