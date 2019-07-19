'''
:author: Nicolas Strike
:date: Mid 2019
'''
import sys
from _warnings import warn

from netCDF4._netCDF4 import Dataset

from pyplotgen.Panel import Panel
from pyplotgen.DataReader import DataReader, NetCdfVariable
from pyplotgen.Line import Line


class VariableGroup:
    '''
    This is the parent PanelGroup class. All other panel groups
    should be created as a subclass of this one. Properties and
    methods common to each PanelGroup can be found here.


    A PanelGroup child defines Lines, Panels, and is responsible for
    calculating any 'calculated' variables from netcdf
    '''

    def __init__(self, ncdf_datasets, case, sam_file=None):
        print("\tGenerating variable-group data")
        self.variables = []
        self.panels = []
        self.panel_type = Panel.TYPE_PROFILE
        self.sam_file = sam_file
        self.ncdf_files = ncdf_datasets
        self.casename = case.name
        self.start_time = case.start_time
        self.end_time = case.end_time
        self.height_min_value = case.height_min_value
        self.height_max_value = case.height_max_value
        self.default_line_format = 'b-'

        for variable in self.variable_definitions:
            print("\tProcessing ", variable['clubb_name'])
            if variable['clubb_name'] not in case.blacklisted_variables:
                # Skip this variable if it's blacklisted for the case
                self.addClubbVariable(variable)
        self.generatePanels()

    def addClubbVariable(self, variable):
        '''
        Given basic details about a variable like name,
        create plots/panels for the variable
        :return:
        '''
        data_reader = DataReader()
        clubb_name = variable['clubb_name']
        fill_zeros = False
        sam_name = None
        sam_file = self.sam_file
        sam_conv_factor = 1
        if 'fill_zeros' in variable.keys():
            fill_zeros = variable['fill_zeros']
        if 'sam_calc' in variable.keys():
            sam_file = None  # don't try to autoplot sam if sam is a calculated value
        if 'sam_name' in variable.keys():
            sam_name = variable['sam_name']
            sam_file = self.sam_file  # redefine sam_file incase sam_calc wiped it
        if 'sam_conv_factor' in variable.keys():
            sam_conv_factor = variable['sam_conv_factor']

        panel_type = self.panel_type
        fallback = None
        if 'fallback_func' in variable.keys():
            fallback = variable['fallback_func']
        if 'type' in variable.keys():
            panel_type = variable['type']
        plots = self.__getVarLines__(clubb_name, self.ncdf_files, start_time=self.start_time,
                                     end_time=self.end_time, sam_name=sam_name, sam_file=sam_file,
                                     sam_conv_factor=sam_conv_factor, label="current clubb",
                                     line_format=self.default_line_format, fill_zeros = fill_zeros,
                                     override_panel_type=panel_type, fallback_func=fallback)
        variable['plots'] = plots
        if 'title' not in variable.keys():
            imported_title = data_reader.getLongName(self.ncdf_files, clubb_name)
            variable['title'] = imported_title
        if 'axis_title' not in variable.keys():
            imported_axis_title = data_reader.getAxisTitle(self.ncdf_files, clubb_name)
            variable['axis_title'] = imported_axis_title
        if 'sam_calc' in variable.keys() and sam_file is not None:
            samplot = variable['sam_calc']()
            plots.append(samplot)
        self.variables.append(variable)

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
            plotset = variable['plots']
            panel_type = self.panel_type
            if 'type' in variable.keys():
                panel_type = variable['type']
            panel = Panel(plotset, title=title, dependant_title=axis_label, panel_type=panel_type)
            self.panels.append(panel)

    def __getVarLines__(self, varname, ncdf_datasets, label="", line_format="", avg_axis=0, override_panel_type=None,
                        start_time=0, end_time=-1, sam_name=None, sam_file=None, conversion_factor=1,
                        sam_conv_factor=1, fallback_func=None, fill_zeros=False):
        '''
        Get a list of Line objects for a specific clubb variable. If sam_file is specified it will also
        attempt to generate Lines for the SAM equivalent variables, using the name conversions found in
        VarnameConversions.py. If a SAM variable needs to be calculated (uses an equation) then it will have
        to be created within that variable group's dataset and not here.

        :param varname: str name of the clubb variable to be plotted, case sensitive
        :param ncdf_datasets: List of Dataset objects containing clubb or sam netcdf data
        :param label: Label to give the base-plotAll on the legend. This is normally 'current clubb', but not provided as default to help avoid debugging confusion.
        :param line_format: Line formatting string used by matplotlib's PyPlot
        :param avg_axis: Axis over which to average values. 0 - time average, 1 - height average
        :param override_panel_type: Override the VariableGroup's default panel type
        :param start_time: Beginning period of the averaging interval. Give a time VALUE, e.g. 240
        :param end_time: Ending period of the averaging interval. Give a time VALUE, e.g. 120
        :param sam_file: Dataset object containing SAM plotAll data
        :param conversion_factor: A multiplying factor used to scale clubb output. Defaults to 1.
        :param sam_conv_factor: A multiplying factor used to scale sam output. Defaults to 1.
        :return: A list of Line objects containing clubb and (if requested) sam data. Returns None if requested variable is not found.
        '''

        if override_panel_type is not None:
            panel_type = override_panel_type
        else:
            panel_type = self.panel_type

        all_lines = []

        if sam_file is not None and sam_name is not None:
            # clubb_sec_to_sam_min = 1 / 60
            sam_plot = self.__getVarLines__(sam_name, {'sam': sam_file}, label="LES output",
                                            line_format="k-", avg_axis=avg_axis, conversion_factor=sam_conv_factor,
                                            start_time=start_time,  #* clubb_sec_to_sam_min,
                                            end_time=end_time,  # * clubb_sec_to_sam_min,
                                            override_panel_type=panel_type, fallback_func=fallback_func, fill_zeros=fill_zeros)
            all_lines.extend(sam_plot)

        if isinstance(ncdf_datasets, Dataset):
            ncdf_datasets = {'converted_to_dict': ncdf_datasets}

        line = None
        for dataset in ncdf_datasets.values():

            if varname not in dataset.variables.keys():
                continue  # Skip loop if varname isn't in the dataset

            variable = NetCdfVariable(varname, ncdf_datasets, start_time=start_time,
                                      end_time=end_time, avg_axis=avg_axis,
                                      conversion_factor=conversion_factor, fill_zeros=fill_zeros)
            if panel_type is Panel.TYPE_PROFILE:
                line = self.__get_profile_line__(variable, dataset, label, line_format)
                break
            elif panel_type is Panel.TYPE_BUDGET:
                pass
            elif panel_type is Panel.TYPE_TIMESERIES:
                # TODO redeclaring variable is a temp fix until timeseries is auto-discovered
                variable = NetCdfVariable(varname, ncdf_datasets, start_time=0,
                                          end_time=end_time, avg_axis=1,
                                          conversion_factor=conversion_factor)
                time = NetCdfVariable('time', dataset)
                variable.constrain(0, variable.end_time, data=time.data)
                time.constrain(0, variable.end_time)
                line = Line(time, variable, label=label, line_format=line_format)
            else:
                raise ValueError('Invalid panel type ' + panel_type + '. Valid options are profile, budget, timeseries')
        if line is None:
            warn("\tFailed to find variable " + varname + " in case " + self.casename +
                 ". Attempting to use fallback function.")
            line = self.__getVarFromFallback__(fallback_func, varname)
        all_lines.append(line)
        return all_lines

    def __get_profile_line__(self, variable, dataset, label, line_format):
        '''
        Assumes variable can be plotted as a profile and returns a Line object
        representing the given variable for a profile plot. 
        
        :param variable: 
        :return: Line object representing the given variable for a profile plot
        '''
        if 'altitude' in dataset.variables.keys():
            z = NetCdfVariable('altitude', dataset, start_time=self.start_time, end_time=self.end_time)
        else:
            z = NetCdfVariable('z', dataset, start_time=self.start_time, end_time=self.end_time)
        variable.constrain(self.height_min_value, self.height_max_value, data=z.data)
        z.constrain(self.height_min_value, self.height_max_value)
        line = Line(variable, z, label=label, line_format=line_format)
        return line

    def __getVarFromFallback__(self, fallback, varname):
        '''
        Some older model output doesn't contain all variables, e.g cgils_s12 doesn't contain WPTHLP.
        This method attempts to use a function under the name 'fallback_func' to generate the requested
        data similarly to how 'sam_calc' calculates sam output into a clubb format. If successful,
        returns varline data, otherwise raises error
        :return: None
        '''
        if fallback is None:
            raise TypeError("Failed to find variable " + varname + " in clubb output for case " +
                            self.casename + " and there is no fallback function specified. If this is expected "
                                            "(e.g. this model doesn't output the " + varname +
                            " variable) then please add a fallback function to the variable's definition or "
                            "   allow allow the variable to be filled with zeros. "
                            "If it is not expected, please  make sure the correct .nc files are being loaded.")
        varline = fallback()
        print("\tFallback for ", varname, " successful")
        return varline
