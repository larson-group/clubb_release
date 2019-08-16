'''
:author: Nicolas Strike
:date: Mid 2019
'''
import sys
from _warnings import warn

from netCDF4._netCDF4 import Dataset

from Panel import Panel
from DataReader import DataReader, NetCdfVariable
from Line import Line


class VariableGroup:
    '''
    This is the parent PanelGroup class. All other panel groups
    should be created as a subclass of this one. Properties and
    methods common to each PanelGroup can be found here.


    A PanelGroup child defines Lines, Panels, and is responsible for
    calculating any 'calculated' variables from netcdf
    '''

    def __init__(self, ncdf_datasets, case, sam_file=None, coamps_file=None):
        print("\tGenerating variable-group data")
        self.variables = []
        self.panels = []
        self.defualt_panel_type = Panel.TYPE_PROFILE
        self.sam_file = sam_file
        self.coamps_file = coamps_file
        self.ncdf_files = ncdf_datasets
        self.case = case
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

    def addClubbVariable(self, variable_def_dict):
        '''
        Given basic details about a variable like name,
        create lines/panels for the variable
        :return:
        '''
        data_reader = DataReader()
        clubb_name = variable_def_dict['clubb_name']
        fill_zeros = False
        lines = None

        sam_name = None
        sam_file = self.sam_file
        sam_conv_factor = 1
        if 'fill_zeros' in variable_def_dict.keys():
            fill_zeros = variable_def_dict['fill_zeros']
        if 'sam_calc' in variable_def_dict.keys():
            sam_file = None  # don't try to autoplot sam if sam is a calculated value
        if 'sam_name' in variable_def_dict.keys():
            sam_name = variable_def_dict['sam_name']
            sam_file = self.sam_file  # redefine sam_file incase sam_calc wiped it
        if 'sam_conv_factor' in variable_def_dict.keys():
            sam_conv_factor = variable_def_dict['sam_conv_factor']

        coamps_file = self.coamps_file
        coamps_name = None
        coamps_conv_factor = 1
        if 'fill_zeros' in variable_def_dict.keys():
            fill_zeros = variable_def_dict['fill_zeros']
        if 'coamps_calc' in variable_def_dict.keys():
            coamps_file = None  # don't try to autoplot coamps if coamps is a calculated value
        if 'coamps_name' in variable_def_dict.keys():
            coamps_name = variable_def_dict['coamps_name']
            coamps_file = self.coamps_file  # redefine coamps_file incase coamps_calc wiped it
        if 'coamps_conv_factor' in variable_def_dict.keys():
            coamps_conv_factor = variable_def_dict['coamps_conv_factor']
            
        if 'lines' in variable_def_dict.keys():
            lines = variable_def_dict['lines']
        panel_type = self.defualt_panel_type
        fallback = None
        if 'fallback_func' in variable_def_dict.keys():
            fallback = variable_def_dict['fallback_func']
        if 'type' in variable_def_dict.keys():
            panel_type = variable_def_dict['type']

        lines = self.__getVarLines__(clubb_name, self.ncdf_files, start_time=self.start_time,
                                     end_time=self.end_time,  label="current clubb",
                                     line_format=self.default_line_format, fill_zeros = fill_zeros,
                                     override_panel_type=panel_type, fallback_func=fallback, lines=lines,
                                     sam_name=sam_name, sam_file=sam_file, sam_conv_factor=sam_conv_factor,
                                     coamps_name=coamps_name, coamps_file=coamps_file, coamps_conv_factor=coamps_conv_factor)
        variable_def_dict['plots'] = lines
        if 'title' not in variable_def_dict.keys():
            if panel_type == Panel.TYPE_BUDGET:
                variable_def_dict['title'] = variable_def_dict['clubb_name']
            else:
                imported_title = data_reader.getLongName(self.ncdf_files, clubb_name)
                variable_def_dict['title'] = imported_title
        if 'axis_title' not in variable_def_dict.keys():
            if panel_type == Panel.TYPE_BUDGET:
                variable_def_dict['axis_title'] = "["+data_reader.__getUnits__(self.ncdf_files, clubb_name)+"]"
            else:
                imported_axis_title = data_reader.getAxisTitle(self.ncdf_files, clubb_name)
                variable_def_dict['axis_title'] = imported_axis_title
        if 'sam_calc' in variable_def_dict.keys() and sam_file is not None:
            samplot = variable_def_dict['sam_calc']()
            lines.append(samplot)
        self.variables.append(variable_def_dict)

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
            panel_type = self.defualt_panel_type
            if 'type' in variable.keys():
                panel_type = variable['type']
            panel = Panel(plotset, title=title, dependant_title=axis_label, panel_type=panel_type)
            self.panels.append(panel)

    def __getVarLines__(self, varname, ncdf_datasets, label="", line_format="", avg_axis=0, override_panel_type=None,
                        start_time=0, end_time=-1, fallback_func=None, fill_zeros=False, lines=None, conversion_factor=1,
                        sam_name=None, sam_file=None, sam_conv_factor=1,
                        coamps_name=None, coamps_file=None, coamps_conv_factor=1):
        """
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
        """
        panel_type = self.defualt_panel_type
        if override_panel_type is not None:
            panel_type = override_panel_type
        all_lines = []
        if sam_file is not None and sam_name is not None:
            sam_plot = self.__getVarLines__(sam_name, {'sam': sam_file}, label="SAM-LES",
                                            line_format="k-", avg_axis=avg_axis, conversion_factor=sam_conv_factor,
                                            start_time=start_time, end_time=end_time, override_panel_type=panel_type,
                                            fallback_func=fallback_func, fill_zeros=fill_zeros)
            all_lines.extend(sam_plot)

        if coamps_file is not None and coamps_name is not None:
            coamps_plot = self.__getVarLines__(coamps_name, {'coamps': coamps_file}, label="COAMPS-LES",
                                            line_format="k-", avg_axis=avg_axis, conversion_factor=coamps_conv_factor,
                                            start_time=start_time, end_time=end_time, override_panel_type=panel_type,
                                            fallback_func=fallback_func, fill_zeros=fill_zeros)
            all_lines.extend(coamps_plot)

        if isinstance(ncdf_datasets, Dataset):
            ncdf_datasets = {'converted_to_dict': ncdf_datasets}
        line = None
        for dataset in ncdf_datasets.values():
            if varname not in dataset.variables.keys():
                continue  # Skip loop if varname isn't in the dataset
            if panel_type is Panel.TYPE_PROFILE:
                line = self.__get_profile_line__(varname, dataset, label, line_format, conversion_factor, avg_axis, fill_zeros)
                break
            elif panel_type is Panel.TYPE_BUDGET:
                # for overline in lines:
                    if varname in dataset.variables.keys():
                        budget_lines = self.__get_budget_lines__(lines, dataset, fill_zeros)
                        all_lines.extend(budget_lines)
            elif panel_type is Panel.TYPE_TIMESERIES:
                line = self.__get_timeseries_line__(varname, dataset, end_time, conversion_factor, label, line_format)
            else:
                raise ValueError('Invalid panel type ' + panel_type + '. Valid options are profile, budget, timeseries')
        if line is None and panel_type != Panel.TYPE_BUDGET:
            print("\tFailed to find variable " + varname + " in case " + self.casename +
                 ". Attempting to use fallback function.")
            line = self.__getVarFromFallback__(fallback_func, varname)
        if panel_type != Panel.TYPE_BUDGET:
            all_lines.append(line)
        return all_lines

    def __get_profile_line__(self, varname, dataset, label, line_format, conversion_factor, avg_axis, fill_zeros):
        """
        Assumes variable can be plotted as a profile and returns a Line object
        representing the given variable for a profile plot. 
        
        :param variable: 
        :return: Line object representing the given variable for a profile plot
        """
        variable = NetCdfVariable(varname, dataset, start_time=self.start_time, end_time=self.end_time, avg_axis=avg_axis,
                                  conversion_factor=conversion_factor, fill_zeros=fill_zeros)
        if 'altitude' in dataset.variables.keys():
            z = NetCdfVariable('altitude', dataset, start_time=self.start_time, end_time=self.end_time)
        elif 'z' in dataset.variables.keys():
            z = NetCdfVariable('z', dataset, start_time=self.start_time, end_time=self.end_time)
        else:
            z = NetCdfVariable('lev', dataset, start_time=self.start_time, end_time=self.end_time)

        variable.constrain(self.height_min_value, self.height_max_value, data=z.data)
        z.constrain(self.height_min_value, self.height_max_value)
        line = Line(variable, z, label=label, line_format=line_format)
        return line

    def __get_timeseries_line__(self, varname, dataset, end_time, conversion_factor, label, line_format):
        '''

        :param variable:
        :param dataset:
        :param label:
        :param line_format:
        :return:
        '''
        # TODO redeclaring variable is a temp fix until timeseries is auto-discovered
        variable = NetCdfVariable(varname, dataset, start_time=0, end_time=end_time, avg_axis=1, conversion_factor=conversion_factor)
        time = NetCdfVariable('time', dataset)
        variable.constrain(0, variable.end_time, data=time.data)
        time.constrain(0, variable.end_time)
        line = Line(time, variable, label=label, line_format=line_format)
        return line

# TODO convert to add_overline()
    def __get_budget_lines__(self, lines, dataset, fill_zeros):
        '''

        :param variable:
        :param dataset:
        :param label:
        :return:
        '''
        output_lines = []
        for line_definition in lines:
            varname = line_definition['clubb_name']
            if varname not in dataset.variables.keys():
                print("\tFailed to find variable " + line_definition['clubb_name'] + " in case " + self.casename +
                     ". Attempting to use fallback function.")
                if 'fallback_func' in line_definition.keys():
                    fallback = line_definition['fallback_func']
                    fallback_output = self.__getVarFromFallback__(fallback, line_definition['clubb_name'])
                    output_lines.append(fallback_output)

                    return output_lines

                elif not fill_zeros:
                    raise TypeError("Failed to find variable " + varname + " in clubb output for case " +
                                    self.casename + " and there is no fallback function specified.\nIf this is expected "
                                                    "(e.g. this model doesn't output the " + varname +
                                    " variable) then please add a fallback function to the variable's definition to manually calculate it or"
                                    " allow allow the variable to be filled with zeros.\n"
                                    "If it is not expected, please  make sure the correct .nc files are being loaded; ncdump is the recommended "
                                    "tool for verifying a variable/data exists in a file.")
                else:
                    print("\t", varname, " fallback function not found. Filling values with zeros instead.")

            variable = NetCdfVariable(line_definition['clubb_name'], dataset, start_time=self.start_time, end_time=self.end_time, fill_zeros=fill_zeros)
            if 'altitude' in dataset.variables.keys():
                z = NetCdfVariable('altitude', dataset, start_time=self.start_time, end_time=self.end_time)
            else:
                z = NetCdfVariable('z', dataset, start_time=self.start_time, end_time=self.end_time)
            variable.constrain(self.height_min_value, self.height_max_value, data=z.data)
            z.constrain(self.height_min_value, self.height_max_value)
            line_definition = Line(variable, z, label=line_definition['label'], line_format="")  # uses auto-generating line format
            output_lines.append(line_definition)
        return output_lines

    def __getVarFromFallback__(self, fallback, varname):
        '''
        Some older model output doesn't contain all variables, e.g cgils_s12 doesn't contain WPTHLP.
        This method attempts to use a function under the name 'fallback_func' to generate the requested
        data similarly to how 'sam_calc' calculates sam output into a clubb format. If successful,
        returns varline data, otherwise raises error
        :return: None
        '''
        varline = fallback()
        print("\tFallback for ", varname, " successful")
        return varline

    def __getFallbackVar__(self, varname, dataset, conversion_factor = 1, fill_zeros = False):
        '''
        This function is used within a fallback function to get the data of a certain variable,
        constrained between a min/max height.

        :param varname:
        :return:
        '''
        if 'z' in dataset.variables.keys():
            z_ncdf = NetCdfVariable('z', dataset, 1)
        else:
            z_ncdf = NetCdfVariable('altitude', dataset, 1)
        var_ncdf = NetCdfVariable(varname, dataset, conversion_factor, start_time=self.start_time, end_time=self.end_time, fill_zeros=fill_zeros)
        var_ncdf.constrain(self.height_min_value, self.height_max_value, data=z_ncdf.data)
        var_data = var_ncdf.data

        return var_data