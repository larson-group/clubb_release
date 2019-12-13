"""
:author: Nicolas Strike
:date: Mid 2019
"""
import os
from warnings import warn

import numpy as np
from netCDF4._netCDF4 import Dataset

from config import Case_definitions
from config import Style_definitions
from src.DataReader import DataReader, NetCdfVariable
from src.Line import Line
from src.Panel import Panel


class VariableGroup:
    """
    This is the parent VariableGroup class. All other panel groups
    should be created as a subclass of this one. Properties and
    methods common to each PanelGroup can be found here.

    A VariableGroup child defines Lines, Panels, etc., and is responsible for
    calculating any 'calculated' variables from netcdf
    """

    def __init__(self, ncdf_datasets, case, sam_file=None, coamps_file=None, r408_dataset=None, hoc_dataset=None,
                 e3sm_dataset=None):
        """
        Initialize common VariableGroup parameters
        :param ncdf_datasets: A dictionary or Netcdf dataset containing the data to be plotted
        :param case: An instance of a Case object
        :param sam_file: NetCDF4 Dataset object containing sam ouput
        :param coamps_file: NetCDF4 Dataset object containing coamps ouput
        :param r408_dataset: NetCDF4 Dataset object containing clubb r408 ('chris golaz') output
        """
        print("\tGenerating variable-group data")
        self.variables = []
        self.panels = []
        self.defualt_panel_type = Panel.TYPE_PROFILE
        self.ncdf_files = ncdf_datasets
        self.case = case
        self.sam_file = sam_file
        self.e3sm_dataset = e3sm_dataset
        self.coamps_file = coamps_file
        self.r408_dataset = r408_dataset
        self.hoc_dataset = hoc_dataset
        self.casename = case.name
        self.start_time = case.start_time
        self.end_time = case.end_time
        self.height_min_value = case.height_min_value
        self.height_max_value = case.height_max_value

        for variable in self.variable_definitions:
            print("\tProcessing ", variable['aliases'])

            # Only add variable if none of the aliases are blacklisted
            if len(list(set(variable['aliases']).intersection(
                    case.blacklisted_variables))) == 0:  # variable['aliases'] not in case.blacklisted_variables:
                # Skip this variable if it's blacklisted for the case
                self.addVariable(variable)
        self.generatePanels()

    def addVariable(self, variable_def_dict):
        """
        Given basic details about a variable, this
        creates lines/panels for the variable and then adds that
        variable to the list of all variables for this VariableGroup

        :param variable_def_dict: A dict containing the information defining a variable. These definitions are declared
        inside a VariableGroup child class, e.g. VariableGroupBase.

        Valid dict keys/options include:

        *aliases*: A list of names various models refer to this variable as. E.g. ['wprtp', 'WPRTP', 'wpqtp']

        *sam_calc*: (optional) A functional reference to a method that calculates a sam variable. This is given as the name of the
        function *without*  a set of parenthesis () after the name. E.g. self.getThlmSamCalc

        *coamps_calc*: (optional) A functional reference to a method that calculates a coamps variable. This is given as the name of the
        function *without* a set of parenthesis () after the name. E.g. self.getThlmCoampsCalc

        *r408_calc*: (optional) A functional reference to a method that calculates a r408 variable. This is given as the name of the
        function *without*  a set of parenthesis () after the name. E.g. self.getR408SamCalc

        *hoc_calc*: (optional) A functional reference to a method that calculates a hoc-2005 variable. This is given as the name of the
        function *without* a set of parenthesis () after the name. E.g. self.getThlmHocCalc

        *e3sm_calc*: (optional)

        *sam_conv_factor*: (optional) Numeric value to scale a sam variable by. E.g. 1/1000, or 100

        *coamps_conv_factor*: (optional) Numeric value to scale a coamps variable by. E.g. 1/1000, or 100

        *r408_conv_factor*: (optional) Numeric value to scale a clubb r408 variable by. E.g. 1/1000, or 100

        *e3sm_conv_factor*: (optional) Numeric value to scale an e3sm variable by. E.g. 1/1000, or 100

        *type*: (optional) Override the default type 'profile' with either 'budget' or 'timeseries'

        *fallback_func*: (optional) If a variable is not found within a dataset and a fallback_func is specified, this
        function will be called to attempt retrieving the variable (before filling zeros if fill_zeros=True). Like the
        model_calc option, this is a functional reference to a method that calculates the given variable. E.g. self.getWpthlpFallback

        *title*: (optional) Override the default panel title, or provide one if it's not specified in the netcdf file.

        *axis_title*: (optional) Override the default dependant axis title, or provide one if it's not specified in the netcdf file.

        *fill_zeros*: (optional) If a variable isn't found in netcdf output, and there is no fallback_func (or the fallback failed),
        setting this to True allows PyPlotgen to fill all datapoints with 0 and continue plotting.

        *lines*: * (budget variables only) Defines lines to plot for budget cases. Passed seperately because it's a lot of text.
        This is given in the form of a list of lines, here's an example:

        .. code-block:: python
            :linenos:

            thlm_lines = [
            {'aliases': ['thlm_bt'], 'legend_label': 'thlm_bt'},
            {'aliases': ['thlm_ma'], 'legend_label': 'thlm_ma'},
            {'aliases': ['thlm_ta'], 'legend_label': 'thlm_ta'},
            {'aliases': ['thlm_mc'], 'legend_label': 'thlm_mc'},
            {'aliases': ['thlm_clipping'], 'legend_label': 'thlm_bt', 'fallback_func': self.getThlmClipping},
            {'aliases': ['radht'], 'legend_label': 'radht'},
            {'aliases': ['lsforcing'], 'legend_label': 'lsforcing', 'fallback_func': self.getLsforcing},
            {'aliases': ['thlm_residual'], 'legend_label': 'thlm_residual', 'fallback_func': self.getThlmResidual},
            ]

        Here are a couple examples of other (non-budget) variable definitions:

        .. code-block:: python
            :linenos:

            {
            'aliases': ['Skrt_zt'],
            'sam_calc': self.getSkrtZtLesCalc,
            'coamps_calc': self.getSkrtZtLesCalc,
            'fill_zeros': True
            }

            {
            'aliases': ['lwp', 'CWP'],
            'type': Panel.TYPE_TIMESERIES,
            'sam_conv_factor': 1/1000
            },

        :return: None

        """

        data_reader = DataReader()
        aliases = variable_def_dict['aliases']
        fill_zeros = False
        lines = None

        if 'fill_zeros' in variable_def_dict.keys():
            fill_zeros = variable_def_dict['fill_zeros']

        # TODO refactor these chunks into method arguments
        sam_file = self.sam_file
        sam_conv_factor = 1
        if 'sam_calc' in variable_def_dict.keys():
            sam_file = None  # don't try to autoplot sam if sam is a calculated value
        if 'sam_conv_factor' in variable_def_dict.keys():
            sam_conv_factor = variable_def_dict['sam_conv_factor']

        coamps_file = self.coamps_file
        coamps_conv_factor = 1
        if 'coamps_calc' in variable_def_dict.keys():
            coamps_file = None  # don't try to autoplot coamps if coamps is a calculated value
        if 'coamps_conv_factor' in variable_def_dict.keys():
            coamps_conv_factor = variable_def_dict['coamps_conv_factor']

        r408_dataset = self.r408_dataset
        r408_conv_factor = 1
        if 'r408_calc' in variable_def_dict.keys():
            r408_dataset = None  # don't try to autoplot if calculated value
        if 'r408_conv_factor' in variable_def_dict.keys():
            r408_conv_factor = variable_def_dict['r408_conv_factor']

        hoc_dataset = self.hoc_dataset
        hoc_conv_factor = 1
        if 'hoc_calc' in variable_def_dict.keys():
            hoc_dataset = None  # don't try to autoplot if calculated value
        if 'hoc_conv_factor' in variable_def_dict.keys():
            hoc_conv_factor = variable_def_dict['hoc_conv_factor']

        e3sm_dataset = self.e3sm_dataset
        e3sm_conv_factor = 1
        if 'e3sm_calc' in variable_def_dict.keys():
            e3sm_dataset = None  # don't try to autoplot e3sm if e3sm is a calculated value
        if 'e3sm_conv_factor' in variable_def_dict.keys():
            e3sm_conv_factor = variable_def_dict['e3sm_conv_factor']

        if 'lines' in variable_def_dict.keys():
            lines = variable_def_dict['lines']
        panel_type = self.defualt_panel_type
        fallback = None
        if 'fallback_func' in variable_def_dict.keys():
            fallback = variable_def_dict['fallback_func']
        if 'type' in variable_def_dict.keys():
            panel_type = variable_def_dict['type']

        all_lines = []
        if sam_file is not None:
            all_lines.extend(self.__getVarLines__(aliases, sam_file, conversion_factor=sam_conv_factor,
                                                  label=Style_definitions.SAM_LABEL,
                                                  line_format=Style_definitions.LES_LINE_STYLE, fill_zeros=fill_zeros,
                                                  override_panel_type=panel_type,
                                                  fallback_func=fallback, lines=lines))

        if coamps_file is not None:
            all_lines.extend(self.__getVarLines__(aliases, coamps_file, conversion_factor=coamps_conv_factor,
                                                  label=Style_definitions.COAMPS_LABEL,
                                                  line_format=Style_definitions.LES_LINE_STYLE, fill_zeros=fill_zeros,
                                                  override_panel_type=panel_type,
                                                  fallback_func=fallback, lines=lines))

        if r408_dataset is not None:
            all_lines.extend(self.__getVarLines__(aliases, r408_dataset, conversion_factor=r408_conv_factor,
                                                  label=Style_definitions.GOLAZ_LABEL,
                                                  line_format=Style_definitions.GOLAZ_BEST_R408_LINE_STYLE,
                                                  fill_zeros=fill_zeros,
                                                  override_panel_type=panel_type, fallback_func=fallback, lines=lines))

        if hoc_dataset is not None:
            all_lines.extend(self.__getVarLines__(aliases, hoc_dataset, conversion_factor=hoc_conv_factor,
                                                  label=Style_definitions.HOC_LABEL,
                                                  line_format=Style_definitions.HOC_LINE_STYLE, fill_zeros=fill_zeros,
                                                  override_panel_type=panel_type, fallback_func=fallback, lines=lines))

        if e3sm_dataset is not None:
            all_lines.extend(self.__getVarLines__(aliases, e3sm_dataset, conversion_factor=e3sm_conv_factor,
                                                  label=Style_definitions.E3SM_LABEL,
                                                  line_format=Style_definitions.E3SM_LINE_STYLE, fill_zeros=fill_zeros,
                                                  override_panel_type=panel_type, fallback_func=fallback, lines=lines))


        num_folders_plotted = 0
        if self.ncdf_files is not None:
            for ncdf_files_subfolder in self.ncdf_files:
                datasets = self.ncdf_files[ncdf_files_subfolder]
                label = os.path.basename(ncdf_files_subfolder)
                if num_folders_plotted < len(Style_definitions.CLUBB_LABEL_OVERRIDE):
                    label = Style_definitions.CLUBB_LABEL_OVERRIDE[num_folders_plotted]
                all_lines.extend(self.__getVarLines__(aliases, datasets, label=label, fill_zeros=fill_zeros,
                                                      override_panel_type=panel_type, fallback_func=fallback, lines=lines))
                num_folders_plotted += 1
        else:
            label = ""

        aliases = variable_def_dict['aliases']
        clubb_name = aliases[0]
        variable_def_dict['plots'] = all_lines
        first_input_datasets = self.getTextDefiningDataset()
        if 'title' not in variable_def_dict.keys():
            if panel_type == Panel.TYPE_BUDGET:
                variable_def_dict['title'] = label + ' ' + clubb_name
            else:
                imported_title = data_reader.getLongName(first_input_datasets, aliases)
                variable_def_dict['title'] = imported_title
        if 'axis_title' not in variable_def_dict.keys():
            if panel_type == Panel.TYPE_BUDGET:
                any_varname_with_budget_units = all_lines[0].label
                variable_def_dict['axis_title'] = "[" + data_reader.__getUnits__(first_input_datasets, any_varname_with_budget_units) + "]"
            else:
                imported_axis_title = data_reader.getAxisTitle(first_input_datasets, clubb_name)
                variable_def_dict['axis_title'] = imported_axis_title

        if 'sam_calc' in variable_def_dict.keys() and self.sam_file is not None and data_reader.getNcdfSourceModel(
                self.sam_file) == 'sam':
            samplot_data, z = variable_def_dict['sam_calc']()
            samplot = Line(samplot_data, z, line_format=Style_definitions.LES_LINE_STYLE,
                           label=Style_definitions.SAM_LABEL)
            variable_def_dict['plots'].append(samplot)

        if 'coamps_calc' in variable_def_dict.keys() and self.coamps_file is not None and data_reader.getNcdfSourceModel(
                self.coamps_file) == 'coamps':
            coampsplot_data, z = variable_def_dict['coamps_calc']()
            coampsplot = Line(coampsplot_data, z, line_format=Style_definitions.LES_LINE_STYLE,
                              label=Style_definitions.COAMPS_LABEL)
            variable_def_dict['plots'].append(coampsplot)

        if 'e3sm_calc' in variable_def_dict.keys() and self.e3sm_dataset is not None and data_reader.getNcdfSourceModel(
                self.e3sm_dataset) == 'e3sm':
            e3smplot_data, z = variable_def_dict['e3sm_calc']()
            e3smplot = Line(e3smplot_data, z, line_format=Style_definitions.E3SM_LINE_STYLE,
                           label=Style_definitions.E3SM_LABEL)
            variable_def_dict['plots'].append(e3smplot)

        if len(variable_def_dict['plots']) > 0:
            self.variables.append(variable_def_dict)


    def getTextDefiningDataset(self):
        """
        Finds a dataset from which text names and descriptions can be pulled from.
        By default, it is always clubb, however if pyplotgen is not plotting clubb then
        it will return the first dataset it can find.
        :return:
        """

    #     self.coamps_file = coamps_file
    # self.r408_dataset = r408_dataset
    # self.hoc_dataset = hoc_dataset

        datasets = []
        if self.ncdf_files is not None:
            for input_folder in self.ncdf_files.values():
                for dataset in input_folder.values():
                    datasets.append(dataset)
        if self.e3sm_dataset is not None:
            datasets.append(self.e3sm_dataset)
        if self.sam_file is not None:
            datasets.append(self.sam_file)
        if self.coamps_file is not None:
            datasets.extend(self.coamps_file.values())
        if self.r408_dataset is not None:
            datasets.extend(self.r408_dataset.values())
        if self.hoc_dataset is not None:
            datasets.extend(self.hoc_dataset.values())
        if len(datasets) == 0:
            raise FileExistsError("No dataset could be found to pull text descriptions from")
        return datasets

    def generatePanels(self):
        """
        Generates a set of panels from the plots stored in self.
        Does not return anything, simply assigns the panels into
        self.panels
        :return: None
        """
        for variable in self.variables:
            title = variable['title']
            axis_label = variable['axis_title']
            plotset = variable['plots']
            panel_type = self.defualt_panel_type
            if 'type' in variable.keys():
                panel_type = variable['type']
            panel = Panel(plotset, title=title, dependant_title=axis_label, panel_type=panel_type)
            self.panels.append(panel)

    def __getVarLines__(self, aliases, ncdf_datasets, label="", line_format="", avg_axis=0, override_panel_type=None,
                        fallback_func=None, fill_zeros=False, lines=None, conversion_factor=1):
        """
        Get a list of Line objects for a specific clubb variable. If sam_file is specified it will also
        attempt to generate Lines for the SAM equivalent variables, using the name conversions found in
        VarnameConversions.py. If a SAM variable needs to be calculated (uses an equation) then it will have
        to be created within that variable group's dataset and not here.

        :param aliases: list of strings. Each string is the name of the variable to be plotted, case sensitive
        :param ncdf_datasets: List of Dataset objects containing clubb or sam netcdf data
        :param legend_label: Label to give the base-plotAll on the legend. This is normally Style_definitions.CLUBB_LABEL, but not provided as default to help avoid debugging confusion.
        :param line_format: Line formatting string used by matplotlib's PyPlot
        :param avg_axis: Axis over which to average values. 0 - time average, 1 - height average
        :param override_panel_type: Override the VariableGroup's default panel type
        :param fallback_func: If pyplotgen fails to find a variable by any of its aliases, this function will be called
            (if the fuction is provided) as a second attempt at getting the variable. fallback_func takes in a function
            as its data type (i.e. functional programming) and expects it to return data in the format (dependent_data, independent_data)
            where the data is a list of numeric values.
        :param fill_zeros: If a variable can't be found via it's aliases and any fallback or calc functions defined for that variable have failed,
        fill_zeros is an optional last resort to simply use a list of 0s in place of the variable for calculations.
        :param lines: Some plots require many lines to be plotted, such as budget plots. The lines parameter allows
        someone to write a dictionary containing a description of each line that needs to be plotted.
        For information on how to format the lines parameter, please see the config/VariableGroupBaseBudgets.py file and docs.
        :param conversion_factor: A multiplying factor used to scale clubb output. Defaults to 1.
        :return: A list of Line objects containing model data from whatever models were passed in.
        Returns None if requested variable is not found.
        """
        panel_type = self.defualt_panel_type
        if override_panel_type is not None:
            panel_type = override_panel_type
        all_lines = []

        if isinstance(ncdf_datasets, Dataset):
            ncdf_datasets = {'converted_to_dict': ncdf_datasets}
        line = None
        for dataset in ncdf_datasets.values():
            for varname in aliases:
                if varname not in dataset.variables.keys():
                    continue  # Skip loop if varname isn't in the dataset, this avoids some crashes
                if panel_type is Panel.TYPE_PROFILE:
                    line = self.__get_profile_line__(varname, dataset, label, line_format, conversion_factor, avg_axis,
                                                     fill_zeros)
                    break
                elif panel_type is Panel.TYPE_BUDGET:
                    # for overline in lines:
                    if varname in dataset.variables.keys():
                        budget_lines = self.__get_budget_lines__(lines, dataset, fill_zeros, label, line_format)
                        all_lines.extend(budget_lines)
                elif panel_type is Panel.TYPE_TIMESERIES:
                    line = self.__get_timeseries_line__(varname, dataset, self.end_time, conversion_factor, label,
                                                        line_format, fill_zeros)
                else:
                    raise ValueError(
                        'Invalid panel type ' + panel_type + '. Valid options are profile, budget, timeseries')
        if line is None and panel_type != Panel.TYPE_BUDGET and fill_zeros is False:
            datareader = DataReader()
            src_model = datareader.getNcdfSourceModel(dataset)
            print("\tFailed to find variable " + str(aliases) + " in " + src_model + " output for " + " case " + str(
                self.casename) +
                  ". Attempting to use fallback function.")
            line = self.__getVarDataFromFallback__(fallback_func, aliases, ncdf_datasets, label, line_format)
        if panel_type != Panel.TYPE_BUDGET and line is not None:
            all_lines.append(line)
        return all_lines

    def __get_profile_line__(self, varname, dataset, label, line_format, conversion_factor, avg_axis, fill_zeros):
        """
        Assumes variable can be plotted as a profile and returns a Line object
        representing the given variable for a profile plot. Profile plots are plots that show Height vs Varname

        :param varname: The name of the variable as a string
        :param dataset: NetCdf4 Dataset object containing the model output being plotted
        :param label: This is the name that will be shown on the legend for this line
        :param line_format: The pyplot line formated string representing how this line.
        Recommended value: emtpy string ""
        :param conversion_factor: This is a numerical value that will be multiplied element-wise to the variable. It's useful for doing
        basic model to model conversions, e.g. SAM -> CLUBB.
        :param avg_axis: Values will be averaged along this axis. This is basically only used for time-averaging profile plots.
        :param fill_zeros: If a variable can't be found via it's aliases and any fallback or calc functions defined for that variable have failed,
        fill_zeros is an optional last resort to simply use a list of 0s in place of the variable for calculations.
        :return: Line object representing the given variable for a profile plot
        """
        variable = NetCdfVariable(varname, dataset, start_time=self.start_time, end_time=self.end_time,
                                  avg_axis=avg_axis,
                                  conversion_factor=conversion_factor, fill_zeros=fill_zeros)
        if 'altitude' in dataset.variables.keys():
            z = NetCdfVariable('altitude', dataset, start_time=self.start_time, end_time=self.end_time)
        elif 'z' in dataset.variables.keys():
            z = NetCdfVariable('z', dataset, start_time=self.start_time, end_time=self.end_time)
        elif 'Z3' in dataset.variables.keys():
            z = NetCdfVariable('Z3', dataset, start_time=self.start_time, end_time=self.end_time)
        else:
            z = NetCdfVariable('lev', dataset, start_time=self.start_time, end_time=self.end_time)

        variable.constrain(self.height_min_value, self.height_max_value, data=z.data)
        z.constrain(self.height_min_value, self.height_max_value)
        line = Line(variable, z, label=label, line_format=line_format)
        return line

    def __get_timeseries_line__(self, varname, dataset, end_time, conversion_factor, label, line_format, fill_zeros):
        """

        :param varname: The name of the variable as a string
        :param dataset: NetCdf4 Dataset object containing the model output being plotted
        :param end_time: This is the last time value to plot. Timeseries plots always start at 0.
        :param label: This is the name that will be shown on the legend for this line
        :param line_format: The pyplot line formated string representing how this line.
        Recommended value: emtpy string ""
        :param conversion_factor: This is a numerical value that will be multiplied element-wise to the variable. It's useful for doing
        basic model to model conversions, e.g. SAM -> CLUBB.
        :param fill_zeros: If a variable can't be found via it's aliases and any fallback or calc functions defined for that variable have failed,
        fill_zeros is an optional last resort to simply use a list of 0s in place of the variable for calculations.
        :return: Line object representing the given variable for a timeseries plot
        """

        variable = NetCdfVariable(varname, dataset, start_time=0, end_time=end_time, avg_axis=1,
                                  conversion_factor=conversion_factor, fill_zeros=fill_zeros)
        time = NetCdfVariable('time', dataset)
        variable.constrain(0, variable.end_time, data=time.data)
        time.constrain(0, variable.end_time)
        line = Line(time, variable, label=label, line_format=line_format)
        return line

    def __get_budget_lines__(self, lines, dataset, fill_zeros, label, line_format):
        """

        :param variable:
        :param dataset:
        :param label:
        :return:
        """
        output_lines = []
        for line_definition in lines:
            line_added = False
            varnames = line_definition['aliases']
            for varname in varnames:
                if varname in dataset.variables.keys():
                    variable = NetCdfVariable(varname, dataset, start_time=self.start_time, end_time=self.end_time,
                                              fill_zeros=fill_zeros)
                    # if 'altitude' in dataset.variables.keys():
                    #     z = NetCdfVariable('altitude', dataset, start_time=self.start_time, end_time=self.end_time)
                    # else:
                    z = NetCdfVariable(['z', 'lev', 'altitude'], dataset, start_time=self.start_time, end_time=self.end_time)
                    variable.constrain(self.height_min_value, self.height_max_value, data=z.data)
                    z.constrain(self.height_min_value, self.height_max_value)
                    line_definition = Line(variable, z, label=line_definition['legend_label'],
                                           line_format="")  # uses auto-generating line format
                    output_lines.append(line_definition)
                    line_added = True
                    break

            if not line_added:
                print("\tFailed to find variable " + varname + " in case " + self.casename +
                      ". Attempting to use fallback function.")
                if 'fallback_func' in line_definition.keys():
                    fallback = line_definition['fallback_func']
                    fallback_output = self.__getVarDataFromFallback__(fallback, varname, {'budget': dataset}, label,
                                                                      line_format)
                    fallback_output = Line(fallback_output.x, fallback_output.y, line_format="", label=line_definition['legend_label'])
                    output_lines.append(fallback_output)

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
        return output_lines

    def __getVarDataFromFallback__(self, fallback, varname, datasets, label, line_format):
        """
        Some older model output doesn't contain all variables, e.g cgils_s12 doesn't contain WPTHLP.
        This method attempts to use a function under the name 'fallback_func' to generate the requested
        data similarly to how 'sam_calc' calculates sam output into a clubb format. If successful,
        returns varline data, otherwise raises error
        :return:
        """
        for dataset in datasets.values():
            try:
                varline, z = fallback(dataset_override=dataset)
                varline = Line(varline, z, line_format=line_format, label=label)
                print("\tFallback for ", varname, " successful")
                return varline
            except ValueError as e:
                # variable wasn't contained in dataset, try next dataset
                continue
            except TypeError as e:
                warn("Fallback failed for variable " + str(varname) + ". Skipping it.\n" + str(e))
                return None
        warn("Fallback failed for variable " + str(varname) + ". Skipping it.\n")
        return None

    def getVarForCalculations(self, varname, datasets, conversion_factor=1, fill_zeros=False):
        """
        This function is used within a fallback function to get the data of a certain variable,
        constrained between a min/max height.

        :param varname:
        :param include_z: If set to True, getVarForCalculations will return a tuple containing the data for both the varname
        variable and the height variable in a tuple ordered as (varname, z)

        :return:

        """
        if isinstance(datasets, Dataset):
            datasets = {'auto': datasets}
        for dataset in datasets.values():
            z_ncdf = NetCdfVariable(Case_definitions.HEIGHT_VAR_NAMES, dataset, 1)
            var_ncdf = NetCdfVariable(varname, dataset, conversion_factor, start_time=self.start_time,
                                      end_time=self.end_time, fill_zeros=fill_zeros)
            var_ncdf.constrain(self.height_min_value, self.height_max_value, data=z_ncdf.data)
            var_data = var_ncdf.data
            if np.amax(var_data) != 0:
                break  # Avoid overwriting vardata with auto-filled zeros
        return var_data
