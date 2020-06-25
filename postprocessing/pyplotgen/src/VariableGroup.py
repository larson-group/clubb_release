"""
:author: Nicolas Strike
:date: Mid 2019
"""
import math
import os
from pathlib import Path
from statistics import mean
from warnings import warn

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

    def __init__(self, case, clubb_datasets=None, les_dataset=None, coamps_dataset=None, r408_dataset=None,
                 hoc_dataset=None, cam_datasets=None,
                 e3sm_datasets=None, sam_datasets=None, wrf_datasets=None):
        """
        Initialize a VariableGroup object with the passed parameters

        :param case: An instance of a Case object
        :param clubb_datasets: A dictionary of or a single NetCDF4 Dataset object
            containing the dependent_data to be plotted
        :param les_dataset: NetCDF4 Dataset object containing sam ouput
        :param coamps_dataset: NetCDF4 Dataset object containing coamps ouput
        :param r408_dataset: NetCDF4 Dataset object containing clubb r408 ('chris golaz') output
        :param hoc_dataset: NetCDF4 Dataset object containing hoc output
        :param cam_datasets: A dictionary of or a single NetCDF4 Dataset object containing cam output
        :param e3sm_datasets: A dictionary of or a single NetCDF4 Dataset object containing e3sm output
        :param sam_datasets: A dictionary of or a single NetCDF4 Dataset object containing sam output
        :param wrf_datasets: A dictionary of or a single NetCDF4 Dataset object containing wrf output
        """
        self.variables = []
        self.panels = []
        self.default_panel_type = Panel.TYPE_PROFILE
        self.clubb_datasets = clubb_datasets
        self.case = case
        self.les_dataset = les_dataset
        self.e3sm_datasets = e3sm_datasets
        self.sam_datasets = sam_datasets
        self.cam_datasets = cam_datasets
        self.wrf_datasets = wrf_datasets
        self.coamps_dataset = coamps_dataset
        self.r408_datasets = r408_dataset
        self.hoc_datasets = hoc_dataset
        self.casename = case.name
        self.start_time = case.start_time
        self.end_time = case.end_time
        self.height_min_value = case.height_min_value
        self.height_max_value = case.height_max_value

        # Loop over the list self.variable_definitions which is only defined in the subclasses
        # that can be found in the config folder such as VariableGroupBase
        for variable in self.variable_definitions:
            print("\tProcessing ", variable['var_names']['clubb'])

            # Only add variable if none of the var_names are blacklisted
            all_var_names = []
            for model_var_names in variable['var_names'].values():
                all_var_names.extend(model_var_names)
            variable_is_blacklisted = len(list(set(all_var_names).intersection(case.blacklisted_variables))) != 0

            if not variable_is_blacklisted:
                self.addVariable(variable)
            else:
                print('\tVariable {} is blacklisted and will therefore not be plotted.'.format(variable))

        self.generatePanels()

    def addVariable(self, variable_def_dict):
        """
        Given basic details about a variable, this
        creates lines/panels for the variable and then adds that
        variable to the list of all variables for this VariableGroup

        :param variable_def_dict: A dict containing the information defining a variable. These definitions are declared
            inside a VariableGroup child class, e.g. VariableGroupBase.

        Valid dict keys/options include:

        *var_names*: A list of names various models refer to this variable as. E.g. ['wprtp', 'WPRTP', 'wpqtp']

        *sam_calc*: (optional) A functional reference to a method that calculates a sam variable. This is given as the
            name of the function *without*  a set of parenthesis () after the name. E.g. self.getThlmSamCalc

        *coamps_calc*: (optional) A functional reference to a method that calculates a coamps variable. This is given
            as the name of the function *without* a set of parenthesis () after the name. E.g. self.getThlmCoampsCalc

        *r408_calc*: (optional) A functional reference to a method that calculates a r408 variable. This is given as
            the name of the function *without*  a set of parenthesis () after the name. E.g. self.getTHlmR408Calc

        *hoc_calc*: (optional) A functional reference to a method that calculates a hoc-2005 variable. This is given
            as the name of the function *without* a set of parenthesis () after the name. E.g. self.getThlmHocCalc

        *e3sm_calc*: (optional) A functional reference to a method that calculates an e3sm variable. This is given
            as the name of the function *without* a set of parenthesis () after the name. E.g. self.getThlmE3smCalc

        *sam_conv_factor*: (optional) Numeric value to scale a sam variable by. E.g. 1/1000, or 100

        *coamps_conv_factor*: (optional) Numeric value to scale a coamps variable by. E.g. 1/1000, or 100

        *r408_conv_factor*: (optional) Numeric value to scale a clubb r408 variable by. E.g. 1/1000, or 100

        *e3sm_conv_factor*: (optional) Numeric value to scale an e3sm variable by. E.g. 1/1000, or 100

        *type*: (optional) Override the default type 'profile' with either 'budget' or 'timeseries'

        *fallback_func*: (optional) (DEPRECATED: Use the calc parameters instead)
            If a variable is not found within a dataset_name and a fallback_func is specified, this
            function will be called to attempt retrieving the variable. Like the
            model_calc option, this is a functional reference to a method that calculates the given variable.
            E.g. self.getWpthlpFallback

        *title*: (optional) Override the default panel title, or provide one if it's not specified in the netcdf file.

        *axis_title*: (optional) Override the default dependent axis title, or provide one if it's not
            specified in the netcdf file.

        *lines*: * (budget variables only) Defines lines to plot for budget cases. Passed seperately because
            it's a lot of text.
            This is given in the form of a list of lines, here's an example:

        .. code-block:: python
            :linenos:

            thlm_lines = [
            {'var_names': ['thlm_bt'], 'legend_label': 'thlm_bt'},
            {'var_names': ['thlm_ma'], 'legend_label': 'thlm_ma'},
            {'var_names': ['thlm_ta'], 'legend_label': 'thlm_ta'},
            {'var_names': ['thlm_mc'], 'legend_label': 'thlm_mc'},
            {'var_names': ['thlm_clipping'], 'legend_label': 'thlm_bt', 'clubb_calc': self.getThlmClipping},
            {'var_names': ['radht'], 'legend_label': 'radht'},
            {'var_names': ['lsforcing'], 'legend_label': 'lsforcing', 'clubb_calc': self.getLsforcing},
            {'var_names': ['thlm_residual'], 'legend_label': 'thlm_residual', 'clubb_cal': self.getThlmResidual},
            ]

        Here are a couple examples of other (non-budget) variable definitions:

        .. code-block:: python
            :linenos:

            {
            'var_names': ['Skrt_zt'],
            'sam_calc': self.getSkrtZtLesCalc,
            'coamps_calc': self.getSkrtZtLesCalc,
            }

            {
            'var_names': ['lwp', 'CWP'],
            'type': Panel.TYPE_TIMESERIES,
            'sam_conv_factor': 1/1000
            },

        :return: None
        """
        var_names = variable_def_dict['var_names']

        plot_les = self.les_dataset is not None
        plot_coamps = self.coamps_dataset is not None
        plot_r408 = self.r408_datasets is not None
        plot_hoc = self.hoc_datasets is not None

        plot_sam = self.sam_datasets is not None and len(self.sam_datasets) > 0 \
                   and (len(var_names['sam']) > 0 or 'sam_calc' in variable_def_dict.keys())
        plot_e3sm = self.e3sm_datasets is not None and len(self.e3sm_datasets) > 0 \
                    and (len(var_names['e3sm']) > 0 or 'e3sm_calc' in variable_def_dict.keys())
        plot_clubb = self.clubb_datasets is not None and len(self.clubb_datasets) > 0 \
                     and (len(var_names['clubb']) > 0 or 'clubb_calc' in variable_def_dict.keys())
        plot_wrf = self.wrf_datasets is not None and len(self.wrf_datasets) > 0 \
                   and (len(var_names['wrf']) > 0 or 'wrf_calc' in variable_def_dict.keys())
        plot_cam = self.cam_datasets is not None and len(self.cam_datasets) > 0 \
                   and (len(var_names['cam']) > 0 or 'cam_calc' in variable_def_dict.keys())

        all_lines = []
        # Plot benchmarks
        if plot_les:
            all_lines.extend(self.__getVarLinesForModel__('sam', variable_def_dict, self.les_dataset))
        if plot_coamps:
            all_lines.extend(self.__getVarLinesForModel__('coamps', variable_def_dict, self.coamps_dataset))
        if plot_r408:
            all_lines.extend(self.__getVarLinesForModel__('r408', variable_def_dict, self.r408_datasets))
        if plot_hoc:
            all_lines.extend(self.__getVarLinesForModel__('hoc', variable_def_dict, self.hoc_datasets))

        # Plot input folders
        if plot_clubb:
            for input_folder in self.clubb_datasets:  # TODO this loop is causing extra budget lines
                folder_name = os.path.basename(input_folder)
                all_lines.extend(
                    self.__getVarLinesForModel__('clubb', variable_def_dict, self.clubb_datasets[input_folder],
                                                 label=folder_name))

        if plot_sam:
            for input_folder in self.sam_datasets:
                folder_name = os.path.basename(input_folder)
                all_lines.extend(self.__getVarLinesForModel__('sam', variable_def_dict, self.sam_datasets[input_folder],
                                                              label=folder_name))
        if plot_e3sm:
            for input_folder in self.e3sm_datasets:
                folder_name = os.path.basename(input_folder)
                all_lines.extend(
                    self.__getVarLinesForModel__('e3sm', variable_def_dict, self.e3sm_datasets[input_folder],
                                                 label=folder_name))

        if plot_cam:
            for input_folder in self.cam_datasets:
                folder_name = os.path.basename(input_folder)
                all_lines.extend(
                    self.__getVarLinesForModel__('cam', variable_def_dict, self.cam_datasets[input_folder],
                                                 label=folder_name))

        if plot_wrf:
            for input_folder in self.wrf_datasets:
                folder_name = os.path.basename(input_folder)
                all_lines.extend(self.__getVarLinesForModel__('wrf', variable_def_dict, self.wrf_datasets[input_folder],
                                                              label=folder_name))

        variable_def_dict['plots'] = all_lines

        plotted_models_varname = self.__getRelativeVarName__(var_names)
        title, axis_title = self.__getTitles__(variable_def_dict, plotted_models_varname)
        variable_def_dict['title'] = title
        variable_def_dict['axis_title'] = axis_title

        if len(variable_def_dict['plots']) > 0:
            self.variables.append(variable_def_dict)

    def __getVarLinesForModel__(self, model_name, variable_def_dict, dataset, label="no label found"):
        """
        Generates and returns lines for plotting on a Panel for a single model (e.g. sam, clubb, e3sm).
        This function can produce lines for either benchmarks or user-inputted folders.
        To get lines for a benchmark, pass in the 3 required parameters. To plot an clubb folder, the label
        parameter must also be specified, as this is the label that will appear on the panel's legend. Benchmarks
        do not need the label parameter because their labels are specified in the Style_definitions.py file.

        :param model_name: The name of the model to be plotted.
            This needs to be the same as the name used to specify
            the model throughout the various dictionary keys in pyplotgen.
            E.g. sam output is specified with 'sam', where sam calclulation functions are denoted by 'sam_calc'
            and Style_definitions.py' uses the key 'sam' to refer to sam configurations such as BENCHMARK_LABELS['sam']
        :param variable_def_dict: This is the variable defining dict, e.g. from VariableGroupBase.py
        :param dataset: Either a NetCDF Dataset object or a dict of the format {'some key', Dataset}
        :param label: The label to be used on the panel's legend.
            This is a required field when getting lines for a user-inputted folder,
            but not a required field for processing benchmark output.
        :return: An array of lines to be added to a panel's plots (via the extend() function)
        """
        all_model_var_names = variable_def_dict['var_names']
        conv_factors = self.__getConvFactors__(variable_def_dict)

        lines = None
        if 'lines' in variable_def_dict.keys():
            lines = variable_def_dict['lines']

        panel_type = self.default_panel_type
        if 'type' in variable_def_dict.keys():
            panel_type = variable_def_dict['type']

        if not isinstance(dataset, dict):
            dataset = {model_name: dataset}

        line_style = ""
        if label == "no label found":
            label = Style_definitions.BENCHMARK_LABELS[model_name]
            line_style = Style_definitions.BENCHMARK_LINE_STYLES[model_name]

        all_lines = []

        data_was_calculated = False
        if (model_name + '_calc') in variable_def_dict.keys() and not self.__varnamesInDataset__(
                all_model_var_names[model_name], dataset):
            plot_data, z = variable_def_dict[(model_name + '_calc')](dataset_override=dataset)
            plot = Line(plot_data, z, line_format=line_style, label=label)
            all_lines.append(plot)
            data_was_calculated = True

        if lines is not None:
            for line in lines:
                if (model_name + '_calc') in line.keys() and not self.__varnamesInDataset__(line['var_names'], dataset):
                    plot_data, z = line[(model_name + '_calc')](dataset_override=dataset)
                    plot = Line(plot_data, z, line_format=line_style, label=line['legend_label'])
                    all_lines.append(plot)

        if len(all_model_var_names[model_name]) > 0 and not data_was_calculated:
            all_lines.extend(self.__getVarLines__(all_model_var_names[model_name], dataset,
                                                  conversion_factor=conv_factors[model_name],
                                                  label=label,
                                                  line_format=line_style,
                                                  override_panel_type=panel_type,
                                                  lines=lines))
        return all_lines

    def getTextDefiningDataset(self):
        """
        Finds a dataset from which text names and descriptions can be pulled from.
        By default, it is always clubb, however if pyplotgen is not plotting clubb then
        it will return the first dataset it can find.

        :return: A list of all Dataset objects for all models
        """
        datasets = []
        if self.clubb_datasets is not None:
            for i in self.clubb_datasets.values():
                if isinstance(i, dict):
                    for dataset in i.values():
                        datasets.append(dataset)
                elif isinstance(i, Dataset):
                    datasets.append(i)
                else:
                    raise TypeError("getTextDefiningDataset received an unexpected format of datasets")
        if self.e3sm_datasets is not None:
            for dataset in self.e3sm_datasets.values():
                datasets.append(dataset)
        if self.cam_datasets is not None:
            for dataset in self.cam_datasets.values():
                datasets.append(dataset)
        if self.sam_datasets is not None:
            for dataset in self.sam_datasets.values():
                datasets.append(dataset)
        if self.wrf_datasets is not None:
            for i in self.wrf_datasets.values():
                if isinstance(i, dict):
                    datasets.extend(i.values())
                elif isinstance(i, Dataset):
                    datasets.append(i)
                else:
                    raise TypeError("getTextDefiningDataset recieved an unexpected format of datasets")
        if self.les_dataset is not None:
            datasets.append(self.les_dataset)
        if self.coamps_dataset is not None:
            datasets.extend(self.coamps_dataset.values())
        if self.r408_datasets is not None:
            datasets.extend(self.r408_datasets.values())
        if self.hoc_datasets is not None:
            datasets.extend(self.hoc_datasets.values())
        if len(datasets) == 0:
            raise FileExistsError("No dataset could be found to pull text descriptions from")
        return datasets

    def generatePanels(self):
        """
        Generates a set of panels from the plots stored in self.
        Does not return anything, simply assigns the panels into self.panels

        :return: None
        """
        for variable in self.variables:
            title = variable['title']
            axis_label = variable['axis_title']
            plotset = variable['plots']
            panel_type = self.default_panel_type
            centered = False
            if 'type' in variable.keys():
                panel_type = variable['type']
            if 'sci_scale' in variable.keys():
                sci_scale = variable['sci_scale']
            else:
                sci_scale = None
            if 'centered' in variable.keys():
                centered = variable['centered']
            panel = Panel(plotset, title=title, dependent_title=axis_label,
                          panel_type=panel_type, sci_scale=sci_scale, centered=centered)
            self.panels.append(panel)

    def __getTitles__(self, variable_def_dict, plotted_models_varname):
        """
        Creates and returns figure and axis titles of the panel based on panel type and
        panel definition in variable_def_dict

        :param variable_def_dict: Definition of the variable being used
        :param plotted_models_varname: String containing the name of the plotted variable
        :return: title, axis_title
        """
        title = "Title not found"
        axis_title = "axis_title not found"

        data_reader = DataReader()
        var_names = variable_def_dict['var_names']
        panel_type = self.default_panel_type
        all_lines = variable_def_dict['plots']
        if 'type' in variable_def_dict.keys():
            panel_type = variable_def_dict['type']
        try:
            first_input_datasets = self.getTextDefiningDataset()
        except FileExistsError as e:
            warn(str(e))
            return title, axis_title

        if 'title' not in variable_def_dict.keys():
            # No title given so it must be generated from given information
            if panel_type == Panel.TYPE_BUDGET:
                source_folder = os.path.basename(Path(first_input_datasets[0].filepath()).parent)
                title = source_folder + ' ' + plotted_models_varname
            else:
                all_var_names = []
                for model_var_names in var_names.values():
                    all_var_names.extend(model_var_names)
                # Get long name for any of the var_names given in variable_def_dict
                imported_title = data_reader.getLongName(first_input_datasets, all_var_names)
                title = imported_title
        else:
            # A title is manually defined in variable_def_dict
            title = variable_def_dict['title']

        if 'axis_title' not in variable_def_dict.keys():
            # No axis_title given so it must be generated from given information
            if panel_type == Panel.TYPE_BUDGET and len(all_lines) > 0:
                any_varname_with_budget_units = [var.label for var in all_lines]
                axis_title = "[" + data_reader.__getUnits__(first_input_datasets, any_varname_with_budget_units) + "]"
            else:
                # Get axis title for any of the var_names given in variable_def_dict
                imported_axis_title = data_reader.getAxisTitle(first_input_datasets, plotted_models_varname)
                axis_title = imported_axis_title
        else:
            # An axis_title is manually defined in variable_def_dict
            axis_title = variable_def_dict['axis_title']

        return title, axis_title

    def __varnamesInDataset__(self, varnames, datasets):
        """
        Returns True if the dataset contains a variable with one of the given varnames, False otherwise

        :param varnames: A list of possible variable names
        :param datasets: Either the dataset being investigated (NetCDF Dataset object) or
            a dict of datasets to be compared
        :return: True if a name is found, False otherwise
        """
        if not isinstance(datasets, dict):
            datasets = {'auto': datasets}
        for name in varnames:
            for dataset in datasets.values():
                if name in dataset.variables.keys():
                    return True
        return False

    def __getRelativeVarName__(self, var_names):
        """
        Returns the varname for a variable relative to the models that were being plotted

        :param var_names: The dict of various names for the given variable
        :return: A relevant name of the given variable
        """
        plotted_models_varname = "unknown_model_var"
        if self.clubb_datasets is not None and len(var_names['clubb']) > 0:
            plotted_models_varname = var_names['clubb'][0]
        elif self.clubb_datasets is None and self.wrf_datasets is not None and len(var_names['wrf']) > 0:
            plotted_models_varname = var_names['wrf'][0]
        elif self.clubb_datasets is None and self.sam_datasets is not None and len(var_names['sam']) > 0:
            plotted_models_varname = var_names['sam'][0]
        elif self.clubb_datasets is None and self.e3sm_datasets is not None and len(var_names['e3sm']) > 0:
            plotted_models_varname = var_names['e3sm'][0]
        return plotted_models_varname

    def __getConvFactors__(self, variable_def_dict):
        """
        This is a helper method that loads parameters from the variables
        definition and assigns them to python variables accordingly.

        :param variable_def_dict: Definition of the variable being used
        :return: A dict containing the various conversion factors for the inputted variable
        """
        conv_factors = {
            'les': 1,
            'sam': 1,
            'wrf': 1,
            'coamps': 1,
            'r408': 1,
            'hoc': 1,
            'e3sm': 1,
            'cam': 1,
            'clubb': 1
        }

        if 'sam_conv_factor' in variable_def_dict.keys():
            conv_factors['les'] = variable_def_dict['sam_conv_factor']
            conv_factors['sam'] = variable_def_dict['sam_conv_factor']

        if 'coamps_conv_factor' in variable_def_dict.keys():
            conv_factors['coamps'] = variable_def_dict['coamps_conv_factor']

        if 'r408_conv_factor' in variable_def_dict.keys():
            conv_factors['r408'] = variable_def_dict['r408_conv_factor']

        if 'hoc_conv_factor' in variable_def_dict.keys():
            conv_factors['hoc'] = variable_def_dict['hoc_conv_factor']

        if 'e3sm_conv_factor' in variable_def_dict.keys():
            conv_factors['e3sm'] = variable_def_dict['e3sm_conv_factor']

        if 'cam_conv_factor' in variable_def_dict.keys():
            conv_factors['cam'] = variable_def_dict['cam_conv_factor']

        if 'wrf_conv_factor' in variable_def_dict.keys():
            conv_factors['wrf'] = variable_def_dict['wrf_conv_factor']

        return conv_factors

    def __getVarLines__(self, var_names, ncdf_datasets, label="", line_format="", avg_axis=0, override_panel_type=None,
                        lines=None, conversion_factor=1):
        """
        Get a list of Line objects for a specific clubb variable. If les_dataset is specified it will also
        attempt to generate Lines for the SAM equivalent variables, using the name conversions found in
        VarnameConversions.py. If a SAM variable needs to be calculated (uses an equation) then it will have
        to be created within that variable group's dataset and not here.

        :param var_names: A string or list of strings.
            Each string is the name of the variable to be plotted, case sensitive
        :param ncdf_datasets: A list of Dataset objects containing clubb or sam netcdf dependent_data
        :param label: Label to give the base-plotAll on the legend. This is normally
            Style_definitions.CLUBB_LABEL, but not provided as default to help avoid debugging confusion.
        :param line_format: Line formatting string used by matplotlib's PyPlot
        :param avg_axis: Axis over which to average values. 0 - time average, 1 - height average
        :param override_panel_type: Override the VariableGroup's default panel type
        :param lines: Some plots require many lines to be plotted, such as budget plots. The lines parameter allows
            someone to write a dictionary containing a description of each line that needs to be plotted.
            For information on how to format the lines parameter,
            please see the config/VariableGroupBaseBudgets.py file and docs.
        :param conversion_factor: A multiplying factor used to scale clubb output. Defaults to 1.
        :return: A list of Line objects containing model dependent_data from whatever models were passed in.
            Returns None if requested variable is not found.
        """
        panel_type = self.default_panel_type
        if override_panel_type is not None:
            panel_type = override_panel_type
        all_lines = []

        if isinstance(ncdf_datasets, Dataset):
            ncdf_datasets = {'converted_to_dict': ncdf_datasets}
        line = None
        if panel_type is Panel.TYPE_PROFILE:
            profile_lines = self.__get_profile_line__(var_names, ncdf_datasets, label, line_format, conversion_factor,
                                                      avg_axis, lines=lines)
            if profile_lines is not None:
                all_lines.extend(profile_lines)
        elif panel_type is Panel.TYPE_BUDGET:
            # for dataset in ncdf_datasets.values():
            budget_lines = self.__get_budget_lines__(lines, ncdf_datasets)
            all_lines.extend(budget_lines)
        elif panel_type is Panel.TYPE_TIMESERIES:
            line = self.__get_timeseries_line__(var_names, ncdf_datasets, self.end_time, conversion_factor, label,
                                                line_format)
            if line is not None:
                all_lines.append(line)
        elif panel_type:
            raise ValueError('Invalid panel type ' + panel_type + '. Valid options are profile, budget, timeseries')

        return all_lines

    def __get_profile_line__(self, varname, dataset, label, line_format, conversion_factor, avg_axis, lines=None):
        """
        Assumes variable can be plotted as a profile and returns a Line object
        representing the given variable for a profile plot. Profile plots are plots that show Height vs Varname

        :param varname: The name of the variable as a string
        :param dataset: NetCdf4 Dataset object containing the model output being plotted
        :param label: This is the name that will be shown on the legend for this line
        :param line_format: A string representing how this line should be formatted using pyplot formatting.
            Recommended value: emtpy string ""
        :param conversion_factor: This is a numerical value that will be multiplied element-wise to the variable.
            It's useful for doing basic model to model conversions, e.g. SAM -> CLUBB.
        :param avg_axis: Values will be averaged along this axis. This is basically only used for time-averaging
            profile plots.
        :return: Line object representing the given variable for a profile plot
        """
        output_lines = []
        variable = NetCdfVariable(varname, dataset, independent_var_names=Case_definitions.HEIGHT_VAR_NAMES,
                                  start_time=self.start_time, end_time=self.end_time,
                                  avg_axis=avg_axis,
                                  conversion_factor=conversion_factor)
        variable.trimArray(self.height_min_value, self.height_max_value, data=variable.independent_data.data)

        if lines is not None:
            additional_lines = self.__processLinesParameter__(lines, dataset, line_format=line_format, label_suffix=label)
            output_lines.extend(additional_lines)
        else:
            line = Line(variable, variable.independent_data, label=label, line_format=line_format)
            output_lines.append(line)
        return output_lines

    def __get_timeseries_line__(self, varname, dataset, end_time, conversion_factor, label, line_format):
        """
        Returns a Line object for a timeseries plot of the variable varname in dataset.

        :param varname: The name of the variable as a string
        :param dataset: NetCdf4 Dataset object containing the model output being plotted
        :param end_time: This is the last time value to plot. Timeseries plots always start at 0.
        :param conversion_factor: This is a numerical value that will be multiplied element-wise to the variable.
            It's useful for doing basic model to model conversions, e.g. SAM -> CLUBB.
        :param label: This is the name that will be shown on the legend for this line
        :param line_format: A string representing how this line should be formatted using pyplot formatting.
            Recommended value: emtpy string ""
        :return: Line object representing the given variable for a timeseries plot
        """
        variable = NetCdfVariable(varname, dataset, independent_var_names=Case_definitions.TIME_VAR_NAMES, start_time=0,
                                  end_time=end_time, avg_axis=1,
                                  conversion_factor=conversion_factor)
        variable.trimArray(0, variable.end_time, data=variable.independent_data.data)
        line = Line(variable.independent_data, variable, label=label, line_format=line_format)
        return line

    def __get_budget_lines__(self, lines, dataset):
        """
        Returns a list of Line objects for a budget plot for each variable defined in lines

        :param lines: A list of line definitions containing the variable names, legend labels etc.
            See the addVariable documentation above for a complete description
        :param dataset: NetCdf4 Dataset object containing the model output being plotted
        :param label: This is the name that will be shown on the legend for this line (currently not used here)
        :param line_format: A string representing how this line should be formatted using pyplot formatting.
            Recommended value: emtpy string "" (currently not used here)
        :return: A list of Line objects for a budget plot derived from lines
        """
        output_lines = self.__processLinesParameter__(lines, dataset)
        return output_lines

    def getVarForCalculations(self, varname, datasets, conversion_factor=1):
        """
        This function is used within model_calc functions to get the dependent_data of a certain variable,
        constrained between a min/max height.

        :param varname: The name of the variable as a string.
        :param datasets: A netCDF4 Dataset object or dict of such containing the model output being plotted.
        :param conversion_factor: This is a numerical value that will be multiplied element-wise to the variable.
            It's useful for doing basic model to model conversions, e.g. SAM -> CLUBB.
        :return: A tuple containing the dependent_data for the variable, the height data, and the datasets
        """
        var_ncdf = NetCdfVariable(varname, datasets, independent_var_names=Case_definitions.HEIGHT_VAR_NAMES,
                                  conversion_factor=conversion_factor, start_time=self.start_time,
                                  end_time=self.end_time)
        var_ncdf.trimArray(self.height_min_value, self.height_max_value, data=var_ncdf.independent_data)
        var_data = var_ncdf.dependent_data
        z_data = var_ncdf.independent_data
        return var_data, z_data, datasets

    def isSurfaceData(self, dataset):
        """
        Returns False if the given dataset contains a height var, True otherwise

        :param dataset: A netCDF4 Dataset object
        :return: False if a dataset contains a height var, True otherwise
        """
        for height_var in Case_definitions.HEIGHT_VAR_NAMES:
            if height_var in dataset.variables.keys():
                return False
        return True

    def pickNonZeroOutput(self, output1, output2):
        """
        Sometimes there are more than 1 ways to calculate a variable, and
        it is impossible to know which equation to use before hand. In these
        cases, all equations should be used and their output given to this
        function. This function will then attempt to determine which of the
        outputs from the equations is most likely to be the expected answer.

        This function currenlty works by determining picking the non-nan and non-zero array.
        In the event that both arrays contain non-nan and non-zero output, it will return output1.

        Example usage:
        def get_rc_coef_zm_X_wprcp_coamps_calc(self, dataset_override=None):
            wprlp, z, dataset = self.getVarForCalculations(['thlpqcp', 'wpqcp', 'wprlp'], dataset)
            ex0, z, dataset = self.getVarForCalculations(['ex0'], dataset)
            p, z, dataset = self.getVarForCalculations('p', dataset)
            thvm, z, dataset = self.getVarForCalculations('thvm', dataset)
            output1 = wprlp * (2.5e6 / (1004.67 * ex0) - 1.61 * thvm)
            output2 = wprlp * (2.5e6 / (1004.67 * ((p / 1.0e5) ** (287.04 / 1004.67))) - 1.61 * thvm)
            output = self.pickMostLikelyOutputList(output1, output2)
            return output, z


        :param output1: An arraylike list of numbers outputted from a variable calculation
            (e.g.     output1 = wprlp * (2.5e6 / (1004.67 * ex0) - 1.61 * thvm)
        :param output2: An arraylike list of numbers outputted from a variable calculation
            (e.g.     output2 = wprlp * (2.5e6 / (1004.67 * ((p / 1.0e5) ** (287.04 / 1004.67))) - 1.61 * thvm)
        :return: The listlike array of the two datasets most likely to be the correct answer.
        """
        output1_mean = mean(output1)
        output2_mean = mean(output2)

        if math.isnan(output1_mean):
            return output2
        if math.isnan(output2_mean):
            return output1
        if output1_mean == 0:
            return output2
        elif output2_mean == 0:
            return output1
        raise UserWarning("Failed to find an easy answer to the best output.")
        # return output1

    def __processLinesParameter__(self, lines, dataset, label_suffix="", line_format = ""):
        """
        This method processes a 'lines' parameter from a variable definition and translates
        it into a list of Line objects for plotting.

        :param lines: a lines definition. See VariableGroupBaseBudgets.py for examples
        :param dataset: A netcdf Dataset object or a dict of datasets
        :return: list of Line objects for plotting
        """
        output_lines = []
        for line_definition in lines:
            varnames = line_definition['var_names']
            label = line_definition['legend_label'] + " " + label_suffix
            variable = NetCdfVariable(varnames, dataset, independent_var_names=Case_definitions.HEIGHT_VAR_NAMES,
                                      start_time=self.start_time, end_time=self.end_time)
            variable.trimArray(self.height_min_value, self.height_max_value, data=variable.independent_data)
            line_definition = Line(variable, variable.independent_data, label=label,
                                   line_format=line_format)  # uses auto-generating line format
            output_lines.append(line_definition)
        return output_lines
