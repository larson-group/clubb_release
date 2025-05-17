"""
:author: Nicolas Strike
:date: Mid 2019
"""
import math
import os
from pathlib import Path
from statistics import mean

import numpy as np
from netCDF4._netCDF4 import Dataset

from config import Case_definitions
from config import Style_definitions
from src.Contour import Contour
from src.ContourPanel import ContourPanel
from src.DataReader import DataReader, NetCdfVariable
from src.Line import Line
from src.Panel import Panel
from src.AnimationPanel import AnimationPanel
from src.OutputHandler import logToFile, logToFileAndConsole

class VariableGroup:
    """
    This is the parent VariableGroup class. All other panel groups
    should be created as a subclass of this one. Properties and
    methods common to each PanelGroup can be found here.

    A VariableGroup child defines Lines, Panels, etc., and is responsible for
    calculating any 'calculated' variables from netcdf

    For information on the input parameters of this class, please see the documentation for the
    ``__init__()`` method.
    """

    def __init__(self, case, clubb_datasets=None, sam_benchmark_dataset=None, coamps_benchmark_dataset=None,
                 wrf_benchmark_dataset=None, r408_dataset=None,
                 hoc_dataset=None, cam_datasets=None, e3sm_datasets=None, sam_datasets=None, wrf_datasets=None,
                 priority_vars=False, background_rcm=False, background_rcm_folder=None):
        """
        Initialize a VariableGroup object with the passed parameters

        :param case: An instance of a CaseGallerySetup object
        :param clubb_datasets: A dictionary of or a single NetCDF4 Dataset object
            containing the dependent_data to be plotted
        :param sam_benchmark_dataset: NetCDF4 Dataset object containing sam ouput
        :param coamps_benchmark_dataset: NetCDF4 Dataset object containing coamps ouput
        :param wrf_benchmark_dataset: NetCDF4 Dataset object containing WRF-LASSO output
        :param r408_dataset: NetCDF4 Dataset object containing clubb r408 ('chris golaz') output
        :param hoc_dataset: NetCDF4 Dataset object containing hoc output
        :param cam_datasets: A dictionary of or a single NetCDF4 Dataset object containing cam output
        :param e3sm_datasets: A dictionary of or a single NetCDF4 Dataset object containing e3sm output
        :param sam_datasets: A dictionary of or a single NetCDF4 Dataset object containing sam output
        :param wrf_datasets: A dictionary of or a single NetCDF4 Dataset object containing wrf output
        :param background_rcm: Show a height-based "contour" plot of time-averaged rcm behind CLUBB profiles.
        :param background_rcm_folder: Folder of CLUBB output data from which to pull the background rcm (if background_rcm is enabled).
        """
        self.variables = []
        self.panels = []
        self.default_panel_type = Panel.TYPE_PROFILE
        self.clubb_datasets = clubb_datasets
        self.case = case
        self.sam_benchmark_dataset = sam_benchmark_dataset
        self.coamps_benchmark_dataset = coamps_benchmark_dataset
        self.wrf_benchmark_dataset = wrf_benchmark_dataset
        self.e3sm_datasets = e3sm_datasets
        self.sam_datasets = sam_datasets
        self.cam_datasets = cam_datasets
        self.wrf_datasets = wrf_datasets
        self.r408_datasets = r408_dataset
        self.hoc_datasets = hoc_dataset
        self.casename = case.name
        self.start_time = case.start_time
        self.end_time = case.end_time
        self.height_min_value = case.height_min_value
        self.height_max_value = case.height_max_value
        self.time_height = case.time_height
        self.animation = case.animation
        self.priority_vars = priority_vars
        self.bkgrnd_rcm_flag = background_rcm
        self.bkgrnd_rcm_folder = background_rcm_folder

        # Loop over the list self.variable_definitions which is only defined in the subclasses
        # that can be found in the config folder such as VariableGroupBase
        for variable in self.variable_definitions:
            logToFile("\tProcessing {}".format(variable['var_names']['clubb']))
            # Only add variable if none of the var_names are blacklisted
            all_var_names = []
            for model_var_names in variable['var_names'].values():
                all_var_names.extend(model_var_names)
            variable_is_blacklisted = len(list(set(all_var_names).intersection(case.blacklisted_variables))) != 0

            if not variable_is_blacklisted:
                self.addVariable(variable)
            else:
                logToFile('\tVariable {} is blacklisted and will therefore not be plotted.'.format(variable))

        if self.bkgrnd_rcm_flag:
            if self.clubb_datasets is not None and len(self.clubb_datasets) != 0:
                # Extract rcm from the zt NetCDF file. Also extract the time and height values to which the
                # rcm data points correspond.
                bkgrnd_rcm = np.squeeze( self.clubb_datasets[self.bkgrnd_rcm_folder]['zt'].variables['rcm'] )
                self.altitude_bkgrnd_rcm = np.squeeze( self.clubb_datasets[self.bkgrnd_rcm_folder]['zt'].variables['altitude'] )
                self.time_bkgrnd_rcm = np.squeeze( self.clubb_datasets[self.bkgrnd_rcm_folder]['zt'].variables['time'] )
                # Find the indices in the rcm data that correspond to the start time and end time requested as the
                # time-averaging interval for the case, as well as the minimum height and maximum height requested
                # for the plots.
                start_time_seconds = 60.0 * self.start_time # self.start_time is in minutes, while time_bkgrnd_rcm is in seconds.
                end_time_seconds = 60.0 * self.end_time # self.end_time is in minutes, while time_bkgrnd_rcm is in seconds.
                start_time_idx, end_time_idx = DataReader.__getStartEndIndex__(self.time_bkgrnd_rcm, start_time_seconds, end_time_seconds)
                self.start_alt_idx, self.end_alt_idx = DataReader.__getStartEndIndex__(self.altitude_bkgrnd_rcm, self.height_min_value, self.height_max_value)
                if self.animation:
                    self.bkgrnd_rcm = bkgrnd_rcm
                else:
                    # Calculate the time-averaged vertical profile of rcm for use as contours in the background of plots
                    # of CLUBB time-averaged vertical profiles of various model fields.
                    nzt = np.shape(bkgrnd_rcm)[1]
                    self.bkgrnd_rcm = np.zeros(nzt)
                    for z_indx in range(nzt):
                        lev_sum = 0
                        count = 0
                        for t_indx in range(start_time_idx, end_time_idx):
                            lev_sum = lev_sum + bkgrnd_rcm[t_indx,z_indx]
                            count = count + 1
                        self.bkgrnd_rcm[z_indx] = lev_sum / float(count)
            else:
                # Set the relevant "output" variables to a value just to have them set.
                # They will be flagged out of interacting with the code.
                self.bkgrnd_rcm = 0
                self.altitude_bkgrnd_rcm = 0
                self.time_bkgrnd_rcm = 0
                self.start_alt_idx = 0
                self.end_alt_idx = 0
        else:
            # Set the relevant "output" variables to a value just to have them set.
            # They will be flagged out of interacting with the code.
            self.bkgrnd_rcm = 0
            self.altitude_bkgrnd_rcm = 0
            self.time_bkgrnd_rcm = 0
            self.start_alt_idx = 0
            self.end_alt_idx = 0

        self.generatePanels()

    def addVariable(self, variable_def_dict):
        """
        Given basic details about a variable, this
        creates lines/panels for the variable and then adds that
        variable to the list of all variables for this VariableGroup

        **Valid dict keys/options include:**

        * *var_names*: A list of names various models refer to this variable as. E.g. ['wprtp', 'WPRTP', 'wpqtp']. This
          parameter can also take in references to calc functions. Names/functions are prioritized left -> right.

        * *[model name]_conv_factor*: (optional) Numeric value to scale a model's variable by. E.g. 1/1000, or 100
          Valid models are: clubb, hoc, r408, sam, cam, wrf, e3sm. E.g. "clubb_conv_factor"

        * *type*: (optional) Override the default type 'profile' with either 'budget' or 'timeseries'

        * *title*: (optional) Override the default panel title, or provide one if it's not specified in the netcdf file.

        * *sci_scale*: (optional) A numerical power of ten to scale the panel axis to. E.g. a value of "4" results in the
          axis being scaled to x 1e4

        * *axis_title*: (optional) Override the default dependent axis title, or provide one if it's not
          specified in the netcdf file.

        * *lines*: (optional) Defines lines to plot for budget cases. Passed separately because
          it's a lot of text. The lines parameter currently does not support defining model-specific names/functions
          This is given in the form of a list of lines, here's an example:

        .. code-block:: python
            :linenos:

            thlm_budget_lines = [
                {'var_names': ['thlm_bt'], 'legend_label': 'thlm_bt'},
                {'var_names': ['thlm_ma'], 'legend_label': 'thlm_ma'},
                {'var_names': ['thlm_ta'], 'legend_label': 'thlm_ta'},
                {'var_names': ['thlm_mc'], 'legend_label': 'thlm_mc'},
                {'var_names': ['thlm_clipping',self.getThlmClipping], 'legend_label': 'thlm_clipping'},
                {'var_names': ['radht'], 'legend_label': 'radht'},
                {'var_names': ['ls_forcing', self.getThlmLsforcing], 'legend_label': 'thlm_ls_forcing'},
                {'var_names': ['thlm_residual',self.getThlmResidual], 'legend_label': 'thlm_residual'},
            ]

        :param variable_def_dict: A dict containing the information defining a variable. These definitions are declared
            inside a VariableGroup child class, e.g. VariableGroupBase.

        :return: None
        """
        var_names = variable_def_dict['var_names']

        plot_sam_benchmark = self.sam_benchmark_dataset is not None and len(self.sam_benchmark_dataset) > 0 \
                             and len(var_names['sam']) > 0
        plot_coamps_benchmark = self.coamps_benchmark_dataset is not None and len(self.coamps_benchmark_dataset) > 0 \
                      and len(var_names['coamps']) > 0
        plot_wrf_benchmark = self.wrf_benchmark_dataset is not None and len(self.wrf_benchmark_dataset) > 0 \
                      and len(var_names['wrf']) > 0
        plot_r408 = self.r408_datasets is not None and len(self.r408_datasets) > 0 \
                    and len(var_names['r408']) > 0
        plot_hoc = self.hoc_datasets is not None and len(self.hoc_datasets) > 0 \
                   and len(var_names['hoc']) > 0

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
        if plot_sam_benchmark:
            all_lines.extend(self.__getVarLinesForModel__('sam', variable_def_dict, self.sam_benchmark_dataset))
        if plot_coamps_benchmark:
            all_lines.extend(self.__getVarLinesForModel__('coamps', variable_def_dict, self.coamps_benchmark_dataset))
        if plot_wrf_benchmark:
            all_lines.extend(self.__getVarLinesForModel__('wrf', variable_def_dict, self.wrf_benchmark_dataset))
        if plot_r408:
            all_lines.extend(self.__getVarLinesForModel__('r408', variable_def_dict, self.r408_datasets))
        if plot_hoc:
            all_lines.extend(self.__getVarLinesForModel__('hoc', variable_def_dict, self.hoc_datasets))

        # Plot input folders
        if plot_clubb:
            for input_folder in self.clubb_datasets:
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
        conversion_factors = self.__getConversionFactors__(variable_def_dict)

        lines = None
        if 'lines' in variable_def_dict.keys():
            lines = variable_def_dict['lines']

        panel_type = self.default_panel_type
        if self.time_height:
            panel_type = Panel.TYPE_TIMEHEIGHT
        elif 'type' in variable_def_dict.keys():
            panel_type = variable_def_dict['type']

        if not isinstance(dataset, dict):
            dataset = {model_name: dataset}

        line_style = ""
        if label == "no label found":
            label = Style_definitions.BENCHMARK_LABELS[model_name]
            line_style = Style_definitions.BENCHMARK_LINE_STYLES[model_name]

        all_lines = []

        data_was_calculated = False

        # Process extra variable lines (e.g. budget lines)
        # Time-height plots can only display one variable each.
        # So a line parameter is incompatible.
        if lines is not None:
            if panel_type == Panel.TYPE_TIMEHEIGHT:
                logToFile('Warning. Panel type is time-height but a lines argument was found in for variable ' +
                     str(all_model_var_names[model_name]) + '. Those are not compatible. Ignoring lines.')
            else:
                for line in lines:
                    line['calculated'] = False
                    if (model_name + '_calc') in line.keys() and \
                            not self.__varnamesInDataset__(line['var_names'], dataset):
                        plot_data, z = line[(model_name + '_calc')](dataset_override=dataset)

                        #kludgy trimming for these variables since they are processed here but never trimmed
                        if np.any(z<self.height_min_value):
                            z_min_index=np.max(np.argwhere(z<self.height_min_value))+1
                        else:
                            z_min_index=0
                        if np.any(z>self.height_max_value):
                            z_max_index=np.min(np.argwhere(z>self.height_max_value))-1
                        else:
                            z_max_index=len(z)-1
                        plot_data=plot_data[z_min_index:z_max_index]
                        z=z[z_min_index:z_max_index]
                        #end of kludgy trimming

                        if self.animation is None:
                            plot = Line(plot_data, z, line_format=line_style, label=line['legend_label'])
                        else:
                            plot = Contour(x_data=z['time'],y_data=z['height'],c_data=plot_data,
                                           colors=Style_definitions.CONTOUR_CMAP_GENERAL,
                                           label=line['legend_label'], line_format=line_style)
                        all_lines.append(plot)
                        line['calculated'] = True

        if len(all_model_var_names[model_name]) > 0 and not data_was_calculated:
            all_lines.extend(self.__getVarLines__(all_model_var_names[model_name], dataset,
                                                  conversion_factor=conversion_factors[model_name],
                                                  label=label,
                                                  line_format=line_style,
                                                  override_panel_type=panel_type,
                                                  lines=lines, model_name=model_name))
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
            for dict_of_datasets in self.e3sm_datasets.values():
                datasets.extend(dict_of_datasets.values())
        if self.cam_datasets is not None:
            for dict_of_datasets in self.cam_datasets.values():
                datasets.extend(dict_of_datasets.values())
        if self.sam_datasets is not None:
            for dict_of_datasets in self.sam_datasets.values():
                datasets.extend(dict_of_datasets.values())
        if self.wrf_datasets is not None:
            for i in self.wrf_datasets.values():
                if isinstance(i, dict):
                    datasets.extend(i.values())
                elif isinstance(i, Dataset):
                    datasets.append(i)
                else:
                    raise TypeError("getTextDefiningDataset recieved an unexpected format of datasets")
        if self.sam_benchmark_dataset is not None:
            datasets.extend(self.sam_benchmark_dataset.values())
        if self.coamps_benchmark_dataset is not None:
            datasets.extend(self.coamps_benchmark_dataset.values())
        if self.wrf_benchmark_dataset is not None:
            datasets.extend(self.wrf_benchmark_dataset.values())
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
            # Only display a variable if we are either not doing a priority run or the priority flag has been set and is True
            if not self.priority_vars or (self.priority_vars and 'priority' in variable.keys() and variable['priority']):
                title = variable['title']
                axis_label = variable['axis_title']
                plotset = variable['plots']
                centered = False
                panel_type = self.default_panel_type
                if self.time_height:
                    panel_type = Panel.TYPE_TIMEHEIGHT
                elif 'type' in variable.keys():
                    panel_type = variable['type']
                if 'sci_scale' in variable.keys():
                    sci_scale = variable['sci_scale']
                else:
                    sci_scale = None
                if 'centered' in variable.keys():
                    centered = variable['centered']

                if panel_type == Panel.TYPE_TIMEHEIGHT:
                    if 'var_min' in variable and 'var_max' in variable:
                        panel = ContourPanel(plotset, title=title, dependent_title=axis_label, panel_type=panel_type, var_min=variable['var_min'], var_max=variable['var_max'], file_identifier=variable['file_identifier'])
                    else:
                        panel = ContourPanel(plotset, title=title, dependent_title=axis_label, panel_type=panel_type)
                elif self.animation is not None:
                    panel = AnimationPanel(plotset, self.bkgrnd_rcm, self.altitude_bkgrnd_rcm, self.time_bkgrnd_rcm,
                                           self.start_alt_idx, self.end_alt_idx, panel_type=panel_type, title=title,
                                           dependent_title=axis_label, sci_scale=sci_scale, centered=centered,
                                           bkgrnd_rcm_flag=self.bkgrnd_rcm_flag)
                else:
                    if 'file_identifier' in variable:
                        panel = Panel(plotset, self.bkgrnd_rcm, self.altitude_bkgrnd_rcm, self.start_alt_idx, self.end_alt_idx,
                                      panel_type=panel_type, title=title, dependent_title=axis_label, sci_scale=sci_scale,
                                      centered=centered, bkgrnd_rcm_flag=self.bkgrnd_rcm_flag, file_identifier=variable['file_identifier'])
                    else:
                        panel = Panel(plotset, self.bkgrnd_rcm, self.altitude_bkgrnd_rcm, self.start_alt_idx, self.end_alt_idx,
                                      panel_type=panel_type, title=title, dependent_title=axis_label, sci_scale=sci_scale,
                                      centered=centered, bkgrnd_rcm_flag=self.bkgrnd_rcm_flag)
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
        all_lines = variable_def_dict['plots']
        panel_type = self.default_panel_type
        if self.time_height:
            panel_type = Panel.TYPE_TIMEHEIGHT
        elif 'type' in variable_def_dict.keys():
            panel_type = variable_def_dict['type']
        try:
            first_input_datasets = self.getTextDefiningDataset()
        except FileExistsError as e:
            logToFile(str(e))
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
                all_var_names = []
                for model_var_names in var_names.values():
                    all_var_names.extend(model_var_names)
                # Get axis title for any of the var_names given in variable_def_dict
                imported_axis_title = data_reader.getAxisTitle(first_input_datasets, all_var_names)
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
        if isinstance(datasets, list):
            for name in varnames:
                for dataset in datasets:
                    if name in dataset.variables.keys():
                        return True
        elif isinstance(datasets, dict):
            for name in varnames:
                for dataset in datasets.values():
                    if name in dataset.variables.keys():
                        return True
        elif isinstance(datasets, Dataset):
            for name in varnames:
                if name in datasets.variables.keys():
                    return True
        else:
            raise ValueError("Unknown data type for dataset")

        return False

    def __getRelativeVarName__(self, var_names):
        """
        Returns the varname for a variable relative to the models that were being plotted

        :param var_names: The dict of various names for the given variable. More info found in the addVariable() docstring

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

    def __getConversionFactors__(self, variable_def_dict):
        """
        This is a helper method that loads parameters from the variables
        definition and assigns them to python variables accordingly.

        :param variable_def_dict: Definition of the variable being used. More info found in the addVariable() docstring
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
                        lines=None, conversion_factor=1, model_name="unknown"):
        """
        Get a list of Line objects for a specific clubb variable. If sam_benchmark_dataset is specified it will also
        attempt to generate Lines for the SAM equivalent variables, using the name conversions found in
        VarnameConversions.py. If a SAM variable needs to be calculated (uses an equation) then it will have
        to be created within that variable group's dataset and not here.

        :param var_names: A string or list of strings. More info found in the addVariable() docstring
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
        :param model_name: The name of the model being plotted.
        :return: A list of Line objects containing model dependent_data from whatever models were passed in.
            Returns None if requested variable is not found.
        """
        panel_type = self.default_panel_type
        if override_panel_type is not None:
            panel_type = override_panel_type
        all_lines = []

        if isinstance(ncdf_datasets, Dataset):
            ncdf_datasets = {'converted_to_dict': ncdf_datasets}

        if self.animation is not None:
            if panel_type is Panel.TYPE_PROFILE:
                contours = self.__getProfileLineForAnim__(var_names, ncdf_datasets, label, line_format, conversion_factor,
                                                          lines=lines,model_name=model_name)
                if contours is not None:
                    all_lines.extend(contours)
            elif panel_type is Panel.TYPE_SUBCOLUMN:
                lines = self.__getSubcolumnLinesForAnim__(var_names, ncdf_datasets, label, line_format, conversion_factor,
                                                    avg_axis, lines=lines, model_name=model_name)
                all_lines.extend(lines)
            elif panel_type is Panel.TYPE_BUDGET:
                budget_lines = self.__getBudgetLines__(lines, ncdf_datasets, model_name=model_name)
                all_lines.extend(budget_lines)
        elif panel_type is Panel.TYPE_PROFILE:
            profile_lines = self.__getProfileLine__(var_names, ncdf_datasets, label, line_format, conversion_factor,
                                                    avg_axis, lines=lines, model_name=model_name)
            if profile_lines is not None:
                all_lines.extend(profile_lines)
        elif panel_type is Panel.TYPE_BUDGET:
            budget_lines = self.__getBudgetLines__(lines, ncdf_datasets, model_name=model_name)
            all_lines.extend(budget_lines)
        elif panel_type is Panel.TYPE_TIMESERIES:
            line = self.__getTimeseriesLine__(var_names, ncdf_datasets, self.end_time, conversion_factor, label,
                                              line_format)
            if line is not None:
                all_lines.append(line)
        elif panel_type is Panel.TYPE_TIMEHEIGHT:
            contour = self.__getTimeHeightContours__(var_names, ncdf_datasets, label, conversion_factor)
            if contour is not None:
                all_lines.append(contour)
        elif panel_type is Panel.TYPE_SUBCOLUMN:
            lines = self.__getSubcolumnLines__(var_names, ncdf_datasets, label, line_format, conversion_factor,
                                                    avg_axis, lines=lines, model_name=model_name)
            all_lines.extend(lines)
        elif panel_type:
            raise ValueError('Invalid panel type ' + panel_type + '. Valid options are profile, budget, timeseries')

        return all_lines

    def __getProfileLine__(self, varnames, dataset, label, line_format, conversion_factor, avg_axis, lines=None,
                           model_name="unknown"):
        """
        Assumes variable can be plotted as a profile and returns a Line object
        representing the given variable for a profile plot. Profile plots are plots that show Height vs Varname

        :param varnames: A list of variable names
        :param dataset: NetCdf4 Dataset object containing the model output being plotted
        :param label: This is the name that will be shown on the legend for this line
        :param line_format: A string representing how this line should be formatted using pyplot formatting.
            Recommended value: emtpy string ""
        :param conversion_factor: This is a numerical value that will be multiplied element-wise to the variable.
            It's useful for doing basic model to model conversions, e.g. SAM -> CLUBB.
        :param avg_axis: Values will be averaged along this axis. This is basically only used for time-averaging
            profile plots.
        :param lines: A lines parameter definition as found in a VariableGroup___.py. See addVariable() for more details.
        :return: Line object representing the given variable for a profile plot
        """
        output_lines = []
        variable = NetCdfVariable(varnames, dataset, independent_var_names=Case_definitions.HEIGHT_VAR_NAMES,
                                  start_time=self.start_time, end_time=self.end_time, min_height=self.height_min_value,
                                  max_height=self.height_max_value, avg_axis=avg_axis,
                                  conversion_factor=conversion_factor, model_name=model_name)
        variable.trimArray(self.height_min_value, self.height_max_value, data=variable.independent_data)

        if lines is not None:
            additional_lines = self.__processLinesParameter__(lines, dataset, line_format=line_format,
                                                              label_suffix=label, model_name=model_name)
            output_lines.extend(additional_lines)
        else:
            line = Line(variable, variable.independent_data, label=label, line_format=line_format)
            output_lines.append(line)
        return output_lines

    def __getProfileLineForAnim__(self, varnames, dataset, label, line_format, conversion_factor, lines=None,
                                 model_name="unknown"):
        """
        Performs the same function as __getProfileLine__() except keeps 2D data for animations.
        """
        output_contours=[]
        variable = NetCdfVariable(varnames, dataset,
                                      independent_var_names={'time': Case_definitions.TIME_VAR_NAMES,
                                                             'height': Case_definitions.HEIGHT_VAR_NAMES},
                                  start_time=self.start_time, end_time=self.end_time, avg_axis=2,
                                  conversion_factor=conversion_factor)
        if variable.dependent_data.ndim != 2:
            logToFile("Warning: Variable {} can not be plotted as an animation ".format(varnames) +
                 "array's dimension is {} and not 2. Returning no data.".format(variable.dependent_data.ndim))
            return None
        # From here on: variable.dependent_data.ndim == 2
        variable.trimArray(self.start_time, self.end_time, data=variable.independent_data['time'], axis=0)
        # Do we want to trim height?
        variable.trimArray(self.height_min_value, self.height_max_value, data=variable.independent_data['height'],
                           axis=1)

        if lines is not None:
            contour = self.__processLinesParamForAnim__(lines, dataset, line_format=line_format,
                                                              label_suffix=label, model_name=model_name)
            output_contours.extend(contour)
        else:
            contour = Contour(x_data=variable.independent_data['time'], y_data=variable.independent_data['height'],
                              c_data=variable, colors=Style_definitions.CONTOUR_CMAP_GENERAL, label=label,line_format=line_format)
            output_contours.append(contour)

        return output_contours

    def __getSubcolumnLines__(self, varnames, dataset, label, line_format, conversion_factor, avg_axis, lines=None,
                           model_name="unknown"):
        """
        Generate Line objects for subcolumn panels.

        :param varnames: A list of variable names
        :param dataset: NetCdf4 Dataset object containing the model output being plotted
        :param label: This is the name that will be shown on the legend for this line
        :param line_format: A string representing how this line should be formatted using pyplot formatting.
            Recommended value: emtpy string ""
        :param conversion_factor: This is a numerical value that will be multiplied element-wise to the variable.
            It's useful for doing basic model to model conversions, e.g. SAM -> CLUBB.
        :param avg_axis: Values will be averaged along this axis. This is basically only used for time-averaging
            profile plots.
        :param lines: A lines parameter definition as found in a VariableGroup___.py. See addVariable() for more details.
        :return: Line objects representing the given variable for a subcolumn plot
        """
        output_lines = []
        variable = NetCdfVariable(varnames, dataset["subcolumns"], independent_var_names=Case_definitions.HEIGHT_VAR_NAMES,
                                  start_time=self.start_time, end_time=self.end_time, min_height=self.height_min_value,
                                  max_height=self.height_max_value, avg_axis=avg_axis,
                                  conversion_factor=conversion_factor, model_name=model_name)
        variable.trimArray(self.height_min_value, self.height_max_value, data=variable.independent_data, axis=0)
        # this if statement accommodates cases with some but not all subcolumn vars (e.g. RICO_SILHS)
        if len(variable.dependent_data.shape) == 2:
            for i in range(len(variable.dependent_data[0])):
                label_with_offset = label + "_" + str(i+1)
                x_data_i = variable.dependent_data[:,i]
                line = Line(x_data_i, variable.independent_data, label=label_with_offset)
                output_lines.append(line)

        if lines is not None:
            additional_lines = self.__processLinesParameter__(lines, dataset, line_format=line_format,
                                                              label_suffix=label, model_name=model_name)
            output_lines.extend(additional_lines)

        return output_lines

    def __getSubcolumnLinesForAnim__(self, varnames, dataset, label, line_format, conversion_factor, avg_axis, lines=None,
                           model_name="unknown"):
        """
        Same as __getSubcolumnLines__() but retains 2D data for animations.
        """
        output_lines = []
        variable = NetCdfVariable(varnames, dataset["subcolumns"],
                                  independent_var_names={'time': Case_definitions.TIME_VAR_NAMES,
                                                         'height': Case_definitions.HEIGHT_VAR_NAMES},
                                  start_time=self.start_time, end_time=self.end_time, avg_axis=2,
                                  conversion_factor=conversion_factor, model_name=model_name)
        variable.trimArray(self.start_time, self.end_time, data=variable.independent_data['time'], axis=0)
        variable.trimArray(self.height_min_value, self.height_max_value, data=variable.independent_data['height'],axis=1)
        # this if statement accommodates cases with some but not all subcolumn vars (e.g. RICO_SILHS)
        if len(variable.dependent_data.shape) == 3:
            for i in range(len(variable.dependent_data[0,0])):
                label_with_offset = label + "_" + str(i+1)
                data_i = variable.dependent_data[:,:,i]
                contour = Contour(x_data=variable.independent_data['time'], y_data=variable.independent_data['height'],
                                  c_data=data_i, colors=Style_definitions.CONTOUR_CMAP_GENERAL, label=label_with_offset,
                                  line_format=line_format)
                output_lines.append(contour)

        if lines is not None:
            additional_lines = self.__processLinesParamForAnim__(lines, dataset, line_format=line_format,
                                                              label_suffix=label, model_name=model_name)
            output_lines.extend(additional_lines)

        return output_lines

    def __getTimeseriesLine__(self, varname, dataset, end_time, conversion_factor, label, line_format):
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
        variable.trimArray(0, variable.end_time, data=variable.independent_data)
        line = Line(variable.independent_data, variable, label=label, line_format=line_format)
        return line

    def __getBudgetLines__(self, lines, dataset, model_name="unknown"):
        """
        Returns a list of Line objects (or 2D objects in case of animations) for a budget plot for each
        variable defined in lines

        :param lines: A list of line definitions containing the variable names, legend labels etc.
            See the addVariable documentation above for a complete description
        :param dataset: NetCdf4 Dataset object containing the model output being plotted
        :param model_name: A string name of the model being plotted. I.e. clubb, hoc, r408, sam, cam, wrf, e3sm
        :return: A list of Line objects for a budget plot derived from lines
        """
        if self.animation is not None:
            output_lines = self.__processLinesParamForAnim__(lines, dataset, model_name=model_name)
        else:
            output_lines = self.__processLinesParameter__(lines, dataset, model_name=model_name)
        return output_lines

    def __getTimeHeightContours__(self, varname, dataset, label, conversion_factor):
        """
        Return a Contour object for a time-height plot,
        plotting a 2-dimensional contour plot with time as x axis and height as y axis.

        :param varname: name of the variable being plotted
        :param dataset: dataset containing the variable
        :param label: The legend label for this variable
        :param conversion_factor: A conversion factor to apply to the variable (usually 1.0)
        :return: A Contour object that can be used to create a time-height plot.
        """
        # This does not yet work since we need two different independent_data arrays
        variable = NetCdfVariable(varname, dataset,
                                  independent_var_names={'time': Case_definitions.TIME_VAR_NAMES,
                                                         'height': Case_definitions.HEIGHT_VAR_NAMES},
                                  start_time=self.start_time, end_time=self.end_time, avg_axis=2,
                                  conversion_factor=conversion_factor)
        if variable.dependent_data.ndim == 3:
            variable.dependent_data = variable.dependent_data[:,:,0]
            logToFile("Warning: assuming {} is a SILHS subcolumn variable. Plotting only the ".format(varname) +
                 "first subcolumn for time-height plots.  If {} is not a subcolumn variable, ".format(varname) +
                 "please update src/VariableGroup.py.")
        elif variable.dependent_data.ndim == 1 or variable.dependent_data.ndim > 3:
            logToFile("Warning: Variable {} can not be plotted as time-height plot as the data ".format(varname) +
                 "array's dimension is {} and not 2. Returning no data.".format(variable.dependent_data.ndim))
            return None
        # From here on: variable.dependent_data.ndim == 2
        variable.trimArray(self.start_time, self.end_time, data=variable.independent_data['time'], axis=0)
        # Do we want to trim height?
        variable.trimArray(self.height_min_value, self.height_max_value, data=variable.independent_data['height'],
                           axis=1)

        # Create Contour instance to return
        contour = Contour(x_data=variable.independent_data['time'], y_data=variable.independent_data['height'],
                          c_data=variable, colors=Style_definitions.CONTOUR_CMAP_GENERAL, label=label)
        return contour

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
        if self.time_height or self.animation is not None:
            var_ncdf = NetCdfVariable(varname, datasets,
                                      independent_var_names={'time': Case_definitions.TIME_VAR_NAMES,
                                                             'height': Case_definitions.HEIGHT_VAR_NAMES},
                                      conversion_factor=conversion_factor, start_time=self.start_time,
                                      end_time=self.end_time, avg_axis=2)
            var_ncdf.trimArray(self.start_time, self.end_time, data=var_ncdf.independent_data['time'], axis=0)
            # Do we want to trim height?
            var_ncdf.trimArray(self.height_min_value, self.height_max_value, data=var_ncdf.independent_data['height'],
                               axis=1)
        else:
            var_ncdf = NetCdfVariable(varname, datasets, independent_var_names=Case_definitions.HEIGHT_VAR_NAMES,
                                      conversion_factor=conversion_factor, start_time=self.start_time,
                                      end_time=self.end_time)
            # I think trimming height is redundant with later processes, and this line can cause a problem with
            # data going outside start/end limits. Commenting for now, maybe eventually delete. BAS 11/20
            #var_ncdf.trimArray(self.height_min_value, self.height_max_value, data=var_ncdf.independent_data)
        var_data = var_ncdf.dependent_data
        indep_data = var_ncdf.independent_data
        return var_data, indep_data, var_ncdf.ncdf_data  # changed datasets to var_ncdf.ncdf_data to fix budget plots

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

    def pickNonZeroOutput(self, output1, output2, favor_output=None):
        """
        Sometimes there are more than 1 ways to calculate a variable, and
        it is impossible to know which equation to use before hand. In these
        cases, all equations should be used and their output given to this
        function. This function will then attempt to determine which of the
        outputs from the equations is most likely to be the expected answer.

        This function currently works by determining picking the non-nan and non-zero array.
        In the event that both arrays contain non-nan and non-zero output, it will return output1.
        In the event both arrays are identical w/in tolerance, it will return output1.

        If both arrays are zero or nan, then it will return the first array arbitrarily.

        Example usage:

        .. code-block:: python
            :linenos:
#
#            def get_rc_coef_zm_X_wprcp_coamps_calc(self, dataset_override=None):
#                wprlp, z, dataset = self.getVarForCalculations(['thlpqcp', 'wpqcp', 'wprlp'], dataset)
#                ex0, z, dataset = self.getVarForCalculations(['ex0'], dataset)
#                p, z, dataset = self.getVarForCalculations('p', dataset)
#                thvm, z, dataset = self.getVarForCalculations('thvm', dataset)
#                output1 = wprlp * (2.5e6 / (1004.67 * ex0) - 1.61 * thvm)
#                output2 = wprlp * (2.5e6 / (1004.67 * ((p / 1.0e5) ** (287.04 / 1004.67))) - 1.61 * thvm)
#                output = self.pickMostLikelyOutputList(output1, output2)
#                return output, z
#
#
        :param output1: An arraylike list of numbers outputted from a variable calculation
            (e.g.     output1 = wprlp * (2.5e6 / (1004.67 * ex0) - 1.61 * thvm)
        :param output2: An arraylike list of numbers outputted from a variable calculation
            (e.g.     output2 = wprlp * (2.5e6 / (1004.67 * ((p / 1.0e5) ** (287.04 / 1004.67))) - 1.61 * thvm)
        :param favor_output: (optional) If an array cannot easily be picked by looking at whether they are zeros/NaNs,
            use this return the data from this array. I.e. if both output1 and output2 contain unique data, the data in
            favor_output will be returned. This array should be identical to either output1 or output2.
        :return: The listlike array of the two datasets most likely to be the correct answer.
        """
        out1_is_nan = math.isnan(mean(np.ndarray.flatten(output1)))
        out2_is_nan = math.isnan(mean(np.ndarray.flatten(output2)))

        out1_is_zero = np.all(np.isclose(output1, 0))
        out2_is_zero = np.all(np.isclose(output2, 0))

        return_out_1 = (out2_is_zero or out2_is_nan) and not out1_is_nan
        return_out_2 = (out1_is_zero or out1_is_nan) and not out2_is_nan

        output1_close_to_output2 = np.allclose(output1, output2)
        return_any = (out1_is_zero and out2_is_zero) or (out1_is_nan and out2_is_nan) or output1_close_to_output2

        if return_out_2:
            return output2
        elif return_out_1:
            return output1
        elif return_any:
            return output1
        elif favor_output is not None:
            return favor_output
        raise UserWarning("Failed to find an easy answer to the best output.")

    def __processLinesParameter__(self, lines, dataset, label_suffix="", line_format="", model_name="unknown"):
        """
        This method processes a 'lines' parameter from a variable definition and translates
        it into a list of Line objects for plotting.

        :param lines: a lines definition. See VariableGroupBaseBudgets.py for examples
        :param dataset: A netcdf Dataset object or a dict of datasets
        :param label_suffix: String. Optional suffix to add to the end of legend label.
        :param line_format: matplotlib line formating string. A blank "" string results in auto formatting
        :param model_name: String name of the model. I.e.  clubb, hoc, r408, sam, cam, wrf, e3sm
        :return: list of Line objects for plotting
        """
        if isinstance(lines, dict):
            lines = lines[model_name]
        output_lines = []
        for line_definition in lines:
            if line_definition['calculated'] is True:
                continue
            varnames = line_definition['var_names']
            label = line_definition['legend_label']
            if label_suffix != "":
                #label = line_definition['legend_label'] + " " + label_suffix
                label = label
            variable = NetCdfVariable(varnames, dataset, independent_var_names=Case_definitions.HEIGHT_VAR_NAMES,
                                      start_time=self.start_time, end_time=self.end_time)
            variable.trimArray(self.height_min_value, self.height_max_value, data=variable.independent_data)
            line_definition = Line(variable, variable.independent_data, label=label,
                                   line_format=line_format)  # uses auto-generating line format
            output_lines.append(line_definition)
        return output_lines

    def __processLinesParamForAnim__(self, lines, dataset, label_suffix="", line_format="", model_name="unknown"):
        """
        This method is the same as __processLinesParameter__() except it keeps 2D data for animations.

        :params: same as __processLinesParameter__()
        :return: list of Contour objects for plotting
        """
        if isinstance(lines, dict):
            lines = lines[model_name]
        output_lines = []
        for line_definition in lines:
            if line_definition['calculated'] is True:
                continue
            varnames = line_definition['var_names']
            label = line_definition['legend_label']
            if label_suffix != "":
                label = line_definition['legend_label'] + " " + label_suffix
            variable = NetCdfVariable(varnames, dataset,
                                  independent_var_names={'time': Case_definitions.TIME_VAR_NAMES,
                                                         'height': Case_definitions.HEIGHT_VAR_NAMES},
                                  start_time=self.start_time, end_time=self.end_time, avg_axis=2)
            if len(variable.dependent_data.shape) > 1:
                variable.trimArray(self.start_time, self.end_time, data=variable.independent_data['time'], axis=0)
                variable.trimArray(self.height_min_value, self.height_max_value, data=variable.independent_data['height'],axis=1)
                contour = Contour(x_data=variable.independent_data['time'], y_data=variable.independent_data['height'],
                              c_data=variable, colors=Style_definitions.CONTOUR_CMAP_GENERAL, label=label,line_format=line_format)
                output_lines.append(contour)
        return output_lines