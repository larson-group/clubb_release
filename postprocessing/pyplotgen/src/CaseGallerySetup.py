"""
:author: Nicolas Strike
:date: 2019
"""
import os
from warnings import warn

import numpy as np

from config.VariableGroupBaseBudgets import VariableGroupBaseBudgets
from config.VariableGroupSamBudgets import VariableGroupSamBudgets
from config.VariableGroupSubcolumns import VariableGroupSubcolumns
from src.DataReader import DataReader
from src.Panel import Panel


class CaseGallerySetup:
    """
    A case is basically a collection of variable groups and
    panels that share attributes like start/end time and get plotted together.
    This is a generic class used by Case_Definitions.py to provide functionality.
    In order to create a new case, please add the case's definition to Case_Definitions.py (don't forget to add the
    definition to the list CASES_TO_PLOT = [...] at the bottom of the Case_Definitions.py file).

    For information on the input parameters of this class, please see the documentation for the
    ``__init__()`` method.
    """

    def __init__(self, case_definition, clubb_folders=[], diff_datasets=None, sam_folders=[""], wrf_folders=[""],
                 plot_les=False, plot_budgets=False, plot_r408=False, plot_hoc=False, e3sm_dirs=[], cam_folders=[],
                 time_height=False, animation=None, plot_subcolumns=False):
        """
        Initialize a CaseGallerySetup object with the passed parameters
        :param case_definition: dict containing case specific elements. These are pulled in from Case_definitions.py,
            see Case_definitions.py for details on how to structure the dict
        :param clubb_folders: dict containing Dataset objects holding the dependent_data needed for the case.
            The key for each value/Dataset in the dict is set to the ext provided in the filename (e.g. sfc, zt, zm)
        :param diff_datasets: Takes in netcdf datasets. If datasets are passed in, pyplotgen will plot the numeric
            difference between the folder passed in an the clubb folder.
        :param sam_folders: List of foldernames containing sam netcdf files to be plotted
        :param wrf_folders: List of foldernames containing wrf netcdf files to be plotted
        :param plot_les: If True pyplotgen plots LES lines, if False pyplotgen does not plot LES lines
        :param plot_budgets: If True pyplotgen will plot Budgets in addition to the other plots
            If False (default), pyplotgen will not plot budgets
        :param plot_r408: If True, pyplotgen will plot the Chris Golaz 'best ever' clubb r408 dependent_data lines
            If False (default), pyplotgen will not plot the Chris Golaz 'best ever' clubb r408 dependent_data lines
        :param plot_hoc: If True, pyplotgen will plot the HOC 2005 dependent_data lines
            If False (default), pyplotgen will not plot the HOC 2005 dependent_data lines
        :param e3sm_dirs: List of foldernames containing e3sm netcdf files to be plotted
        :param cam_folders: List of foldernames containing cam netcdf files to be plotted
        :param time_height: TODO
        :param animation: TODO
        """
        self.name = case_definition['name']
        self.start_time = case_definition['start_time']
        self.end_time = case_definition['end_time']
        self.var_groups = case_definition['var_groups']
        self.height_min_value = case_definition['height_min_value']
        self.height_max_value = case_definition['height_max_value']
        self.clubb_datasets = clubb_folders
        self.blacklisted_variables = case_definition['blacklisted_vars']
        self.plot_budgets = plot_budgets
        self.plot_r408 = plot_r408
        self.plot_hoc = plot_hoc
        self.e3sm_folders = e3sm_dirs
        self.sam_folders = sam_folders
        self.cam_folders = cam_folders
        self.wrf_folders = wrf_folders
        self.diff_datasets = diff_datasets
        self.next_panel_alphabetic_id_code = 97
        self.time_height = time_height
        self.animation = animation
        self.panels = []
        self.diff_panels = []
        self.plot_subcolumns = plot_subcolumns

        self.VALID_MODEL_NAMES = ['clubb', 'clubb_hoc','clubb_r408', 'e3sm', 'sam', 'cam', 'wrf', 'coamps']

        if 'disable_budgets' in case_definition.keys() and case_definition['disable_budgets'] is True:
            self.plot_budgets = False

        # Load benchmark files
        self.sam_benchmark_file = self.__loadModelFiles__(None,case_definition,"sam")
        self.coamps_datasets = self.__loadModelFiles__(None, case_definition, "coamps")
        self.r408_datasets = self.__loadModelFiles__(None, case_definition, "clubb_r408")
        self.hoc_datasets = self.__loadModelFiles__(None, case_definition, "clubb_hoc")

        # Load datasets imported via command line parameters
        self.clubb_datasets = self.__loadModelFiles__(clubb_folders, case_definition, "clubb")
        self.sam_datasets = self.__loadModelFiles__(sam_folders, case_definition, "sam")
        self.wrf_datasets = self.__loadModelFiles__(wrf_folders, case_definition, "wrf")
        self.e3sm_file = self.__loadModelFiles__(e3sm_dirs, case_definition, "e3sm")
        self.cam_file = self.__loadModelFiles__(cam_folders, case_definition, "cam")


        self.__generateSubcolumnPanels__()
        self.__generateVariableGroupPanels__()
        self.__generateDiffPanels__()
        self.__generateBudgetPanels__()

    def __generateSubcolumnPanels__(self):
        """
        This function creates the subcolumn panels and adds them into self.panels.
        This function takes no parameters and will only work if both self.plot_subcolumns is True, and
        the model, case, and input folder contain defined data for subcolumn output. Otherwise it will do nothing.

        :return: None. Operates in-place
        """
        # Only attempt subcolumns if enabled and the case defines an output file
        if self.plot_subcolumns and self.clubb_datasets is not None and len(self.clubb_datasets) != 0:
            # for folders_datasets in self.clubb_datasets.values():
            for input_folder in self.clubb_datasets:
                if "subcolumns" in self.clubb_datasets[input_folder].keys():
                    folder_name = os.path.basename(input_folder)
                    subcols_defined_for_this_folder = "subcolumns" in self.clubb_datasets[input_folder]
                    if input_folder in self.clubb_datasets.keys() and subcols_defined_for_this_folder:
                        subcolumn_variables = VariableGroupSubcolumns(self,
                                                    clubb_datasets={folder_name:self.clubb_datasets[input_folder]})
                        self.panels.extend(subcolumn_variables.panels)
                    else:
                        warn("" + folder_name + " does not seem to contain data for case" + self.name)

    def __generateBudgetPanels__(self):
        """
        This function creates the budget panels and adds them into self.panels.
        This function takes no parameters and will only work if self.plot_budgets is True, otherwise it will do nothing

        :return: None. This function operates in-place
        """
        if self.plot_budgets:
            # Create an instance of a budgets VariableGroup. By default, this is VariableGroupBaseBudgets,
            # but for SAM data, use VariableGroupSamBudgets
            if self.clubb_datasets is not None and len(self.clubb_datasets) != 0:
                # for folders_datasets in self.clubb_datasets.values():
                for input_folder in self.clubb_datasets:
                    folder_name = os.path.basename(input_folder)
                    if input_folder in self.clubb_datasets.keys():
                        budget_variables = VariableGroupBaseBudgets(self,
                                                                    clubb_datasets={folder_name:self.clubb_datasets[input_folder]})
                        self.panels.extend(budget_variables.panels)
                    else:
                        warn("" + folder_name + " does not seem to contain data for case" + self.name)
            if self.wrf_datasets is not None and len(self.wrf_datasets) != 0:
                # for folders_datasets in wrf_datasets.values():
                #     budget_variables = VariableGroupBaseBudgets(self, wrf_datasets=folders_datasets)
                for input_folder in self.wrf_datasets:
                    folder_name = os.path.basename(input_folder)
                    budget_variables = VariableGroupBaseBudgets(self,
                                                                wrf_datasets={folder_name:self.wrf_datasets[input_folder]})
                    self.panels.extend(budget_variables.panels)
            if self.e3sm_file is not None and len(self.e3sm_file) != 0:
                for dataset_name in self.e3sm_file:
                    # E3SM dataset must be wrapped in the same form as the clubb datasets
                    e3sm_budgets = VariableGroupBaseBudgets(self,
                                                            e3sm_datasets={dataset_name: self.e3sm_file[dataset_name]})
                    self.panels.extend(e3sm_budgets.panels)
            if self.sam_datasets is not None and len(self.sam_datasets) != 0:
                # for dataset in sam_datasets.values():
                for input_folder in self.sam_datasets:
                    folder_name = os.path.basename(input_folder)
                    budget_variables = VariableGroupSamBudgets(self,
                                                               sam_datasets={folder_name:self.sam_datasets[input_folder]})
                # sam_budgets = VariableGroupSamBudgets(self, sam_datasets=sam_datasets)
                    self.panels.extend(budget_variables.panels)

    def __generateDiffPanels__(self):
        """
        If self.diff_datasets is true (i.e. --diff passed in via command line) then this will generate panels that
        represents the difference of two input folders.

        :return: None. Operates in-place
        """
        # Convert panels to difference panels if user passed in --diff <<folder>>
        if self.diff_datasets is not None:
            # Loop over the VariableGroup classes listed in the 'var_groups' entry
            # for this case in config/Case_definitions.py and create an instance of each of the listed VariableGroups
            for VarGroup in self.var_groups:
                # Call the __init__ function of the VarGroup class and, by doing this, create an instance of it
                diff_group = VarGroup(self, clubb_datasets=self.diff_datasets, sam_file=self.sam_benchmark_file,
                                      coamps_file=self.coamps_datasets, cam_file=self.cam_file, sam_datasets=self.sam_datasets,
                                      r408_file=self.r408_datasets, hoc_dataset=self.hoc_datasets, e3sm_datasets=self.e3sm_file)
                for panel in diff_group.panels:
                    self.diff_panels.append(panel)
            for idx in range(len(self.panels)):
                regular_panel = self.panels[idx]
                diff_panel = self.diff_panels[idx]
                diff_on_y = False
                if regular_panel.panel_type == Panel.TYPE_TIMESERIES:
                    diff_on_y = True
                diff_lines = self.getDiffLinesBetweenPanels(regular_panel, diff_panel, get_y_diff=diff_on_y)
                self.panels[idx].all_plots = diff_lines

    def __generateVariableGroupPanels__(self):
        """
        This function generates the normal profile plots and adds them to self.panels

        :return: None. Operates in-place
        """
        # Loop over the VariableGroup classes listed in the 'var_groups' entry
        for VarGroup in self.var_groups:
            # Calls the __init__ function of the VarGroup class and, by doing this, create an instance of it
            temp_group = VarGroup(self, clubb_datasets=self.clubb_datasets, les_dataset=self.sam_benchmark_file,
                                  coamps_dataset=self.coamps_datasets, sam_datasets=self.sam_datasets,
                                  wrf_datasets=self.wrf_datasets, r408_dataset=self.r408_datasets, hoc_dataset=self.hoc_datasets,
                                  e3sm_datasets=self.e3sm_file, cam_datasets=self.cam_file)
            self.panels.extend(temp_group.panels)

    def __loadModelFiles__(self, folders, case_definition, model_name):
        if model_name not in self.VALID_MODEL_NAMES:
            raise ValueError("Model name " + model_name + " is not a valid model name. Valid model names are: " +
                             str(self.VALID_MODEL_NAMES))
        datareader = DataReader()

        # Load clubb nc files
        model_datasets = {}
        if folders is not None and len(folders) != 0 and case_definition[model_name +'_file'] is not None:
            for foldername in folders:
                files_in_folder = {}
                filenames = case_definition[model_name +'_file']
                for type_ext in filenames:
                    filepath = foldername + filenames[type_ext]
                    ncdf_file = datareader.__loadNcFile__(filepath)
                    if ncdf_file is not None:
                        files_in_folder[type_ext] = ncdf_file
                        model_datasets[foldername] = files_in_folder
        # If is a benchmark
        elif folders is None and case_definition[model_name +'_benchmark_file'] is not None:
            filename_dict = case_definition[model_name +'_benchmark_file']
            for key in filename_dict:
                filename = filename_dict[key]
                ncdf_file = datareader.__loadNcFile__(filename)
                if ncdf_file is not None:
                    model_datasets[key] = ncdf_file

        if len(model_datasets) == 0:
            model_datasets = None
        return model_datasets


    def getDiffLinesBetweenPanels(self, panelA, panelB, get_y_diff=False):
        """
        Given two panels of type Panel, this function calculates the numerical
        difference between each line of the two panels and returns the numerical
        values as a list [] for each line in the panel.

        This function assumes that the dependent_data lines are in the same order on each panel
        and that each panel has the same number of lines

        :param panelA: The first panel for comparison
        :param panelB: The second panel for comparison
        :param get_y_diff: Default behavior is to return the difference along the x-axis,
            setting this to True returns the difference on the y-axis
        :return: A 2D list containing the numerical values for the difference between each line in the given panels.
        """
        linesA = panelA.all_plots
        linesB = panelB.all_plots
        new_lines = linesA
        for idx in range(len(linesA)):
            if get_y_diff:
                a_data = linesA[idx].y
                b_data = linesB[idx].y
                diff_line_data = self.__getArrayDiff__(a_data, b_data)
                new_lines[idx].y = (diff_line_data)
            else:
                a_data = linesA[idx].x
                b_data = linesB[idx].x
                diff_line_data = self.__getArrayDiff__(a_data, b_data)
                new_lines[idx].x = (diff_line_data)
        return new_lines

    def __getArrayDiff__(self, arrA, arrB):
        """
        Returns an array containing the difference between the two arrays.
        arrB is subtracted from arrA and the values are NOT given as absolute value.
        If one array has fewer entries than another, the value of the longer array will be
        appended (i.e. this uses the fill-zeros philosophy).

        :param arrA: The first array (usually contains larger values)
        :param arrB: The second array (usually contains smaller values)
        :return: A numpy array containing |arrA - arrB|
        """
        # Pad shorter array with zeros to match size of longer array and convert to numpy array
        if len(arrA) < len(arrB):
            short = arrA
            long = arrB
        else:
            short = arrB
            long = arrA
        pad_type = type(short[0])
        padding = np.zeros(len(long)-len(short), dtype=pad_type)
        short = np.concatenate((short,padding))

        # Return absolute difference between both arrays
        return np.abs(long-short)

    def plot(self, output_folder, replace_images=False, no_legends=False, thin_lines=False, show_alphabetic_id=False):
        """
        Plot all panels associated with the case, these will be saved to image files in the <<output>>/<<casename>>
        folder

        :param output_folder: Absolute name of the folder to save output into.
        :param replace_images: If True, pyplotgen will overwrite images with the same name.
            If False (default), pyplotgen will add a timestamp to the end of every filename
            (even if there's no filename conflict)
        :param no_legends: If True, pyplotgen will not include a legend on output graphs.
            If False (default), legends will be displayed.
        :param thin_lines: If True, lines plotted will be much thinner than usual.
            If False (default), lines are plotted according to the thickness defined in config/Style_definitions.py
        :param show_alphabetic_id: If True, pyplotgen will add an alphabetic
            label to the top right corner of each plot. These labels will rotate through a-z incrementally.
            If there are more than 26 plots, it will rotate 2 dimensionally,
            e.g. (aa), (ab), (ac),...,(ba),(bb),(bc) and etc.
            If show_alphabetic_id is False (default), this label is not displayed.
            The behavior for rotations greater than zz is not defined,
            and although pyplotgen won't crash it may start to use weird characters.
            The rotation resets between each case,
            e.g. if one case ends on label (ad), the next case will start on (a).
        :return: None
        """
        print("\n\tSaving panels to .png images")
        num_plots = len(self.panels)
        curr_panel_num = 0
        for panel in self.panels:
            print("\r\tplotting ", curr_panel_num, " of ", num_plots, " | ", panel.title)
            if show_alphabetic_id:
                alphabetic_id = self.__getNextAlphabeticID__()
            else:
                alphabetic_id = ""
            plot_paired_lines = True
            if panel.panel_type == panel.TYPE_BUDGET or panel.panel_type == panel.TYPE_SUBCOLUMN:
                plot_paired_lines = False
            panel.plot(output_folder, self.name, replace_images=replace_images, no_legends=no_legends,
                       thin_lines=thin_lines, alphabetic_id=alphabetic_id, paired_plots=plot_paired_lines)
            curr_panel_num += 1
            print("\r\tplotted  ", curr_panel_num, " of ", num_plots, " | ", panel.title, '\n')

    def __getNextAlphabeticID__(self):
        """
        When --show-alphabetic-id is passed in as a run parameter, pyplotgen will add an alphabetic label to each
        plot. These labels are a 1 or 2d rotation through the alphabet (automatically converts to 2d if more than 26
        labels are needed). E.g. this method will return 'a' the first time, 'b' the second time, 'aa' the 27th time,
        'ab' the 28th time, and etc. This function returns the next label as a string, and keeps track of each call
        to this method. Call count tracking is specific per case instance, e.g. if one case plots the label 'bb' last
        the next case will plot 'a' as the first label instead of 'bc'.

        :return: Sequentially next label string for a specific Case object
        """
        a = 97
        z = 122
        num_letters_a_to_z = 26
        if a <= self.next_panel_alphabetic_id_code <= z:
            # Return label with single character
            letter = chr(self.next_panel_alphabetic_id_code)
            # Increment next label code
            self.next_panel_alphabetic_id_code += 1
            return letter
        else:
            # Return label with two characters
            # Find character in first position (Slow rotation)
            first_letter_offset = int((self.next_panel_alphabetic_id_code - a) / num_letters_a_to_z) - 1
            first_letter = chr(a + first_letter_offset)

            # Find character in second position (Fast rotation)
            second_letter_offset = int((self.next_panel_alphabetic_id_code - a * int(
                self.next_panel_alphabetic_id_code / a)) % num_letters_a_to_z)
            second_letter = chr(a + second_letter_offset)

            # Increment next label code
            self.next_panel_alphabetic_id_code += 1
            # Append and return characters
            return first_letter + second_letter
