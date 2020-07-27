"""
:author: Nicolas Strike
:date: 2019
"""
import os
from os import path
from warnings import warn

import numpy as np

from config.VariableGroupBaseBudgets import VariableGroupBaseBudgets
from config.VariableGroupSamBudgets import VariableGroupSamBudgets
from src.DataReader import DataReader
from src.Panel import Panel


class Case:
    """
    A case is basically a collection of variable groups and
    panels that share attributes like start/end time and get plotted together.
    This is a generic class used by Case_Definitions.py to provide functionality.
    In order to create a new case, please add the case's definition to Case_Definitions.py (don't forget to add the
    definition to the list ALL_CASES = [...] at the bottom of the file).
    """

    def __init__(self, case_definition, clubb_folders=[], diff_datasets=None, sam_folders=[""], wrf_folders=[""],
                 plot_les=False, plot_budgets=False, plot_r408=False, plot_hoc=False, e3sm_dirs=[], cam_folders=[],
                 time_height=False, animation=None):
        """
        Initialize a Case object with the passed parameters

        TODO:   - Create function for loading NcFiles to reduce redundant code
                - Bring in line the way CLUBB files are loaded

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
        :param plot_timeheight: TODO
        :param movies: TODO
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

        if 'disable_budgets' in case_definition.keys() and case_definition['disable_budgets'] is True:
            self.plot_budgets = False

        # TODO loading in these files should be refactored

        ## Load nc files for non-clubb models
        # Load single LES dataset
        les_file = None
        if plot_les and case_definition['les_dataset'] is not None:
            datareader = DataReader()
            les_file = datareader.__loadNcFile__(case_definition['les_dataset'])

        # Load SAM nc files from each folder
        sam_datasets = {}
        if sam_folders is not None and len(sam_folders) != 0 and case_definition['sam_file'] is not None:
            datareader = DataReader()
            for foldername in sam_folders:
                sam_filename = foldername + case_definition['sam_file']
                if path.exists(sam_filename):
                    sam_datasets[foldername] = datareader.__loadNcFile__(sam_filename)
                else:
                    warn("Failed to find SAM file " + sam_filename)

        # Load WRF nc files from each folder
        wrf_datasets = {}
        if wrf_folders is not None and len(wrf_folders) != 0 and case_definition['wrf_file'] is not None:
            datareader = DataReader()
            for foldername in wrf_folders:
                files_in_folder = {}
                wrf_filenames = case_definition['wrf_file']
                # Load each different WRF output file from foldername
                for type_ext in wrf_filenames:
                    filepath = foldername + wrf_filenames[type_ext]
                    if path.exists(filepath):
                        files_in_folder[type_ext] = datareader.__loadNcFile__(filepath)
                        wrf_datasets[foldername] = files_in_folder
                    else:
                        warn("Failed to find WRF file " + filepath)

        # Load E3SM nc files from each folder
        e3sm_file = {}
        if e3sm_dirs is not None and len(e3sm_dirs) != 0 and case_definition['e3sm_file'] is not None:
            datareader = DataReader()
            for foldername in e3sm_dirs:
                e3sm_filename = foldername + case_definition['e3sm_file']
                if path.exists(e3sm_filename):
                    e3sm_file[foldername] = datareader.__loadNcFile__(e3sm_filename)
                else:
                    warn("Failed to find E3SM file " + e3sm_filename)

        # Load CAM nc files from each folder
        cam_file = {}
        if cam_folders is not None and len(cam_folders) != 0 and case_definition['cam_file'] is not None:
            datareader = DataReader()
            for foldername in cam_folders:
                cam_filename = foldername + case_definition['cam_file']
                if path.exists(cam_filename):
                    cam_file[foldername] = datareader.__loadNcFile__(cam_filename)
                else:
                    warn("Failed to find CAM file " + cam_filename)

        # Load COAMPS nc files
        coamps_datasets = {}
        if plot_les and case_definition['coamps_dataset'] is not None:
            datareader = DataReader()
            coamps_filenames = case_definition['coamps_dataset']
            # Load the individual COAMPS output files
            for type_ext in coamps_filenames:
                temp_coamps_dataset = datareader.__loadNcFile__(coamps_filenames[type_ext])
                coamps_datasets[type_ext] = temp_coamps_dataset
        else:
            coamps_datasets = None

        # Load clubb nc files
        clubb_datasets = {}
        if clubb_folders is not None and len(clubb_folders) != 0 and case_definition['clubb_file'] is not None:
            datareader = DataReader()
            for foldername in clubb_folders:
                files_in_folder = {}
                clubb_filenames = case_definition['clubb_file']
                # Load each different WRF output file from foldername
                for type_ext in clubb_filenames:
                    filepath = foldername + clubb_filenames[type_ext]
                    if path.exists(filepath):
                        files_in_folder[type_ext] = datareader.__loadNcFile__(filepath)
                        clubb_datasets[foldername] = files_in_folder
                    else:
                        warn("Failed to find CLUBB file " + filepath)

        # Load r408 nc files
        r408_datasets = {}
        if plot_r408 and case_definition['r408_file'] is not None:
            datareader = DataReader()
            r408_filenames = case_definition['r408_file']
            # Load the individual r408 output files
            for type_ext in r408_filenames:
                temp_r408_dataset = datareader.__loadNcFile__(r408_filenames[type_ext])
                r408_datasets[type_ext] = temp_r408_dataset
        else:
            r408_datasets = None

        # Load HOC nc files
        hoc_datasets = {}
        if plot_hoc and case_definition['hoc_file'] is not None:
            datareader = DataReader()
            hoc_filenames = case_definition['hoc_file']
            # Load the individual r408 output files
            for type_ext in hoc_filenames:
                temp_hoc_dataset = datareader.__loadNcFile__(hoc_filenames[type_ext])
                hoc_datasets[type_ext] = temp_hoc_dataset
        else:
            hoc_datasets = None

        self.panels = []
        # Loop over the VariableGroup classes listed in the 'var_groups' entry
        # for this case in config/Case_definitions.py and create an instance of each of the listed VariableGroups
        for VarGroup in self.var_groups:
            # Call the __init__ function of the VarGroup class and, by doing this, create an instance of it
            temp_group = VarGroup(self, clubb_datasets=clubb_datasets, les_dataset=les_file,
                                  coamps_dataset=coamps_datasets, sam_datasets=sam_datasets,
                                  wrf_datasets=wrf_datasets, r408_dataset=r408_datasets, hoc_dataset=hoc_datasets,
                                  e3sm_datasets=e3sm_file, cam_datasets=cam_file,
                                  time_height=self.time_height, anim=self.animation)
            self.panels.extend(temp_group.panels)

        # Convert panels to difference panels if user passed in --diff <<folder>>
        self.diff_panels = []
        if self.diff_datasets is not None:
            # Loop over the VariableGroup classes listed in the 'var_groups' entry
            # for this case in config/Case_definitions.py and create an instance of each of the listed VariableGroups
            for VarGroup in self.var_groups:
                # Call the __init__ function of the VarGroup class and, by doing this, create an instance of it
                diff_group = VarGroup(self, clubb_datasets=self.diff_datasets, sam_file=les_file,
                                      coamps_file=coamps_datasets, cam_file=cam_file, sam_datasets=sam_datasets,
                                      r408_file=r408_datasets, hoc_dataset=hoc_datasets, e3sm_datasets=e3sm_file)
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

        if self.plot_budgets:
            # Create an instance of a budgets VariableGroup. By default, this is VariableGroupBaseBudgets,
            # but for SAM data, use VariableGroupSamBudgets
            if self.clubb_datasets is not None and len(self.clubb_datasets) != 0:
                # for folders_datasets in self.clubb_datasets.values():
                for input_folder in self.clubb_datasets:
                    folder_name = os.path.basename(input_folder)
                    budget_variables = VariableGroupBaseBudgets(self,
                                                                clubb_datasets={folder_name:clubb_datasets[input_folder]}, anim=self.animation)
                    self.panels.extend(budget_variables.panels)
            if wrf_datasets is not None and len(wrf_datasets) != 0:
                # for folders_datasets in wrf_datasets.values():
                #     budget_variables = VariableGroupBaseBudgets(self, wrf_datasets=folders_datasets)
                for input_folder in wrf_datasets:
                    folder_name = os.path.basename(input_folder)
                    budget_variables = VariableGroupBaseBudgets(self,
                                                                wrf_datasets={folder_name:wrf_datasets[input_folder]}, anim=self.animation)
                    self.panels.extend(budget_variables.panels)
            if e3sm_file is not None and len(e3sm_file) != 0:
                for dataset_name in e3sm_file:
                    # E3SM dataset must be wrapped in the same form as the clubb datasets
                    e3sm_budgets = VariableGroupBaseBudgets(self,
                                                            e3sm_datasets={dataset_name: e3sm_file[dataset_name]}, anim=self.animation)
                    self.panels.extend(e3sm_budgets.panels)
            if sam_datasets is not None and len(sam_datasets) != 0:
                # for dataset in sam_datasets.values():
                for input_folder in sam_datasets:
                    folder_name = os.path.basename(input_folder)
                    budget_variables = VariableGroupSamBudgets(self,
                                                               sam_datasets={folder_name:sam_datasets[input_folder]}, anim=self.animation)
                # sam_budgets = VariableGroupSamBudgets(self, sam_datasets=sam_datasets, anim=self.anim)
                    self.panels.extend(budget_variables.panels)

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
            if panel.panel_type == panel.TYPE_BUDGET:
                plot_paired_lines = False
            panel.plot(output_folder, self.name, replace_images=replace_images, no_legends=no_legends,
                       thin_lines=thin_lines, alphabetic_id=alphabetic_id, paired_plots=plot_paired_lines)
            curr_panel_num += 1
            print("\r\tplotted  ", curr_panel_num, " of ", num_plots, " | ", panel.title)

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
