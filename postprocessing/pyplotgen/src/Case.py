"""
:author: Nicolas Strike
:date: 2019
"""

import numpy as np

from config.VariableGroupBaseBudgets import VariableGroupBaseBudgets
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

    def __init__(self, case_definition, ncdf_datasets, diff_datasets=None,
                 plot_les=False, plot_budgets=False, plot_r408=False, plot_hoc=False, e3sm_dirs =[]):
        """
        Initialize a Case

        :param case_definition: dict containing case specific elements. These are pulled in from Case_definitions.py,
            see Case_definitions.py for details on how to structure the dict
        :param ncdf_datasets: dict containing Dataset objects holding the data needed for the case. The key for each
            value/Dataset in the dict is set to the ext provided in the filename (e.g. sfc, zt, zm)
        :param plot_les: If True pyplotgen plots LES lines, if False pyplotgen does not plot LES lines
        :param plot_budgets: If True pyplotgen will plot Budgets in addition to the other plots
                If False (default), pyplotgen will not plot budgets
        :param diff_datasets: Takes in netcdf datasets. If datasets are passed in, pyplotgen will plot the numeric
            difference between the folder passed in an the input folder.
        :param plot_r408: If True, pyplotgen will plot the Chris Golaz 'best ever' clubb r408 data lines
                If False (default), pyplotgen will not plot the Chris Golaz 'best ever' clubb r408 data lines
        :param plot_hoc: If True, pyplotgen will plot the HOC 2005 data lines
                If False (default), pyplotgen will not plot the HOC 2005 data lines
        """
        self.name = case_definition['name']
        self.start_time = case_definition['start_time']
        self.end_time = case_definition['end_time']
        self.var_groups = case_definition['var_groups']
        self.height_min_value = case_definition['height_min_value']
        self.height_max_value = case_definition['height_max_value']
        self.ncdf_datasets = ncdf_datasets
        self.blacklisted_variables = case_definition['blacklisted_vars']
        self.plot_budgets = plot_budgets
        self.plot_r408 = plot_r408
        self.plot_hoc = plot_hoc
        self.plot_e3sm = e3sm_dirs
        self.diff_datasets = diff_datasets
        self.next_panel_alphabetic_id_code = 97
        if 'disable_budgets' in case_definition.keys() and case_definition['disable_budgets'] is True:
            self.plot_budgets = False

        sam_file = None
        if plot_les and case_definition['sam_file'] is not None:
            datareader = DataReader()
            sam_file = datareader.__loadNcFile__(case_definition['sam_file'])

        e3sm_file = {}
        if len(e3sm_dirs) != 0 and e3sm_dirs != None and case_definition['e3sm_file'] != None:
            datareader = DataReader()
            for foldername in e3sm_dirs:
                e3sm_filename = foldername + case_definition['e3sm_file']
                e3sm_file[foldername] = datareader.__loadNcFile__(e3sm_filename)

        coamps_datasets = {}
        if plot_les and case_definition['coamps_file'] is not None:
            datareader = DataReader()
            coamps_filenames = case_definition['coamps_file']
            for type_ext in coamps_filenames:
                temp_coamps_dataset = datareader.__loadNcFile__(coamps_filenames[type_ext])
                coamps_datasets[type_ext] = temp_coamps_dataset
        else:
            coamps_datasets = None

        r408_datasets = {}
        if plot_r408 and case_definition['r408_file'] is not None:
            datareader = DataReader()
            r408_filenames = case_definition['r408_file']
            for type_ext in r408_filenames:
                temp_r408_dataset = datareader.__loadNcFile__(r408_filenames[type_ext])
                r408_datasets[type_ext] = temp_r408_dataset
        else:
            r408_datasets = None

        hoc_datasets = {}
        if plot_hoc and case_definition['hoc_file'] is not None:
            datareader = DataReader()
            hoc_filenames = case_definition['hoc_file']
            for type_ext in hoc_filenames:
                temp_hoc_dataset = datareader.__loadNcFile__(hoc_filenames[type_ext])
                hoc_datasets[type_ext] = temp_hoc_dataset
        else:
            hoc_datasets = None

        self.panels = []
        self.diff_panels = []
        for VarGroup in self.var_groups:
            temp_group = VarGroup(self.ncdf_datasets, self, sam_file=sam_file, coamps_file=coamps_datasets,
                                  r408_dataset=r408_datasets, hoc_dataset=hoc_datasets, e3sm_datasets = e3sm_file)

            for panel in temp_group.panels:
                self.panels.append(panel)

        # Convert panels to difference panels if user passed in --diff <<folder>>
        if self.diff_datasets is not None:
            for VarGroup in self.var_groups:
                diff_group = VarGroup(self.diff_datasets, self, sam_file=sam_file, coamps_file=coamps_datasets,
                                      r408_file=r408_datasets, hoc_dataset=hoc_datasets, e3sm_datasets = e3sm_file)
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
            if self.ncdf_datasets is not None and len(self.ncdf_datasets) is not 0:
                budget_variables = VariableGroupBaseBudgets(self.ncdf_datasets, self)
                self.panels.extend(budget_variables.panels)
            if e3sm_file != None and len(e3sm_file) is not 0:
                for dataset_name in e3sm_file:
                    e3sm_budgets = VariableGroupBaseBudgets({dataset_name:e3sm_file[dataset_name]}, self) # E3SM dataset must be wrapped in the same form as the clubb datasets
                    self.panels.extend(e3sm_budgets.panels)


    def getDiffLinesBetweenPanels(self, panelA, panelB, get_y_diff=False):
        """
        Given two panels of type Panel, this function calculates the numerical
        difference between each line of the two panels and returns the numerical
        values as a list [] for each line in the panel.

        This function assumes that the data lines are in the same order on each panel
        and that each panel has the same number of lines

        :param panelA: The first panel for comparison
        :param panelB: The second panel for comparison
        :param get_y_diff: Default behavior is to return the difference along the x-axis, setting this to True returns
            the difference on the y-axis
        :return: A 2D list containing the numerical values for the difference between each line in the given panels.
        """
        linesA = panelA.all_plots
        linesB = panelB.all_plots
        newLines = linesA
        for idx in range(len(linesA)):
            if get_y_diff:
                a_data = linesA[idx].y
                b_data = linesB[idx].y
                diff_line_data = self.__getArrayDiff__(a_data, b_data)
                newLines[idx].y = (diff_line_data)
            else:
                a_data = linesA[idx].x
                b_data = linesB[idx].x
                diff_line_data = self.__getArrayDiff__(a_data, b_data)
                newLines[idx].x = (diff_line_data)
        return newLines

    def __getArrayDiff__(self, arrA, arrB):
        """
        Returns an array containing the difference between the two arrays.
        arrB is subtracted from arrA and the values are NOT given as absolute value.
        If one array has fewer entries than another, the value of the longer array will be
        appended (i.e. this uses the fill-zeros philosophy).

        :param arrA: The first array (usually contains larger values)
        :param arrB: The second array (usually contains smaller values)
        :return: a numpy array containing arrA - arrB
        """

        # Fill zeros on smallest array
        zerosA = [0 for i in range(len(arrA), len(arrB))]
        zerosB = [0 for i in range(len(arrB), len(arrA))]
        arrA = np.append(arrA, zerosA)
        arrB = np.append(arrB, zerosB)

        diff_line = arrA - arrB
        diff_line = [abs(value) for value in diff_line]  # ensure every difference is positive
        diff_line = np.asarray(diff_line)
        return diff_line

    def plot(self, output_folder, replace_images=False, no_legends=False, thin_lines=False, show_alphabetic_id=False):
        """
        Plot all panels associated with the case, these will be saved to a .jpg file in the <<output>>/<<casename>> folder
        :param output_folder: absolute name of the folder to save output into.
        :param replace_images: If True, pyplotgen will overwrite images with the same name.
            If False (default), pyplotgen will add a timestamp to the end of every filename (even if there's no filename conflict)
        :param no_legends: If True, pyplotgen will not include a legend on output graphs.
        If False (default), legends will be displayed.
        :param thin_lines: If True, lines plotted will be much thinner than usual.
        False (default) lines are plotted according to the thickness defined in config/Style_definitions.py
        :param show_alphabetic_id: If True, pyplotgen will add an alphabetic label to the top right corner of each plot.
        These labels will rotate through a-z incrementally. If there are more than 26 plots, it will rotate 2 dimentionally,
        e.g. (aa), (ab), (ac),...,(ba),(bb),(bc) and etc. If show_alphabetic_id is False (default), this label is not displayed. The
        behavior for rotations greater than zz is not defined, and although pyplotgen won't crash it may start to use
        weird characters. The rotation resets between each case, e.g. if one case ends on label (ad), the next case will start on (a).
        :return:
        """

        print("\n")
        num_plots = len(self.panels)
        curr_panel_num = 0
        for panel in self.panels:
            print("\r\tplotting ", curr_panel_num, " of ", num_plots, " | ", panel.title)
            if show_alphabetic_id:
                alphabetic_id = self.__getNextAlphabeticID()
            else:
                alphabetic_id = ""
            panel.plot(output_folder, self.name, replace_images=replace_images, no_legends=no_legends,
                       thin_lines=thin_lines, alphabetic_id=alphabetic_id)
            curr_panel_num += 1
            print("\r\tplotted  ", curr_panel_num, " of ", num_plots, " | ", panel.title)

    def __getNextAlphabeticID(self):
        """
        When --show-alphabetic-id is passed in as a run parameter, pyplotgen will add an alphabetic label to each plot.
        These labels are a 1 or 2d rotation through the alphabet (automatically converts to 2d if more than 26 labels are needed).
        E.g. this method will return 'a' the first time, 'b' the second time, 'aa' the 27th time, 'ab' the 28th time, and etc.
        This function returns the next label as a string, and keeps track of each call to this method. Call count tracking
        is specific per case instance, e.g. if one case plots the label 'bb' last the next case will plot 'a' as the first
        label instead of 'bc'.
        :return:
        """
        a = 97
        z = 122
        num_letters_a_to_z = 26
        if self.next_panel_alphabetic_id_code >= a and self.next_panel_alphabetic_id_code <= z:
            letter = chr(self.next_panel_alphabetic_id_code)
            self.next_panel_alphabetic_id_code += 1
            return letter
        else:
            first_letter_offset = int((self.next_panel_alphabetic_id_code - a) / num_letters_a_to_z) - 1
            first_letter = chr(a + first_letter_offset)

            second_letter_offset = int((self.next_panel_alphabetic_id_code - a * int(
                self.next_panel_alphabetic_id_code / a)) % num_letters_a_to_z)
            second_letter = chr(a + second_letter_offset)

            self.next_panel_alphabetic_id_code += 1
            return first_letter + second_letter
