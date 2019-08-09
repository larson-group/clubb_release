import math

import numpy as np

from pyplotgen.DataReader import DataReader
from pyplotgen.Panel import Panel
from pyplotgen.VariableGroupBaseBudgets import VariableGroupBaseBudgets


class Case:
    '''
    A case is basically a collection of variable groups and
    panels that share attributes like start/end time and get plotted together.
    '''
    def __init__(self, case_definition, ncdf_datasets, plot_sam = False, plot_budgets = False, diff_datasets=None):
        '''
        Initalize a case
        :param ncdf_datasets: List of Dataset objects containing netcdf data from clubb
        '''
        self.name = case_definition['name']
        self.start_time = case_definition['start_time']
        self.end_time = case_definition['end_time']
        self.enabled = case_definition['enabled']
        self.var_groups = case_definition['var_groups']
        self.height_min_value = case_definition['height_min_value']
        self.height_max_value = case_definition['height_max_value']
        self.ncdf_datasets = ncdf_datasets
        self.blacklisted_variables = case_definition['blacklisted_vars']
        self.plot_budgets = plot_budgets
        self.diff_datasets = diff_datasets
        if 'disable_budgets' in case_definition.keys() and case_definition['disable_budgets'] is True:
            self.plot_budgets = False
        sam_file = None
        if plot_sam and case_definition['sam_file'] is not None:
            datareader = DataReader()
            sam_file = datareader.__loadNcFile__(case_definition['sam_file'])

        self.panels = []
        self.diff_panels = []
        for group in self.var_groups:
            temp_group = group(self.ncdf_datasets, self, sam_file=sam_file)

            for panel in temp_group.panels :
                self.panels.append(panel)
        # Convert panels to difference panels if user passed in --diff <<folder>>
        if self.diff_datasets is not None:
            for group in self.var_groups:
                diff_group = group(self.diff_datasets, self, sam_file=sam_file)
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
            budget_variables = VariableGroupBaseBudgets(self.ncdf_datasets, self)
            self.panels.extend(budget_variables.panels)

    def getDiffLinesBetweenPanels(self, panelA, panelB, get_y_diff = False):
        '''

        :param panelA:
        :param panelB:
        :return:
        '''
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
        '''
        Returns an array containing the difference between the two arrays.
        arrB is subtracted from arrA and the values are NOT given as absolute value.
        If one array has fewer entries than another, the value of the longer array will be
        appended (i.e. this uses the fill-zeros philosophy).

        :param arrA: The first array (usually contains larger values)
        :param arrB: The second array (usually contains smaller values)
        :return: a numpy array containing arrA - arrB
        '''

        # Fill zeros on smallest array
        zerosA = [0 for i in range(len(arrA), len(arrB))]
        zerosB = [0 for i in range(len(arrB), len(arrA))]
        arrA = np.append(arrA, zerosA)
        arrB = np.append(arrB, zerosB)

        diff_line = arrA - arrB
        diff_line = [abs(value) for value in diff_line] # ensure every difference is positive
        diff_line = np.asarray(diff_line)
        return diff_line


    def plot(self, output_folder, replace_images = False, no_legends = False, thin_lines = False):
        '''
        Plot all panels associated with the case
        :param casename: str name of the case
        :return:
        '''
        print("\n")
        num_plots = len(self.panels)
        curr_panel_num = 0
        for panel in self.panels:
            print("\r\tplotting ",  curr_panel_num, " of ", num_plots, " | ", panel.title, end="")
            panel.plot(output_folder, self.name, replace_images=replace_images, no_legends=no_legends, thin_lines=thin_lines)
            curr_panel_num += 1
            print("\r\tplotted  ", curr_panel_num, " of ", num_plots, " | ", panel.title)
