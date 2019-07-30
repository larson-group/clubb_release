from pyplotgen.DataReader import DataReader
from pyplotgen.VariableGroupBase import VariableGroupBase


class Case:
    '''
    A case is basically a collection of variable groups and
    panels that share attributes like start/end time and get plotted together.
    '''
    def __init__(self, case_definition, ncdf_datasets, plot_sam = False, plot_budgets = False):
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
        sam_file = None
        if plot_sam and case_definition['sam_file'] is not None:
            datareader = DataReader()
            sam_file = datareader.__loadNcFile__(case_definition['sam_file'])

        # if self.plot_budgets:
        #     budget_variables = VariableGroupBaseBudgets(ncdf_files, self)
        #     self.panels.extend(budget_variables.panels)

        self.panels = []
        for group in self.var_groups:
            temp_group =   group(self.ncdf_datasets, self, sam_file=sam_file)
            for panel in temp_group.panels :
                self.panels.append(panel)

    def plot(self, output_folder, replace_images = False, no_legends = False):
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
            panel.plot(output_folder, self.name, replace_images=replace_images, no_legends=no_legends)
            curr_panel_num += 1
            print("\r\tplotted  ", curr_panel_num, " of ", num_plots, " | ", panel.title)
