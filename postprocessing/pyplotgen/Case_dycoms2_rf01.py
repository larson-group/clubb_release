from pyplotgen.Case import Case
from pyplotgen.DataReader import DataReader
from pyplotgen.VariableGroupBase import VariableGroupBase
from pyplotgen.VariableGroupBaseBudgets import VariableGroupBaseBudgets
from pyplotgen.VariableGroupWs import VariableGroupWs


class Case_dycoms2_rf01(Case):
    '''

    '''
    name = 'dycoms2_rf01'
    def __init__(self, ncdf_files, plot_sam = True, plot_budgets = False):
        '''

        '''
        self.name = Case_dycoms2_rf01.name
        self.plot_budgets = plot_budgets
        self.start_time = 181
        self.end_time = 240
        self.height_min_value = 0
        self.height_max_value = 1200
        self.enabled = True
        self.ncdf_files = ncdf_files
        self.blacklisted_variables = []
        sam_file = None
        if plot_sam:
            datareader = DataReader()
            sam_file = datareader.__loadNcFile__(
                "/home/nicolas/sam_benchmark_runs/DYCOMS_RF01_96x96x320.nc")
        base_variables = VariableGroupBase(self.ncdf_files, self, sam_file=sam_file)

        w_variables = VariableGroupWs(self.ncdf_files, self, sam_file=sam_file)
        self.panel_groups = [base_variables, w_variables]
        self.panels = []

        if self.plot_budgets:
            budget_variables = VariableGroupBaseBudgets(ncdf_files, self)
            self.panels.extend(budget_variables.panels)

        for panelgroup in self.panel_groups:
            for panel in panelgroup.panels:
                self.panels.append(panel)
