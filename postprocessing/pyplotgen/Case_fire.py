from pyplotgen.Case import Case
from pyplotgen.DataReader import DataReader
from pyplotgen.VariableGroupBase import VariableGroupBase
from pyplotgen.VariableGroupBaseBudgets import VariableGroupBaseBudgets
from pyplotgen.VariableGroupWs import VariableGroupWs


class Case_fire(Case):
    '''

    '''
    name = 'fire'
    def __init__(self, ncdf_files, plot_sam = True, plot_budgets = False):
        '''

        '''
        self.name = Case_fire.name
        self.plot_budgets = False #  plot_budgets
        self.start_time = 61
        self.end_time = 120
        self.height_min_value = 0
        self.height_max_value = 1000
        self.enabled = True
        self.ncdf_files = ncdf_files
        self.blacklisted_variables = []
        sam_file = None
        # TODO LES filename unkown
        # if plot_sam:
        #     datareader = DataReader()
        #     sam_file = datareader.__loadNcFile__(
        #         "/home/nicolas/sam_benchmark_runs/JULY_2017/BOMEX_64x64x75/BOMEX_64x64x75_100m_40m_1s.nc")
        base_variables = VariableGroupBase(self.ncdf_files, self, sam_file=sam_file)
        # TODO budget variables are not in the nc outputut
        # budget_variables = VariableGroupBaseBudgets(ncdf_files, self)
        w_variables = VariableGroupWs(self.ncdf_files, self, sam_file=sam_file)
        self.panel_groups = [base_variables, w_variables]
        self.panels = []

        if self.plot_budgets:
            budget_variables = VariableGroupBaseBudgets(ncdf_files, self)
            self.panels.extend(budget_variables.panels)

        for panelgroup in self.panel_groups:
            for panel in panelgroup.panels:
                self.panels.append(panel)
