from pyplotgen.Case import Case
from pyplotgen.DataReader import DataReader
from pyplotgen.VariableGroupBase import VariableGroupBase
from pyplotgen.VariableGroupBaseBudgets import VariableGroupBaseBudgets
from pyplotgen.VariableGroupCorrelations import VariableGroupCorrelations
from pyplotgen.VariableGroupIceMP import VariableGroupIceMP
from pyplotgen.VariableGroupKKMP import VariableGroupKKMP
from pyplotgen.VariableGroupLiquidMP import VariableGroupLiquidMP
from pyplotgen.VariableGroupWs import VariableGroupWs


class Case_mc3e(Case):
    '''

    '''
    name = 'mc3e'
    def __init__(self, ncdf_files, plot_sam = True, plot_budgets = False):
        '''

        '''
        self.plot_budgets = plot_budgets
        self.start_time = 1
        self.end_time = 64800
        self.height_min_value = 0
        self.height_max_value = 18000
        self.enabled = True
        self.ncdf_files = ncdf_files
        self.blacklisted_variables = ['rtp3', 'thlp3', 'rtpthvp', 'thlpthvp', 'Ngm', 'wprrp', 'wpNrp']
        sam_file = None
        if plot_sam:
            datareader = DataReader()
            sam_file = datareader.__loadNcFile__(
                "/home/nicolas/sam_benchmark_runs/MC3E_r1359_128x128x128_1km_Morrison/MC3E.nc")
        base_variables = VariableGroupBase(self.ncdf_files, self, sam_file=sam_file)
        # w_variables = VariableGroupWs(self.ncdf_files, self, sam_file=sam_file)
        ice_variables = VariableGroupIceMP(self.ncdf_files, self, sam_file=sam_file)
        liquid_variables = VariableGroupLiquidMP(self.ncdf_files, self, sam_file=sam_file)
        # corr_variables = VariableGroupCorrelations(self.ncdf_files, self, sam_file=sam_file)
        # kk_variables = VariableGroupKKMP(self.ncdf_files, self, sam_file=sam_file)

        self.panel_groups = [base_variables, ice_variables, liquid_variables]
        self.panels = []

        if self.plot_budgets:
            budget_variables = VariableGroupBaseBudgets(ncdf_files, self)
            self.panels.extend(budget_variables.panels)

        for panelgroup in self.panel_groups:
            for panel in panelgroup.panels:
                self.panels.append(panel)
