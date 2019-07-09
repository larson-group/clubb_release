from pyplotgen.Case import Case
from pyplotgen.VariableGroupBase import VariableGroupBase
from pyplotgen.VariableGroupIceMP import VariableGroupIceMP
from pyplotgen.VariableGroupLiquidMP import VariableGroupLiquidMP
from pyplotgen.VariableGroupWs import VariableGroupWs


class Case_lba(Case):
    '''

    '''
    def __init__(self, ncdf_files, sam_file = None):
        '''

        '''
        self.name = "lba"
        sec_per_min = 60
        self.averaging_start_time = 300 * sec_per_min
        self.averaging_end_time = 360 * sec_per_min
        self.timeseries_start_time = 0 * sec_per_min
        self.timeseries_end_time = 360 * sec_per_min
        self.height_min_value = 0
        self.height_max_value = 12000
        self.enabled = True
        self.ncdf_files = ncdf_files
        base_variables = VariableGroupBase(self.ncdf_files, self, sam_file=sam_file)
        w_variables = VariableGroupWs(self.ncdf_files, self, sam_file=sam_file)
        ice_variables = VariableGroupIceMP(self.ncdf_files, self, sam_file=sam_file)
        liquid_variables = VariableGroupLiquidMP(self.ncdf_files, self, sam_file=sam_file)

        self.panel_groups = [base_variables, w_variables, ice_variables, liquid_variables]
        self.panels = []

        for panelgroup in self.panel_groups:
            for panel in panelgroup.panels:
                self.panels.append(panel)
