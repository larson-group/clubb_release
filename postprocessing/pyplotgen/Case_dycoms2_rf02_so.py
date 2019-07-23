from pyplotgen.Case import Case
from pyplotgen.DataReader import DataReader
from pyplotgen.VariableGroupBase import VariableGroupBase
from pyplotgen.VariableGroupCorrelations import VariableGroupCorrelations
from pyplotgen.VariableGroupKKMP import VariableGroupKKMP
from pyplotgen.VariableGroupLiquidMP import VariableGroupLiquidMP
from pyplotgen.VariableGroupWs import VariableGroupWs


class Case_dycoms2_rf02_so(Case):
    '''

    '''
    name = 'dycoms2_rf02_so'
    def __init__(self, ncdf_files, plot_sam = True):
        '''

        '''
        self.name = Case_dycoms2_rf02_so.name
        self.start_time = 301
        self.end_time = 360
        self.height_min_value = 0
        self.height_max_value = 1200
        self.enabled = True
        self.ncdf_files = ncdf_files
        self.blacklisted_variables = ['wprrp', 'wpNrp']
        sam_file = None
        if plot_sam:
            datareader = DataReader()
            sam_file = datareader.__loadNcFile__(
                "/home/nicolas/sam_benchmark_runs/JULY_2017/DYCOMS_RF02_128x128x96_nodr_sed/DYCOMS_RF02_128x128x96_nodr_sed.nc")
        base_variables = VariableGroupBase(self.ncdf_files, self, sam_file=sam_file)
        w_variables = VariableGroupWs(self.ncdf_files, self, sam_file=sam_file)
        # ice_variables = VariableGroupIceMP(self.ncdf_files, self, sam_file=sam_file)
        liquid_variables = VariableGroupLiquidMP(self.ncdf_files, self, sam_file=sam_file)
        # corr_variables = VariableGroupCorrelations(self.ncdf_files, self, sam_file=sam_file)
        kk_variables = VariableGroupKKMP(self.ncdf_files, self, sam_file=sam_file)
        self.panel_groups = [base_variables, w_variables, liquid_variables, kk_variables]
        self.panels = []

        for panelgroup in self.panel_groups:
            for panel in panelgroup.panels:
                self.panels.append(panel)
