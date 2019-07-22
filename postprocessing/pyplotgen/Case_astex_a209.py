from pyplotgen.Case import Case
from pyplotgen.DataReader import DataReader
from pyplotgen.VariableGroupBase import VariableGroupBase
from pyplotgen.VariableGroupCorrelations import VariableGroupCorrelations
from pyplotgen.VariableGroupIceMP import VariableGroupIceMP
from pyplotgen.VariableGroupKKMP import VariableGroupKKMP
from pyplotgen.VariableGroupLiquidMP import VariableGroupLiquidMP
from pyplotgen.VariableGroupWs import VariableGroupWs


class Case_astex_a209(Case):
    '''

    '''

    name = 'astex_a209'

    def __init__(self, ncdf_files, plot_sam = True):
        '''

        '''
        self.name = Case_astex_a209.name
        sec_per_min = 60
        self.start_time = 2340 #* sec_per_min
        self.end_time = 2400 #* sec_per_min
        self.height_min_value = 0
        self.height_max_value = 6000
        self.enabled = True
        self.ncdf_files = ncdf_files
        self.blacklisted_variables = []
        sam_file = None
        # Astex_a209 doesn't plot SAM

        # if plot_sam:
        #     datareader = DataReader()
        #     sam_file = datareader.__loadNcFile__(
        #         "/home/nicolas/sam_benchmark_runs/JULY_2017/LBA_128kmx128kmx128_1km_Morrison/LBA_128kmx128kmx128_1km_Morrison.nc")
        base_variables = VariableGroupBase(self.ncdf_files, self, sam_file=sam_file)
        w_variables = VariableGroupWs(self.ncdf_files, self, sam_file=sam_file)
        # ice_variables = VariableGroupIceMP(self.ncdf_files, self, sam_file=sam_file)
        liquid_variables = VariableGroupLiquidMP(self.ncdf_files, self, sam_file=sam_file)
        corr_variables = VariableGroupCorrelations(self.ncdf_files, self, sam_file=sam_file)
        kk_variables = VariableGroupKKMP(self.ncdf_files, self, sam_file=sam_file)


        self.panel_groups = [base_variables, w_variables, liquid_variables, corr_variables, kk_variables]
        self.panels = []

        for panelgroup in self.panel_groups:
            for panel in panelgroup.panels:
                self.panels.append(panel)
