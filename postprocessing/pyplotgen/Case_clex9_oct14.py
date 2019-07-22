from pyplotgen.Case import Case
from pyplotgen.DataReader import DataReader
from pyplotgen.VariableGroupBase import VariableGroupBase
from pyplotgen.VariableGroupCorrelations import VariableGroupCorrelations
from pyplotgen.VariableGroupIceMP import VariableGroupIceMP
from pyplotgen.VariableGroupKKMP import VariableGroupKKMP
from pyplotgen.VariableGroupLiquidMP import VariableGroupLiquidMP
from pyplotgen.VariableGroupWs import VariableGroupWs


class Case_clex9_oct14(Case):
    '''

    '''
    name = 'clex9_oct14'
    def __init__(self, ncdf_files, plot_sam = True):
        '''

        '''
        self.name = Case_clex9_oct14.name
        sec_per_min = 60
        self.start_time = 181 #* sec_per_min
        self.end_time = 240 #* sec_per_min
        self.height_min_value = 2188
        self.height_max_value = 6688
        self.enabled = True
        self.ncdf_files = ncdf_files
        self.blacklisted_variables = []
        sam_file = None
        # TODO Sam filename unknown
        # if plot_sam:
        #     datareader = DataReader()
        #     sam_file = datareader.__loadNcFile__(
        #         "/home/nicolas/sam_benchmark_runs/SAM6.6/CLOUD_FEEDBACK_s12/ctl_s12_96x96x192_25m_DRZ_N100_fixnudge.nc")
        base_variables = VariableGroupBase(self.ncdf_files, self, sam_file=sam_file)
        # w_variables = VariableGroupWs(self.ncdf_files, self, sam_file=sam_file)
        ice_variables = VariableGroupIceMP(self.ncdf_files, self, sam_file=sam_file)
        liquid_variables = VariableGroupLiquidMP(self.ncdf_files, self, sam_file=sam_file)
        # corr_variables = VariableGroupCorrelations(self.ncdf_files, self, sam_file=sam_file)
        # kk_variables = VariableGroupKKMP(self.ncdf_files, self, sam_file=sam_file)

        self.panel_groups = [base_variables, ice_variables, liquid_variables]
        self.panels = []

        for panelgroup in self.panel_groups:
            for panel in panelgroup.panels:
                self.panels.append(panel)
