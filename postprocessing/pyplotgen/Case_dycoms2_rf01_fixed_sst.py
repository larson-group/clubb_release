from pyplotgen.Case import Case
from pyplotgen.DataReader import DataReader
from pyplotgen.VariableGroupBase import VariableGroupBase
from pyplotgen.VariableGroupBaseBudgets import VariableGroupBaseBudgets
from pyplotgen.VariableGroupWs import VariableGroupWs


class Case_dycoms2_rf01_fixed_sst(Case):
    '''

    '''
    name = 'dycoms2_rf01_fixed_sst'

    def __init__(self, ncdf_files, plot_sam = True):
        '''

        '''
        self.name = Case_dycoms2_rf01_fixed_sst.name
        sec_per_min = 60
        self.start_time = 2520 #* sec_per_min
        self.end_time = 2700 #* sec_per_min
        self.height_min_value = 0
        self.height_max_value = 1200
        self.enabled = True
        self.ncdf_files = ncdf_files

        # rtp3,thlp3,rtpthvp are calculable values, but plotgen doesn't output them so they were ignored
        self.blacklisted_variables = ['rtp3', 'thlp3', 'rtpthvp', 'thlpthvp']
        sam_file = None
        if plot_sam:
            datareader = DataReader()
            sam_file = datareader.__loadNcFile__(
                "/home/nicolas/sam_benchmark_runs/SAM6.6/DYCOMS_RF01_fixed_sst/DYCOMS_RF01_96x96x320_LES_fixed_sst.nc")
        base_variables = VariableGroupBase(self.ncdf_files, self, sam_file=sam_file)
        budget_variables = VariableGroupBaseBudgets(ncdf_files, self)
        # w_variables = VariableGroupWs(self.ncdf_files, self, sam_file=sam_file)
        self.panel_groups = [base_variables, budget_variables]
        self.panels = []

        for panelgroup in self.panel_groups:
            for panel in panelgroup.panels:
                self.panels.append(panel)
