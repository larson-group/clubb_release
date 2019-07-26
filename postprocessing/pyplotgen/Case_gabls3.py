from pyplotgen.Case import Case
from pyplotgen.DataReader import DataReader
from pyplotgen.VariableGroupBase import VariableGroupBase
from pyplotgen.VariableGroupWs import VariableGroupWs


class Case_gabls3(Case):
    '''

    '''
    name = 'gabls3'
    def __init__(self, ncdf_files, plot_sam = True):
        '''

        '''
        self.name = "gabls3"
        sec_per_min = 60
        self.start_time = 1081 #* sec_per_min
        self.end_time = 1200 #* sec_per_min
        self.height_min_value = 0
        self.height_max_value = 4970
        self.enabled = True
        self.ncdf_files = ncdf_files
        # tau_zm not outputted to .nc file
        self.blacklisted_variables = []

        sam_file = None
        # TODO LES filename unkown
        # if plot_sam:
        #     datareader = DataReader()
        #     sam_file = datareader.__loadNcFile__(
        #         "/home/nicolas/sam_benchmark_runs/JULY_2017/BOMEX_64x64x75/BOMEX_64x64x75_100m_40m_1s.nc")
        base_variables = VariableGroupBase(self.ncdf_files, self, sam_file=sam_file)
        budget_variables = VariableGroupBaseBudgets(ncdf_files, self)
        # w_vars removed because they were not outputted in .nc file
        # w_variables = VariableGroupWs(self.ncdf_files, self, sam_file=sam_file)
        self.panel_groups = [base_variables, budget_variables]#, w_variables]
        self.panels = []

        for panelgroup in self.panel_groups:
            for panel in panelgroup.panels:
                self.panels.append(panel)
