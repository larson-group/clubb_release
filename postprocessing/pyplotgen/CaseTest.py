from pyplotgen.Case import Case
from pyplotgen.VariableGroupTest import VariableGroupBase


class CaseTest(Case):
    '''

    '''
    def __init__(self, ncdf_files, sam_file = None):
        '''

        '''
        self.name = "Base-variables"

        sec_per_min = 60
        self.averaging_start_time = 181 * sec_per_min
        self.averaging_end_time = 240 * sec_per_min
        self.timeseries_start_time = 0 * sec_per_min
        self.timeseries_end_time = 240 * sec_per_min
        self.height_min_value = 0
        self.height_max_value = 1200

        self.start_time = 0
        self.end_time = -1
        self.enabled = True
        self.ncdf_files = ncdf_files

        testGroup = VariableGroupBase(self.ncdf_files, self, sam_file=sam_file)
        self.panel_groups = [testGroup]
        self.panels = []

        for panelgroup in self.panel_groups:
            for panel in panelgroup.panels:
                self.panels.append(panel)

    def plot(self):
        '''

        :param plotter:
        :return:
        '''
        super().plot()