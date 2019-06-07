from pyplotgen import PanelGroup
from pyplotgen.Case import Case

from pyplotgen.PanelGroupTest import PanelGroupTest


class CaseTest(Case):
    '''

    '''
    def __init__(self, ncdf_file):
        '''

        '''
        self.name = "testcase"
        self.start_time = 0
        self.end_time = -1
        self.enabled = True
        self.ncdf_file = ncdf_file

        testGroup = PanelGroupTest(self.ncdf_file)
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