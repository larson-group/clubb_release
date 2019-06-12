from pyplotgen.VariableGroupTest import VariableGroupBase


class Case:
    '''

    '''
    def __init__(self, ncdf_file):
        '''

        '''
        self.name = "unnamed-case"
        self.start_time = 0
        self.end_time = -1
        self.enabled = True
        self.ncdf_file = ncdf_file

        self.panel_groups = []
        self.panels = []

        for panel in self.panel_groups.panels:
            self.panels.append(panel)

    def plot(self, casename):
        '''

        :param plotter:
        :return:
        '''
        for panel in self.panels:
            panel.plot(casename)
