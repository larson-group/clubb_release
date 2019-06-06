

class Case:
    '''

    '''
    def __init__(self, name, panel_groups, enabled = True, start_time = 0, end_time = None):
        '''

        '''
        self.name = name
        self.enabled = enabled
        self.panel_groups = panel_groups
        self.start_time = start_time
        self.end_time = end_time
        self.panels = []

        for panel in panel_groups:
            self.panels.append(panel)

    def plot(self, plotter, netcdf_data):
        '''

        :param plotter:
        :return:
        '''
        for panel in self.panels:
            panel.plot(plotter, netcdf_data)
