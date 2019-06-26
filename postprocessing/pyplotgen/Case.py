from pyplotgen.VariableGroupTest import VariableGroupBase


class Case:
    '''
    A case is basically a collection of variable groups and
    panels that share attributes like start/end time and get plotted together.
    '''
    def __init__(self, ncdf_datasets):
        '''
        Initalize a case
        :param ncdf_datasets: List of Dataset objects containing netcdf data from clubb
        '''
        self.name = "unnamed-case"
        self.start_time = 0
        self.end_time = -1
        self.enabled = True
        self.ncdf_datasets = ncdf_datasets

        self.panel_groups = []
        self.panels = []

        for panel in self.panel_groups.panels:
            self.panels.append(panel)

    def plot(self):
        '''
        Plot all panels associated with the case
        :param casename: str name of the case
        :return:
        '''
        for panel in self.panels:
            panel.plot(self.name)
