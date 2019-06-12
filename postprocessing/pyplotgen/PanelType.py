'''
:author: Nicolas Strike
:date: Mid 2019
'''

class PanelType:
    '''

    '''
    def __init__(self, dependant_varname, blacklisted_panels = None, dependant_title = "dependant variable"):
        self.TYPE_PROFILE = 'profile'
        self.TYPE_BUDGET = 'budget'
        self.TYPE_TIMESERIES = 'timeseries'
        # self.blacklisted_panels = blacklisted_panels
        self.dependant_varname = dependant_varname
        self.dependant_title = dependant_title

    def plot(self, plotter):
        '''

        :param plotter:
        :return:
        '''

