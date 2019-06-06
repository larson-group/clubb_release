'''
:author: Nicolas Strike
:date: Mid 2019
'''

class PanelType:
    '''

    '''
    def __init__(self, panels = None, blacklisted_panels = None, x_varname = "", y_varname = "", line_format = ""):
        self.panels = panels
        self.blacklisted_panels = blacklisted_panels
        self.x_varname = x_varname
        self.y_varname = y_varname
        self.line_format = line_format

    def plot(self, plotter):
        '''

        :param plotter:
        :return:
        '''

