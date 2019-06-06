from pyplotgen.Lineplot import Lineplot
from pyplotgen.PantelType import PanelType
from pyplotgen.Plotter import Plotter


class Panel:
    '''

    '''
    def __init__(self, base_plot : Lineplot,
                 panel_type : PanelType = None, title : str = "Unnamed panel", overplots = None,
                 x_title = "x axis", y_title = "y axis", line_format = "", blacklisted = False):
        '''

        '''
        # self.x_varname = x_varname
        # self.y_varname = y_varname
        self.base_plot = base_plot
        self.panel_type = panel_type
        self.title = title
        self.overplots = overplots
        self.x_title = x_title
        self.y_title = y_title
        self.blacklisted = blacklisted

    def plot(self, netcdf_data):
        '''

        :param plotter:
        :return:
        '''
        plotter = Plotter()
        plotter.plot(self.base_plot.x.data, self.base_plot.y.data, title=self.title, x_title=self.x_title, y_title=self.y_title)
