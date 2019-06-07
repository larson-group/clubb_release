'''
:author: Nicolas Strike
:date: Mid 2019
'''


from pyplotgen.Lineplot import Lineplot
from pyplotgen.PantelType import PanelType
from pyplotgen.Plotter import Plotter


class Panel:
    '''
    Represents an individual panel/graph. Each graph contains a number of details
    specific to it, such as a title. Each graph can be plotted/saved to a file.
    '''
    def __init__(self, base_plot : Lineplot,panel_type : PanelType = None, title : str = "Unnamed panel",
                 overplots = None, x_title = "x axis", y_title = "y axis", line_format = "", blacklisted = False):
        self.base_plot = base_plot
        self.panel_type = panel_type
        self.title = title
        self.overplots = overplots
        self.x_title = x_title
        self.y_title = y_title
        self.blacklisted = blacklisted

    def plot(self):
        '''
        Save the panel's graphical representation to a file in the 'output' folder from launch parameters

        :param plotter:
        :return:
        '''
        plotter = Plotter()
        plotter.plot(self.base_plot.x.data, self.base_plot.y.data, title=self.title, x_title=self.x_title, y_title=self.y_title)
