"""
:author: Nicolas Strike
:date: Early 2019
"""

from src.DataReader import NetCdfVariable


class Line:
    """
    This class holds information commonly passed around when using plot dependent_data.
    Data included is anything necessary to graph a line using pyplot, including the
    numerical values, line format style, and a label for the plot. In the future we
    would like to seperate the dependent_data from the format.
    """
    def __init__(self, x_data, y_data, line_format = "", label ="Unlabeled plot"):
        """
        Create a new line
        :param x_data: list of values to plot along x axis, must have same length as y
        :param y_data: list of values to plot along y axis, must have same length as x
        :param line_format: A str containing the format for the line plot. See pyplot docs for more info.
        :param label: name of the line (e.g. 'current clubb')
        """
        if len(x_data) != len(y_data):
            raise ValueError("The size of X(" + str(len(x_data)) + ") is not the same as the size of Y(" +
                             str(len(y_data)) + ") for the \"" + label + "\" line.")
        if isinstance(x_data, NetCdfVariable):
            x_data = x_data.dependent_data
        if isinstance(y_data, NetCdfVariable):
            y_data = y_data.dependent_data
        self.x = x_data
        self.y = y_data
        self.line_format = line_format
        self.label = label
