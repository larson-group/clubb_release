# -*- coding: utf-8 -*-
"""
:author: Steffen Domke
:date: June 2020
"""

from src.DataReader import NetCdfVariable


class Contour:
    """
    This class holds any information necessary to create a contour plot using pyplot,
    including the numerical values, format style, and a label for the plot.
    In the future we would like to seperate the dependent_data from the format.

    For information on the input parameters of this class, please see the documentation for the
    ``__init__()`` method.
    """
    def __init__(self, x_data, y_data, c_data, colors = "", label = "Unlabeled plot", line_format = ""):
        """
        Create a new Contour object

        :param x_data: List of values to plot along x axis
        :param y_data: List of values to plot along y axis
        :param c_data: 2-dimensional array of values to plot as contours, must have shape as (len(x_data), len(y_data))
        :param colors: A string specifying the color palette used for the contours. See pyplot docs for more info.
        :param label: The name of the plot (e.g. 'current clubb') to be used in the title
        """
        if isinstance(x_data, NetCdfVariable):
            x_data = x_data.dependent_data
        if isinstance(y_data, NetCdfVariable):
            y_data = y_data.dependent_data
        if isinstance(c_data, NetCdfVariable):
            c_data = c_data.dependent_data
        x_len = len(x_data)
        y_len = len(y_data)
        if c_data.shape != (x_len,y_len):
            raise ValueError("The size of data(" + str(c_data.shape) + ") is not the same as the size of XxY(" +
                             str(len(x_data)) + 'x' + str(len(y_data)) + ") for the \"" + label + "\" contour.")
        self.x = x_data
        self.y = y_data
        self.data = c_data
        self.colors = colors
        self.label = label
        self.line_format = line_format