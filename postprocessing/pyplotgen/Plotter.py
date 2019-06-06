'''
This class is responsible for the actual plotting
of data onto graphs (aka panels). It uses matplotlib's
pyplot to generate images with data read in from
the DataReader class.

:author: Nicolas Strike
:date: Mid 2019
'''
from collections import namedtuple

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


class Plotter:
    PlotValues = namedtuple("PlotValues", "x_values y_values")
    PlotDetails = namedtuple("PlotDetails", "title x_title y_title")



    def __init__(self):
        pass

    def plot_figure(self, case, variable):
        '''
        Plot a single figure/panel. This is 1 graph but
        may include multiple sets of data (e.g. multiple lines)
        :param case:
        :param data:
        :return:
        '''
        pass
        # title = data.
        # figure = plt.figure(1)
        # plt.subplot(111)
        # plt.plot(thlm_time[11])


    def plot(self, x_values, y_values, title = "unnamed plot", x_title = "x axis", y_title = "y title", line_format = "r--"):
        '''
        Saves a single panel/graph to the output directory specified by the pyplotgen launch paramters

        The x_values and y_values lists must be in order, as they are associated by index

        :param x_values: A list containing the x values at which each y value shall be plotted
        :param y_values: A list containing the y to plot at each x value
        :param title: The title of the panel, e.g. Liquid Water Potential Temperature
        :param x_title: Label for x-axis, e.g. thlm [K]
        :param y_title: Label for y-axis, e.g. Height [m]
        :param line_format: Describes how the line shall appear. See https://matplotlib.org/2.1.2/api/_as_gen/matplotlib.pyplot.plot.html
        :return: n/a
        '''
        plt.figure()
        plt.subplot(111)
        plt.plot(x_values, y_values, line_format)
        plt.title(title)
        plt.ylabel(y_title)
        plt.xlabel(x_title)
        # plt.autoscale(tight=True)
        plt.savefig(title.replace('"', '').replace(',', '').replace(' ', '_') + '.png')
