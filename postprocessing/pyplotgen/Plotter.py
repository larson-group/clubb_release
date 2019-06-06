'''
This class is responsible for the actual plotting
of data onto graphs (aka panels). It uses matplotlib's
pyplot to generate images with data read in from
the DataReader class.

:author: Nicolas Strike
:date: Early 2019
'''
import math
from collections import namedtuple

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


class Plotter:
    PlotValues = namedtuple("PlotValues", "x_values y_values")
    PlotDetails = namedtuple("PlotDetails", "title x_title y_title")



    def __init__(self):
        '''
        Initialize a plotter object TODO
        '''
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


    def plot(self, x_values, y_values, title = "unnamed plot", x_title = "x axis", y_title = "y title"):
        '''
        Plots a case onto a set of graphs (panels)
        :param case: The case to be plotted
        :param values: The values to plot with the case formatted as the namedtuple defined in DataReader.py
        :return:
        '''
        # print(values)
        line_format = "r--"

        plt.figure()
        plt.subplot(111)
        plt.plot(x_values, y_values, line_format)
        plt.title(title)
        plt.ylabel(y_title)
        plt.xlabel(x_title)

        # x_min = values.x_values.min()
        # x_max = values.x_values.max()
        #
        # plt.xlim(math.floor(x_min),math.ceil(x_max))
        # plt.ylim(math.floor(x_min),math.ceil(x_max))

        # plt.autoscale(tight=True)
        # print("Saving file")
        plt.savefig(title.replace('"', '').replace(',', '').replace(' ', '_') + '.png')
        # print("Saved file")
