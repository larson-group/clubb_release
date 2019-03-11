'''
This class is responsible for the actual plotting
of data onto graphs (aka panels). It uses matplotlib's
pyplot to generate images with data read in from
the DataReader class.

:author: Nicolas Strike
:date: Early 2019
'''
import math

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


class Plotter:

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


    def plot(self, data):
        '''
        Plots a case onto a set of graphs (panels)
        :param case: The case to be plotted
        :param data: The data to plot with the case formatted as the namedtuple defined in DataReader.py
        :return:
        '''
        print(data)
        line_format = "r--"

        plt.figure(1)
        plt.subplot(111)
        plt.plot(data.x_values, data.y_values, line_format)
        plt.title(data.title)
        plt.ylabel(data.y_title)
        plt.xlabel(data.x_title)

        x_min = data.x_values.min()
        x_max = data.x_values.max()

        # plt.xlim(math.floor(x_min),math.ceil(x_max))
        # plt.autoscale(tight=True)
        # print("Saving file")
        plt.savefig(data.title.replace('"', '').replace(',','').replace(' ', '_') + '.png')
        # print("Saved file")
