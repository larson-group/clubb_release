'''
This class is responsible for the actual plotting
of data onto graphs (aka panels). It uses matplotlib's
pyplot to generate images with data read in from
the DataReader class.

:author: Nicolas Strike
:date: Early 2019
'''

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
        :param data: The data to plot with the case
        :return:
        '''
        line_format = "r--"
        plt.figure(1)
        plt.subplot(111)
        plt.plot(data.x_values, data.y_values, line_format)
        plt.title(data.title)
        plt.ylabel(data.y_title)
        plt.xlabel(data.x_title)
        plt.show()