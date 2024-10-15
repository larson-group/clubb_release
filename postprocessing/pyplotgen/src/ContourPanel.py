# -*- coding: utf-8 -*-
'''
Child class derived from Panel.
Takes the same variable definitions as a Panel for profile plots.
But needs the unaveraged netcdf data.
The changed plot routine outputs contour plots instead of profile plots.

:author: Steffen Domke
:date: July 2020
'''
import os
import warnings
from datetime import datetime

#TODO temporary fix to suppress warnings related to chi/eta corr vars
import logging
logging.captureWarnings(True)

import matplotlib.pyplot as plt
import numpy as np

from config import Style_definitions
from src.Panel import Panel
from src.interoperability import clean_path


class ContourPanel(Panel):
    """
    ContourPanel class derived from Panel
    The difference is that the plot routine will generate a timeheight plot instead of a standard profile plot.

    For information on the input parameters of this class, please see the documentation for the
    ``__init__()`` method.
    """

    def __init__(self, plots, panel_type=Panel.TYPE_TIMEHEIGHT, title="Unnamed panel",
                 dependent_title="dependent variable"):
        """
        Constructor of the ContourPanel subclass. Will call the constructor of the Panel superclass.

        :param plots: List containing exactly one Contour object to plot onto the panel
        :param panel_type: Type of panel being plotted (must be Panel.TYPE_TIMEHEIGHT for this class)
        :param title: The title of this plot (e.g. 'Liquid water potential tempature')
        :param dependent_title: Label of the dependent axis (labels the x-axis for all panel types except timeseries).
        :param sci_scale: The scale at which to display the x axis (e.g. to scale to 1e-3 set sci_scale=-3).
            If not specified, the matplotlib default sci scaling will be used.
        :param centered: If True, the Panel will be centered around 0.
            Profile plots are usually centered, while budget plots are not.
        """

        # Note: Nones added to make contour plots ignore background-rcm option for now
        #       If, at some point, we want this to work with Contours,
        #       we would need to modify the super call here and kind of copy the code from Panel
        #       that handles background-rcm into ContourPanel.
        super().__init__(plots, bkgrnd_rcm=None, altitude_bkgrnd_rcm=None, start_alt_idx=None, end_alt_idx=None,
                         panel_type=panel_type, title=title, dependent_title=dependent_title, sci_scale=None, centered=False)

    def plot(self, output_folder, casename, replace_images = False, no_legends = True, thin_lines = False,
             alphabetic_id = '', paired_plots = True, image_extension=".png"):
        """
        Generate a single contourf plot from the given data

        :param output_folder: String containing path to folder in which the image files should be created
        :param casename: The name of the case that is plotted in this panel
        :param replace_images: Switch to tell pyplotgen if existing files should be overwritten
        :param alphabetic_id: A string printed into the Panel at coordinates (.9,.9) as an identifier.
        :return: None
        """
        # Suppress deprecation warnings
        with warnings.catch_warnings():
            # Create new figure and axis
            warnings.simplefilter("ignore")
            plt.subplot(111)

        # Set font sizes
        plt.rc('font', size=Style_definitions.DEFAULT_TEXT_SIZE)          # controls default text sizes
        plt.rc('axes', titlesize=Style_definitions.AXES_TITLE_FONT_SIZE)     # fontsize of the axes title
        plt.rc('axes', labelsize=Style_definitions.AXES_LABEL_FONT_SIZE)    # fontsize of the x and y labels
        plt.rc('xtick', labelsize=Style_definitions.X_TICKMARK_FONT_SIZE)    # fontsize of the tick labels
        plt.rc('ytick', labelsize=Style_definitions.Y_TICKMARK_FONT_SIZE)    # fontsize of the tick labels
        plt.rc('legend', fontsize=Style_definitions.LEGEND_FONT_SIZE)    # legend fontsize
        plt.rc('figure', titlesize=Style_definitions.TITLE_TEXT_SIZE)  # fontsize of the figure title

        # For each Contour object stored in self.all_plots generate an individual contourf plot
        for var in self.all_plots:
            x_data = var.x
            y_data = var.y
            c_data = var.data
            x_data, y_data = np.meshgrid(x_data, y_data)
            cmap = var.colors
            label = var.label

            # Set graph size
            plt.figure(figsize=(10,6))

            # Prevent x-axis label from getting cut off
            # plt.gcf().subplots_adjust(bottom=0.15)

            cs = plt.contourf(x_data, y_data, c_data.T, cmap=cmap)
            plt.colorbar(cs)
            plt.title(label + ' - ' + self.title, pad=10)
            plt.xlabel(self.x_title)
            plt.ylabel(self.y_title)

            if alphabetic_id != '':
                ax = plt.gca()
                ax.text(0.9, 0.9, '('+alphabetic_id+')', ha='center', va='center', transform=ax.transAxes,
                               fontsize=Style_definitions.LARGE_FONT_SIZE) # Add letter label to panels

            # Create folders
            # Because os.mkdir("output") can fail and prevent os.mkdir("output/" + casename) from being called we must
            # use two separate try blocks
            try:
                os.mkdir(output_folder)
            except FileExistsError:
                pass # do nothing
            try:
                os.mkdir(output_folder + "/" + casename)
            except FileExistsError:
                pass # do nothing

            # Generate image filename
            filename = "timeheight_"+ str(datetime.now())+ "_" + self.title

            filename = self.__removeInvalidFilenameChars__(filename)
            # Concatenate with output foldername
            relative_filename = output_folder + '/' + casename + '/' + filename
            relative_filename = clean_path(relative_filename)
            # Save image file
            plt.savefig(relative_filename+image_extension)
            plt.close()