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
from datetime import datetime

import matplotlib.pyplot as plt

from Panel import Panel
from config import Style_definitions
from src.interoperability import clean_path


class ContourPanel(Panel):
    """
    ContourPanel class derived from Panel
    The difference is that the plot routine will generate a timeheight plot instead of a standard profile plot.
    """

    def __init__(self, plots, panel_type="profile", title="Unnamed panel", dependent_title="dependent variable",
                 sci_scale = None, centered = False):
        """
        Init, pretty much same as Panel init

        :param bla:
        :return: New object
        """
        super(self, plots, Panel.TYPE_TIMEHEIGHT, title, dependent_title, sci_scale, centered)

    def plot(self, output_folder, casename, replace_images = False, no_legends = True, thin_lines = False,
             alphabetic_id = '', paired_plots = True):
        """
        Generate a single contourf plot from the given data

        :param output_folder: String containing path to folder in which the image files should be created
        :param casename: The name of the case that is plotted in this panel
        :param replace_images: Switch to tell pyplotgen if existing files should be overwritten
        :param alphabetic_id: A string printed into the Panel at coordinates (.9,.9) as an identifier.
        :return: None
        """
        plt.subplot(111)

        # Set font sizes
        plt.rc('font', size=Style_definitions.DEFAULT_TEXT_SIZE)          # controls default text sizes
        plt.rc('axes', titlesize=Style_definitions.AXES_TITLE_FONT_SIZE)     # fontsize of the axes title
        plt.rc('axes', labelsize=Style_definitions.AXES_LABEL_FONT_SIZE)    # fontsize of the x and y labels
        plt.rc('xtick', labelsize=Style_definitions.X_TICKMARK_FONT_SIZE)    # fontsize of the tick labels
        plt.rc('ytick', labelsize=Style_definitions.Y_TICKMARK_FONT_SIZE)    # fontsize of the tick labels
        plt.rc('legend', fontsize=Style_definitions.LEGEND_FONT_SIZE)    # legend fontsize
        plt.rc('figure', titlesize=Style_definitions.TITLE_TEXT_SIZE)  # fontsize of the figure title

        # Get Contour object
        var = self.all_plots[0]
        x_data = var.x
        y_data = var.y
        c_data = var.data
        cmap = var.colors

        cs = plt.contourf(x_data, y_data, c_data, cmap=cmap)
        plt.colorbar(cs)
        plt.title(self.title)
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
        filename = "timeheight_"+ str(datetime.now())

        filename = self.__removeInvalidFilenameChars__(filename)
        # Concatenate with output foldername
        rel_filename = output_folder + "/" +casename+'/' + filename
        rel_filename = clean_path(rel_filename)
        # Save image file
        if replace_images is True or not os.path.isfile(rel_filename+Panel.EXTENSION):
            plt.savefig(rel_filename+Panel.EXTENSION)
        else: # os.path.isfile(rel_filename + Panel.EXTENSION) and replace_images is False:
            print("\n\tImage " + rel_filename+Panel.EXTENSION+
                  ' already exists. To overwrite this image during runtime pass in the --replace (-r) parameter.')
        plt.close()