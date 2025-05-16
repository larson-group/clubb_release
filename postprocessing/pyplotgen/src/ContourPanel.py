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

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as mpt
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
            cmap = mpl.colormaps.get_cmap(var.colors)
            vmin, vmax = c_data.min(), c_data.max()
            #vmin, vmax = 0, 0.017 #rtm arm
            #vmin, vmax = 0, 0.0119 #rtm astex
            #vmin, vmax = 0, 0.00296 #rtm gabls2
            
            #vmin, vmax = -0.42337, 0.12668 #wpthlp arm
            #vmin, vmax = -0.49862, 0.06589 #wpthlp astex
            #vmin, vmax = -0.05, 0.008 #wpthlp astex new
            #vmin, vmax = -0.11252, 0.21776 #wpthlp gabls2
            #vmin, vmax = -0.04, 0.07 #wpthlp gabls2 new

            #vmin, vmax = 0, 0.78 #cloud frac arm
            #vmin, vmax = 0, 1.0 #cloud frac astex

            #vmin, vmax = 0, 1.89 #wp2 arm
            #vmin, vmax = 0, 0.49 #wp2 astex
            #vmin, vmax = 0, 1.0 #wp2 gabls2

            #vmin, vmax = -0.26, 0.06 #wpthlp astex noise comparison
            
            # ticks=None means matplotlib will generate ticks automatically
            # For non-correlation variables with positive and negative values, we want to create ticks manually
            # and individually for the positive and negative parts of the colorbar,
            # in case the orders of magnitude are different
            ticks = None
            if "corr" in self.dependent_title.lower() or "corr" in self.title.lower():
                # Variable is correlation data
                # -> normalize colormap to [-1;1] and use two-sided colormap
                vmin = -1
                vmax = 1
                cmap = mpl.colormaps.get_cmap(Style_definitions.CONTOUR_CMAP_CORR)
                norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
            else:
                # Variable is not a correlation
                if vmin<0 and vmax>0:
                    # Positive and negative values -> Use two-sided colormap
                    cmap = mpl.colormaps.get_cmap(Style_definitions.CONTOUR_CMAP_CORR)
                    cmap.set_over("green") # Color for values above vmax
                    cmap.set_under("yellow")  # Color for values below vmin
                    norm = mpl.colors.TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)
                    # Manually create ticks for the positive and negative part individually,
                    # then paste them together
                    ticks = np.append(np.linspace(vmin,0,5), np.linspace(0,vmax,6)[1:])
                else:
                    # Only one-sided values -> Use unicolor colormap
                    cmap = mpl.colormaps.get_cmap(Style_definitions.CONTOUR_CMAP_GENERAL)
                    cmap.set_over("green") # Color for values above vmax
                    cmap.set_under("black")  # Color for values below vmin
                    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
            cont_map = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
            label = var.label

            # Set graph size
            plt.figure(figsize=(10,6))

            # Prevent x-axis label from getting cut off
            # plt.gcf().subplots_adjust(bottom=0.15)

            plt.contourf(x_data, y_data, c_data.T, levels=cmap.N, norm=norm, cmap=cmap)
            plt.contour(x_data, y_data, c_data.T, levels=cmap.N, norm=norm, cmap=cmap)
            # Add colorbar to the current ax (gca)
            # For details on __getFormatter, see below
            #plt.gcf().colorbar(cont_map, ax=plt.gca(), cmap=cmap, norm=norm, label=self.dependent_title,
            #                   orientation='vertical', ticks=ticks, format=ContourPanel.__getFormatter())
            plt.gcf().colorbar(cont_map, ax=plt.gca(), cmap=cmap, norm=norm, label=self.dependent_title,
                               orientation='vertical', format=ContourPanel.__getFormatter())
            #plt.title(label + ' - ' + self.title, pad=10)
            plt.title(self.title + ' - ' + label, pad=10)
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

    @staticmethod
    def __getFormatter():
        # Factory for formatter of the colorbar ticks
        # The "__" means this is a private method that cannot be called from outside.
        #
        # If necessary, this factory could receive arguments and use them within the formatting function.
        # For example, if we wanted to set a limit for exponential notation depending on an input (e.g. the size of the figure):
        # Define __getFormatter(lim)
        # We could do the same split as below, and then check the exponent against lim
        #
        # This function should return a FuncFormatter (or any type of Formatter) that uses the following function for formatting of the ticks labels.
        # For more information, see here: https://www.geeksforgeeks.org/matplotlib-ticker-funcformatter-class-in-python/ (viewed 12/03/2024)
        # Or here: https://matplotlib.org/stable/gallery/ticks/tick-formatters.html (viewed 12/03/2024)
        def format_func(x,pos):
            # Use scientific notation for numbers with too many places (too small/too large)
            num, exp = '{:.3e}'.format(x).split('e')
            if abs(int(exp)) <=4:
                # Exponent small enough
                # -> regular floating point representation, strip all trailing zeros
                ret = '{:.5f}'.format(x).rstrip('0')
                if ret.endswith('.'):
                    # Print at least one decimal place
                    ret += '0'
                return ret
            else:
                # Very large and very small exponents
                # -> Use scientific representation with at most 3 decimal places
                # Strip trailing zeros
                num = num.rstrip('0')
                if num.endswith('.'):
                    # Print at least one decimal place
                    num += '0'
                return num+'e'+str(int(exp))
        # Wrap into Formatter and return
        return mpt.FuncFormatter(format_func)