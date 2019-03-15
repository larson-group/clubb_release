'''
This file is the main driver for PyPlotGen. It handles the overall
process of retrieving the user's input and dictating the loading
of files, generation of plots, and then the export of plots, though
these processes are mostly carried out by other classes/files.

Author: Nicolas Strike
Date: Jan 2019
'''
import argparse
import os

import configparser
import pathlib

from DataReader import DataReader
from Plotter import Plotter


class PyPlotGen:
    def __init__(self, input_folder, output_folder, replace=False, les=False, cgbest=False, hoc=False,plotrefs=False,
                 nightly=False,hq_imgs=False, eps=False, zip=False, thin=False, no_legends=False,ensemble=False,
                 budget_moments=False, bu_morr=False,diff=False):
        '''
        This creates an instance of PyPlotGen. Each parameter is a command line parameter passed in from the argparser
        below.
        :param input_folder:
        :param output_folder:
        :param replace:
        :param les:
        :param cgbest:
        :param hoc:
        :param plotrefs:
        :param nightly:
        :param hq_imgs:
        :param eps:
        :param zip:
        :param thin:
        :param no_legends:
        :param ensemble:
        :param budget_moments:
        :param bu_morr:
        :param diff:
        '''
        self.input_folder = input_folder
        self.output_folder = output_folder
        self.replace = replace
        self.les = les
        self.cgbest = cgbest
        self.hoc = hoc
        self.plotrefs = plotrefs
        self.nightly = nightly
        self.hq_imgs = hq_imgs
        self.eps = eps
        self.zip = zip
        self.thin = thin
        self.no_legends = no_legends
        self.ensemble = ensemble
        self.budget_moments = budget_moments
        self.bu_morr = bu_morr
        self.diff = diff
        self.nc_datasets = None
        self.data_reader = DataReader()


    def run(self):
        '''
        Runs PyPlotGen
        :return: n/a
        '''
        self.nc_datasets = self.data_reader.loadFolder(self.input_folder)
        for dataset in self.nc_datasets:
            dataset_filename = dataset.filepath()
            if ("dycoms_rf01_96x96x320") in dataset_filename.lower():
                casefile = "/home/strike/clubb/postprocessing/pyplotgen/cases/Case_dycoms2_rf01.ini"
                print("Opening " + dataset_filename)
                self.plotCase(dataset, casefile)
                break # TODO TEMP FIX


    def plotCase(self, dataset, casefile):
        '''
        This method plots all the panels relevant to a case
        to png files

        :param dataset: The NetCDF dataset containing the desired data
        :param casefile: The full filename of a casefile describing what panels to plot
        :return: none
        '''
        case_config = configparser.ConfigParser()
        case_config.read(casefile)

        all_plot_values = [] # TODO change None to be the selected case
        all_plot_details = []
        plotter = Plotter()

        panel_groups_str = case_config['defaults']['panel_groups']
        panel_groups = self.data_reader.getArrayFromString(panel_groups_str)
        root_dir = pathlib.Path(__file__).parent
        panels_dir = str(root_dir) + "/cases/panels/"
        for panel_group in panel_groups:
            panel_group_dir =  panels_dir + panel_group + "/"
            panels = self.getFilesInFolder(panel_group_dir)
            for panel in panels:
                panel_config = configparser.ConfigParser()
                panel_config.read(panel)

                title = panel_config['defaults']['title'] #self.get_long_name(netcdf_data, x_variable_name)
                x_axis_title = panel_config['defaults']['x_label'] # x_variable_name + "[K]"
                y_axis_title = panel_config['defaults']['y_label']
                plot_details = Plotter.PlotDetails(title=title, y_title=y_axis_title, x_title=x_axis_title)
                all_plot_details.append(plot_details)

                plot_values = self.data_reader.getPlotsData(dataset, casefile, panel)
                all_plot_values.append(plot_values)

        for i in range(len(all_plot_values)):
            plotter.plot(all_plot_values[i],all_plot_details[i]) # TODO this is a temp test value, loop below is better


    def getFilesInFolder(self, folder):
            '''
            Returns a list of files contained in a folder

            :param folder: The folder containing files
            :return: list of files in the given folder
            '''
            panel_files = []
            for root, dirs, files in os.walk(folder):
                for filename in files:
                    abs_filename = os.path.abspath(os.path.join(root, filename))
                    panel_files.append(abs_filename)

            return panel_files

def process_args():
    '''
    This method takes arguments in from the command line and feeds them into
    a PyPlotGen object

    :return: a PyPlotGen object containing the parameters as given from the commandline.
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument("-r","--replace",help="If the output folder already exists, replace it with the new one.", action="store_true")
    parser.add_argument("-l","--les",help="Plot LES data for comparison.", action="store_true")
    parser.add_argument("-b","--plot-golaz-best",help="Plot Chris Golaz Best Ever data for comparison.", action="store_true")
    parser.add_argument("-d","--plot-hoc-2005",help="Plot !HOC 12/17/2015 data for comparison.", action="store_true")
    parser.add_argument("-a","--all-best",help="Same as -lbd. Plots LES, Golaz Best Ever, and HOC 2005 data for comparison.", action="store_true")
    parser.add_argument("-n","--nightly",help="Run in nightly mode.", action="store_true") # TODO update this description with what it means
    parser.add_argument("-q","--high-quality-imgs",help="Output high quality images (does not auto scale)", action="store_true") #TODO update if this does auto scale
    parser.add_argument("-e","--save-eps",help="Does not delete EPS images after conversion.", action="store_true") #TODO what does this mean?
    parser.add_argument("-z","--zip",help="Output data into a compressed zip file.", action="store_true")
    parser.add_argument("--thin",help="Plot using thin solid lines.", action="store_true")
    parser.add_argument("--no-legends",help="Plot without legend boxes defining the line types.", action="store_true")
    parser.add_argument("--ensemble",help="Plot ensemble tuner runs", action="store_true") #TODO is this needed?
    parser.add_argument("--budget-moments",help="Plot all defined budgets of moments", action="store_true")
    parser.add_argument("--bu-morr",help="For morrison microphysics: breaks microphysical source terms into component processes", action="store_true")
    parser.add_argument("--diff",help="Plot the difference between two input folders", action="store_true")
    parser.add_argument("input",help="Input folder containing grads or netcdf output data.", action="store")
    parser.add_argument("output",help="Name of folder to create and store plots into.", action="store")
    args = parser.parse_args()

    if args.replace:
        print("Replace flag detected, but that feature is not yet implemented")
    if args.les:
        print("LES flag detected, but that feature is not yet implemented")
    if args.plot_golaz_best:
        print("Plot golaz best flag detected, but that feature is not yet implemented")
    if args.plot_hoc_2005:
        print("Plot HOC 2005 flag detected, but that feature is not yet implemented")
    if args.all_best:
        print("Plot all reference plots flag detected, but that feature is not yet implemented")
    if args.nightly:
        print("Nightly flag detected, but that feature is not yet implemented")
    if args.high_quality_imgs:
        print("High quality images flag detected, but that feature is not yet implemented")
    if args.save_eps:
        print("Save EPS images flag detected, but that feature is not yet implemented")
    if args.zip:
        print("Zip flag detected, but that feature is not yet implemented")
    if args.thin:
        print("Thin plot lines flag detected, but that feature is not yet implemented")
    if args.no_legends:
        print("No legends flag detected, but that feature is not yet implemented")
    if args.ensemble:
        print("Ensemble flag detected, but that feature is not yet implemented")
    if args.budget_moments:
        print("Budget moments flag detected, but that feature is not yet implemented")
    if args.bu_morr:
        print("Morrison breakdown flag detected, but that feature is not yet implemented")
    if args.diff:
        print("Diff flag detected, but that feature is not yet implemented")

    pyplotgen = PyPlotGen(args.input, args.output, replace=args.replace, les=args.les, cgbest=args.plot_golaz_best,
                          hoc=args.plot_hoc_2005,plotrefs=args.all_best, nightly=args.nightly,
                          hq_imgs=args.high_quality_imgs, eps=args.save_eps, zip=args.zip, thin=args.thin,
                          no_legends=args.no_legends,ensemble=args.ensemble, budget_moments=args.budget_moments,
                          bu_morr=args.bu_morr,diff=args.diff)
    return pyplotgen
# run main if file is run by interpreter
if __name__ == "__main__":
    pyplotgen = process_args()
    pyplotgen.run()
