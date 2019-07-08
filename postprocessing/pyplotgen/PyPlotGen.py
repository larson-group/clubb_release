'''
This file is the main driver for PyPlotGen. It handles the overall
process of retrieving the user's input and dictating the loading
of files, generation of plots, and then the export of plots, though
these processes are mostly carried out by other classes/files.

Author: Nicolas Strike
Date: Jan 2019
'''
import argparse
import subprocess
import sys

from pyplotgen.DataReader import DataReader
from pyplotgen.CaseTest import CaseTest

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
        # self.output_folder = output_folder
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
        self.sam_data_reader = DataReader()


    def run(self):
        '''
        Runs PyPlotGen
        :return: n/a
        '''
        self.nc_datasets = self.data_reader.loadFolder(self.input_folder)
        # self.sam_datasets = self.sam_data_reader.loadFolder("/home/strike/sam_benchmark_runs")
        cases_plotted = {}
        for case_key in self.nc_datasets.keys():
            ncdf_files = self.nc_datasets[case_key]
            dataset_filenames = [dataset.filepath() for dataset in self.nc_datasets[case_key].values()]

            if "dycoms2_rf01_" in dataset_filenames[0] and "sst" not in dataset_filenames[0] and ".git" not in dataset_filenames:
                sam_file = self.sam_data_reader.__loadNcFile__("/home/nicolas/sam_benchmark_runs/DYCOMS_RF01_96x96x320.nc")
                # sam_file = None
                print("Plotting case ", case_key)
                test_case = CaseTest(ncdf_files, sam_file=sam_file)
                test_case.plot()
                cases_plotted[test_case.name] = test_case
                break # TODO TEMP FIX
        subprocess.run(['sigal', 'build', '-f', 'output/'])  # Use sigal to build html in '_build/'

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
    parser.add_argument("input",help="Input folder containing netcdf output data.", action="store")
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
        print("Thin plotAll lines flag detected, but that feature is not yet implemented")
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

if __name__ == "__main__":
    pyplotgen = process_args()
    pyplotgen.run()
