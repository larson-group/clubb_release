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
import subprocess

import Case_definitions
from Case import Case
from DataReader import DataReader


class PyPlotGen:
    '''

    '''

    def __init__(self, input_folder, output_folder, replace=False, les=False, cgbest=False, hoc=False, plotrefs=False,
                zip=False, thin=False, no_legends=False, ensemble=False,
                 budget_moments=False, bu_morr=False, diff=None):
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
        self.replace_images = replace
        self.les = les
        self.cgbest = cgbest
        self.hoc = hoc
        self.plot_diff = diff
        self.plotrefs = plotrefs
        self.zip = zip
        self.thin = thin
        self.no_legends = no_legends
        self.ensemble = ensemble
        self.plot_budgets = budget_moments
        self.bu_morr = bu_morr
        self.diff = diff
        self.nc_datasets = None
        self.data_reader = DataReader()
        self.diff_files_data_reader = DataReader()
        self.sam_data_reader = DataReader()

    def run(self):
        '''
        Runs all PyPlotGen cases
        :return: n/a
        '''
        self.nc_datasets = self.data_reader.loadFolder(self.input_folder)
        diff_datasets = None
        if self.diff is not None:
            diff_datasets = self.diff_files_data_reader.loadFolder(self.diff)
        all_cases = Case_definitions.ALL_CASES

        # Ensure SAM output is available
        print("Checking for SAM output...")
        if not os.path.isfile(Case_definitions.SAM_OUTPUT_ROOT) and not os.path.islink(Case_definitions.SAM_OUTPUT_ROOT):
            print("Sam output was not found in " + Case_definitions.SAM_OUTPUT_ROOT + ", downloading now.")
            subprocess.run(['git', 'clone', 'https://carson.math.uwm.edu/sam_benchmark_runs.git'])  # Use sigal to build html in '_build/'
        else:
            print("Sam output found in " + Case_definitions.SAM_OUTPUT_ROOT)
        # TODO Handle dataset not found/partial nc output
        for case_def in all_cases:
            print('###########################################')
            print("plotting ", case_def['name'])
            case_diff_datasets = None
            if self.diff is not None:
                case_diff_datasets = diff_datasets[case_def['name']]
            case = Case(case_def,self.nc_datasets[case_def['name']], plot_sam=self.les, plot_budgets=self.plot_budgets, diff_datasets=case_diff_datasets)
            case.plot(self.output_folder, replace_images=self.replace_images, no_legends = self.no_legends, thin_lines=self.thin)

        print('###########################################')
        print("\nGenerating webpage for viewing plots ")
        subprocess.run(['sigal', 'build', '-f', self.output_folder + '/'])  # Use sigal to build html in '_build/'
        print('###########################################')
        print("Output can be viewed at file://" + self.output_folder + "/../_build/index.html with a web browser")
def process_args():
    '''
    This method takes arguments in from the command line and feeds them into
    a PyPlotGen object

    :return: a PyPlotGen object containing the parameters as given from the commandline.
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--replace", help="If the output folder already exists, replace it with the new one.", action="store_true")
    parser.add_argument("-l", "--les", help="Plot LES data for comparison.", action="store_true")
    parser.add_argument("-b", "--plot-golaz-best", help="Plot Chris Golaz Best Ever data for comparison.", action="store_true")
    parser.add_argument("-d", "--plot-hoc-2005", help="Plot !HOC 12/17/2015 data for comparison.", action="store_true")
    parser.add_argument("-a", "--all-best", help="Same as -lbd. Plots LES, Golaz Best Ever, and HOC 2005 data for comparison.", action="store_true")
    parser.add_argument("-z", "--zip", help="Output data into a compressed zip file.", action="store_true")
    parser.add_argument("--thin", help="Plot using thin solid lines.", action="store_true")
    parser.add_argument("--no-legends", help="Plot without legend boxes defining the line types.", action="store_true")
    parser.add_argument("--ensemble", help="Plot ensemble tuner runs", action="store_true")  # TODO is this needed?
    parser.add_argument("--plot-budgets", help="Plot all defined budgets of moments", action="store_true")
    parser.add_argument("--bu-morr", help="For morrison microphysics: breaks microphysical source terms into component processes",action="store_true")
    parser.add_argument("--diff", help="Plot the difference between two input folders", action="store")
    parser.add_argument("input", help="Input folder containing netcdf output data.", action="store")
    parser.add_argument("output", help="Name of folder to create and store plots into.", action="store")
    args = parser.parse_args()

    if args.plot_golaz_best:
        print("Plot golaz best flag detected, but that feature is not yet implemented")
    if args.plot_hoc_2005:
        print("Plot HOC 2005 flag detected, but that feature is not yet implemented")
    if args.all_best:
        print("Plot all reference plots flag detected, but that feature is not yet implemented")
    if args.zip:
        print("Zip flag detected, but that feature is not yet implemented")
    if args.ensemble:
        print("Ensemble flag detected, but that feature is not yet implemented")
    if args.bu_morr:
        print("Morrison breakdown flag detected, but that feature is not yet implemented")

    pyplotgen = PyPlotGen(args.input, args.output, replace=args.replace, les=args.les, cgbest=args.plot_golaz_best,
                          hoc=args.plot_hoc_2005, plotrefs=args.all_best, zip=args.zip, thin=args.thin,
                          no_legends=args.no_legends, ensemble=args.ensemble, budget_moments=args.plot_budgets,
                          bu_morr=args.bu_morr, diff=args.diff)
    return pyplotgen

if __name__ == "__main__":
    pyplotgen = process_args()
    pyplotgen.run()
