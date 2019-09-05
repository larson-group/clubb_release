#!/usr/bin/env python3
"""
:Author: Nicolas Strike
:Date: Jan 2019

This file is the main driver for PyPlotGen. It handles the overall
process of retrieving the user's input and dictating the loading
of files, generation of plots, and then the export of plots, though
these processes are mostly carried out by other classes/files.
"""
import argparse
import os
import subprocess
from warnings import warn

from config import Case_definitions
from python_html_gallery import gallery
from src.Case import Case
from src.DataReader import DataReader


class PyPlotGen:
    """

    """

    def __init__(self, input_folder, output_folder, replace=False, les=False, cgbest=False, hoc=False, plotrefs=False,
                zip=False, thin=False, no_legends=False, ensemble=False,
                 budget_moments=False, bu_morr=False, diff=None):
        """
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
        """
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
        if self.output_folder[0] == '.':
            self.output_folder = os.path.dirname(os.path.realpath(__file__)) + "/" + self.output_folder

    def run(self):
        """
        Runs PyPlotGen
        :return: None
        """
        self.nc_datasets = self.data_reader.loadFolder(self.input_folder)
        diff_datasets = None
        if self.diff is not None:
            diff_datasets = self.diff_files_data_reader.loadFolder(self.diff)
        all_cases = Case_definitions.ALL_CASES

        # Downloads model output (sam, les, clubb) if it doesn't exist
        self.__downloadModelOutputs__()

        # TODO Handle dataset not found/partial nc output
        if self.replace_images:
            print('###########################################')
            print("\nDeleting old plots")
            subprocess.run(['rm', '-rf', self.output_folder + '/'])
        cases_plotted = 0
        for case_def in all_cases:
            if case_def['name'] in self.nc_datasets.keys():
                cases_plotted += 1
                print('###########################################')
                print("plotting ", case_def['name'])
                case_diff_datasets = None
                if self.diff is not None:
                    case_diff_datasets = diff_datasets[case_def['name']]
                case = Case(case_def, self.nc_datasets[case_def['name']], plot_les=self.les, plot_budgets=self.plot_budgets,
                            diff_datasets=case_diff_datasets, plot_r408=self.cgbest, plot_hoc=self.hoc)
                case.plot(self.output_folder, replace_images=self.replace_images, no_legends = self.no_legends, thin_lines=self.thin)
            # else:
            #     print("Clubb input nc files not found. Skipping case.")
        print('###########################################')
        if cases_plotted == 0:
            warn("Warning, no cases were plotted! Make sure the input folder " + self.input_folder +
                 " contains .nc files or use the --input (-i) parameter to manually specify an input folder.")
        print("\nGenerating webpage for viewing plots ")
        if not os.path.exists(self.output_folder):
            os.mkdir(self.output_folder)

        gallery.main(self.output_folder)
        print('###########################################')
        print("Output can be viewed at file://" + self.output_folder + "/index.html with a web browser")

    def __downloadModelOutputs__(self):
        """
        Checks for model output, e.g. sam benchmark runs, and if it
        doesn't exist, downloads it.
        Supports: SAM, COAMPS, CLUBB

        :return:
        """


        # Ensure benchmark output is available
        print("Checking for model benchmark output...")
        if not os.path.isfile(Case_definitions.BENCHMARK_OUTPUT_ROOT) and not os.path.islink(
                Case_definitions.BENCHMARK_OUTPUT_ROOT):
            print("Benchmark output was not found in " + Case_definitions.BENCHMARK_OUTPUT_ROOT + ", downloading now.")
            subprocess.run(['git', 'clone', 'https://carson.math.uwm.edu/les_and_clubb_benchmark_runs.git'])
        else:
            print("Benchmark output found in " + Case_definitions.BENCHMARK_OUTPUT_ROOT)

def __process_args__():
    """
    This method takes arguments in from the command line and feeds them into
    a PyPlotGen object

    :return: a PyPlotGen object containing the parameters as given from the commandline.
    """
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
    parser.add_argument("-i", "--input", help="Input folder containing netcdf output data.", action="store", default="../../output")
    parser.add_argument("-o", "--output", help="Name of folder to create and store plots into.", action="store", default="./output")
    args = parser.parse_args()

    if args.zip:
        print("Zip flag detected, but that feature is not yet implemented")
    if args.ensemble:
        print("Ensemble flag detected, but that feature is not yet implemented")
    if args.bu_morr:
        print("Morrison breakdown flag detected, but that feature is not yet implemented")
    les = args.les
    cgbest = args.plot_golaz_best
    hoc = args.plot_hoc_2005

    if args.all_best:
        les = True
        cgbest = True
        hoc = True

    pyplotgen = PyPlotGen(args.input, args.output, replace=args.replace, les=les, cgbest=cgbest,
                          hoc=hoc, plotrefs=args.all_best, zip=args.zip, thin=args.thin,
                          no_legends=args.no_legends, ensemble=args.ensemble, budget_moments=args.plot_budgets,
                          bu_morr=args.bu_morr, diff=args.diff)
    return pyplotgen

if __name__ == "__main__":
    pyplotgen = __process_args__()
    pyplotgen.run()
