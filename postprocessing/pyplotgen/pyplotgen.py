#!/usr/bin/env python3
"""
:Author: Nicolas Strike
:Date: Jan 2019

This file is the main driver for PyPlotGen. It handles the overall
process of retrieving the user's clubb and dictating the loading
of files, generation of plots, and then the export of plots, though
these processes are mostly carried out by other classes/files.
"""
import argparse
import glob
import os
import re
import shutil
import subprocess
from datetime import datetime
from warnings import warn

from config import Case_definitions
from python_html_gallery import gallery
from src.Case import Case
from src.DataReader import DataReader
from src.interoperability import clean_path


class PyPlotGen:
    """

    """

    def __init__(self, output_folder, clubb_folders=None, replace=False, les=False, cgbest=False, hoc=False,
                 plotrefs=False, benchmark_only=False,
                 zip=False, thin=False, no_legends=False, ensemble=False, plot_e3sm="", sam_folders=[""],
                 wrf_folders=[""], cam_folders=[""],
                 budget_moments=False, bu_morr=False, diff=None, show_alphabetic_id=False):
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
        self.clubb_folders = clubb_folders
        self.output_folder = output_folder
        self.replace_images = replace
        self.les = les
        self.e3sm_dir = plot_e3sm
        self.sam_folders = sam_folders
        self.cam_folders = cam_folders
        self.wrf_folders = wrf_folders
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
        self.cases_plotted = []
        self.clubb_datasets = None
        self.data_reader = DataReader()
        self.diff_files_data_reader = DataReader()
        self.sam_data_reader = DataReader()
        self.show_alphabetic_id = show_alphabetic_id
        self.output_folder = os.path.abspath(self.output_folder)
        self.benchmark_only = benchmark_only

        if os.path.isdir(self.output_folder) and self.replace_images == False:
            self.output_folder = self.output_folder + '_generated_on_' + str(datetime.now())

        self.output_folder = clean_path(self.output_folder)

    def run(self):
        """
        Runs PyPlotGen
        :return: None
        """
        self.clubb_datasets = self.data_reader.loadFolder(self.clubb_folders)
        diff_datasets = None
        if self.diff is not None:
            diff_datasets = self.diff_files_data_reader.loadFolder(self.diff)
        all_cases = Case_definitions.ALL_CASES

        # Downloads model output (sam, les, clubb) if it doesn't exist
        if self.__benchmark_files_needed__():
            self.__downloadModelOutputs__()

        if self.replace_images:
            print('###########################################')
            print("\nDeleting old plots")
            subprocess.run(['rm', '-rf', self.output_folder + '/'])
        num_cases_plotted = 0
        for case_def in all_cases:
            if self.__dataForCaseExists__(case_def):
                num_cases_plotted += 1
                print('###########################################')
                print("plotting ", case_def['name'])
                case_diff_datasets = None
                casename = case_def['name']
                if self.diff is not None:
                    case_diff_datasets = diff_datasets[casename]
                if casename in self.clubb_datasets.keys():
                    clubb_case_datasets = self.clubb_datasets[casename]
                else:
                    clubb_case_datasets = {}
                case = Case(case_def, clubb_folders=clubb_case_datasets, plot_les=self.les,
                            plot_budgets=self.plot_budgets, sam_folders=self.sam_folders, wrf_folders=self.wrf_folders,
                            diff_datasets=case_diff_datasets, plot_r408=self.cgbest, plot_hoc=self.hoc,
                            e3sm_dirs=self.e3sm_dir, cam_folders = self.cam_folders)
                case.plot(self.output_folder, replace_images=self.replace_images, no_legends=self.no_legends,
                          thin_lines=self.thin, show_alphabetic_id=self.show_alphabetic_id)
                self.cases_plotted.append(casename)
        print('###########################################')
        if num_cases_plotted == 0:
            warn(
                "Warning, no cases were plotted! Please either specify an input folder for a supported model "
                "(e.g. using --sam, --clubb, or --e3sm) or make sure the "
                "default clubb output folder contains .nc output. "
                "Please run ./pyplotgen.py -h for more information on parameters.")
        print("\nGenerating webpage for viewing plots ")
        if not os.path.exists(self.output_folder):
            os.mkdir(self.output_folder)

        self.__copySetupFiles__()
        gallery.main(self.output_folder)
        print('###########################################')
        print("Output can be viewed at file://" + self.output_folder + "/index.html with a web browser")

    def __dataForCaseExists__(self, case_def):
        """
        Returns true if there's an clubb nc file for a given case name.

        :param no_clubb: True/false for whether or not clubb lines are to be plotted
        :param case_def: The case definition object
        :param all_case_names: List of all case names that can be plotted
        :return:
        """

        e3sm_given = len(self.e3sm_dir) != 0
        sam_given = len(self.sam_folders) != 0
        wrf_given = len(self.wrf_folders) != 0
        cam_given = len(self.cam_folders) != 0
        clubb_given = self.clubb_datasets is not None

        e3sm_has_case = e3sm_given and case_def['e3sm_file'] != None
        sam_has_case = sam_given and case_def['sam_file'] != None
        wrf_has_case = wrf_given != 0 and case_def['wrf_file'] != None
        cam_has_case = cam_given and case_def['cam_file'] != None
        clubb_has_case = clubb_given and case_def['name'] in self.clubb_datasets.keys()

        return e3sm_has_case or sam_has_case or cam_has_case or wrf_has_case or clubb_has_case or self.benchmark_only

    def __copySetupFiles__(self):
        """
        Copies case setup files from the clubb folder(s)
        into the pyplotgen output folders.

        :return:
        """
        print("Looking for case_setup.txt files")
        for folder in self.clubb_folders:
            setup_file_search_pattern = folder + '/*_setup.txt'
            folder_basename = os.path.basename(folder)

            for file in glob.glob(setup_file_search_pattern):
                file_basename = os.path.basename(file)
                casename = file_basename[:-10]
                copy_dest_folder = self.output_folder + '/' + casename + '/'
                copy_dest_file = copy_dest_folder + file_basename[:-10] + '_' + folder_basename + file_basename[-10:]

                if casename in self.cases_plotted:
                    if not os.path.exists(copy_dest_folder):
                        os.mkdir(copy_dest_folder)
                    shutil.copy(file, copy_dest_file)
                    print("\tFound setup file " + str(file))

    def __benchmark_files_needed__(self):
        """
        Returns true if the user requested to plot model output
        that needs to be downloaded
        :return:
        """
        return self.les or self.hoc or self.cgbest

    def __downloadModelOutputs__(self):
        """
        Checks for model output, e.g. sam benchmark runs, and if it
        doesn't exist, downloads it.
        Supports: SAM, COAMPS, CLUBB

        :return:
        """

        # Ensure benchmark output is available
        print("Checking for model benchmark output...")
        if not os.path.isdir(Case_definitions.BENCHMARK_OUTPUT_ROOT) and \
                not os.path.islink(Case_definitions.BENCHMARK_OUTPUT_ROOT):
            print("\tDownloading the benchmarks to " + Case_definitions.BENCHMARK_OUTPUT_ROOT)
            subprocess.run(['git', 'clone', 'https://carson.math.uwm.edu/les_and_clubb_benchmark_runs.git',
                            Case_definitions.BENCHMARK_OUTPUT_ROOT])
        else:
            print("Benchmark output found in " + Case_definitions.BENCHMARK_OUTPUT_ROOT)

    def __remove_invalid_filename_chars__(self, string):
        """
        Removes characters from a string that are not valid for a filename
        DEPRECATED! Moved to src/interoperability.py
        since this is needed mutiple times throughout a run

        :param string: Filename string to have characters removed
        :return: a character stripped version of the filename
        """
        string = string.replace('.', '')
        string = string.replace(',', '')
        if 'win' in os.name.lower() or 'nt' in os.name.lower():
            string = re.sub(r'(?<![A-Z]):(?!\\)', '-', string)
        else:
            string = string.replace(':', '-')
        return string

def __process_args__():
    """
    This method takes arguments in from the command line and feeds them into
    a PyPlotGen object

    :return: a PyPlotGen object containing the parameters as given from the commandline.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--replace", help="If the output folder already exists, replace it with the new one.",
                        action="store_true")
    parser.add_argument("-l", "--les", help="Plot LES dependent_data for comparison.", action="store_true")
    parser.add_argument("-e", "--e3sm",
                        help="Plot E3SM dependent_data for comparison. Pass a folder in with this option. This "
                             "folder must contain the e3sm nc files in the root directory, "
                             "where each filename is the name of the clubb case to plot it with. E.g. "
                             "name a file dycoms2_rfo2_ds.nc to plot it with that clubb case.",
                        action="store",
                        default=[], nargs='+')
    parser.add_argument("-g", "--plot-golaz-best", help="Plot Chris Golaz Best Ever dependent_data for comparison.",
                        action="store_true")
    parser.add_argument("-d", "--plot-hoc-2005", help="Plot !HOC 12/17/2015 dependent_data for comparison.",
                        action="store_true")
    parser.add_argument("-a", "--all-best",
                        help="Same as -lbd. Plots LES, Golaz Best Ever, and HOC 2005 dependent_data for comparison.",
                        action="store_true")
    parser.add_argument("-z", "--zip", help="Output dependent_data into a compressed zip file.", action="store_true")
    parser.add_argument("--show-alphabetic-id", help="Add an identifying character to the top right of a panel.",
                        action="store_true")
    parser.add_argument("--thin", help="Plot using thin solid lines.", action="store_true")
    parser.add_argument("--no-legends", help="Plot without legend boxes defining the line types.", action="store_true")
    parser.add_argument("--ensemble", help="Plot ensemble tuner runs", action="store_true")  # TODO is this needed?
    parser.add_argument("-b", "--plot-budgets", help="Plot all defined budgets of moments", action="store_true")
    parser.add_argument("--bu-morr",
                        help="For morrison microphysics: breaks microphysical source terms into component processes",
                        action="store_true")
    parser.add_argument("--benchmark-only",
                        help="Prevents autoplotting of clubb's default output folder when no input folders are "
                             "specifed. This results in only plotting the benchmark output, though this "
                             "output doesn't gurantee all text fields or plots are filled.",
                        action="store_true")
    parser.add_argument("--diff", help="Plot the difference between two clubb folders", action="store")
    parser.add_argument("-c", "--clubb", help="Input folder(s) containing clubb netcdf data.", action="store",
                        default=[], nargs='+')
    parser.add_argument("-s", "--sam", help="Input folder(s) containing sam netcdf data.", action="store",
                        default=[], nargs='+')
    parser.add_argument("--cam", help="Input folder(s) containing cam netcdf data.", action="store",
                        default=[], nargs='+')
    parser.add_argument("-w", "--wrf", help="Input folder(s) containing wrf netcdf data.", action="store",
                        default=[], nargs='+')
    parser.add_argument("-o", "--output", help="Name of folder to create and store plots into.", action="store",
                        default="./output")
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
    e3sm = args.e3sm

    # If the last char in folder path is /, remove it
    for i in range(len(args.clubb)):
        if args.clubb[i][-1] == "/":
            args.clubb[i] = args.clubb[i][:-1]

    for i in range(len(args.sam)):
        if args.sam[i][-1] == "/":
            args.sam[i] = args.sam[i][:-1]

    for i in range(len(args.cam)):
        if args.cam[i][-1] == "/":
            args.cam[i] = args.cam[i][:-1]

    for i in range(len(args.e3sm)):
        if args.e3sm[i][-1] == "/":
            args.e3sm[i] = args.e3sm[i][:-1]

    for i in range(len(args.wrf)):
        if args.wrf[i][-1] == "/":
            args.wrf[i] = args.e3sm[i][:-1]

    no_folders_inputed = len(args.e3sm) == 0 and len(args.sam) == 0 and len(args.clubb) == 0 \
                         and len(args.wrf) == 0 and len(args.cam) == 0

    if no_folders_inputed and not args.benchmark_only:
        args.clubb = ["../../output"]

    if args.all_best:
        les = True
        cgbest = True
        hoc = True

    pyplotgen = PyPlotGen(args.output, clubb_folders=args.clubb, replace=args.replace, les=les, plot_e3sm=e3sm,
                          cgbest=cgbest, cam_folders = args.cam,
                          hoc=hoc, plotrefs=args.all_best, zip=args.zip, thin=args.thin, sam_folders=args.sam,
                          wrf_folders=args.wrf, benchmark_only=args.benchmark_only,
                          no_legends=args.no_legends, ensemble=args.ensemble, budget_moments=args.plot_budgets,
                          bu_morr=args.bu_morr, diff=args.diff, show_alphabetic_id=args.show_alphabetic_id)
    return pyplotgen


if __name__ == "__main__":
    pyplotgen = __process_args__()
    pyplotgen.run()
