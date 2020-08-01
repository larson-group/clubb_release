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
    Main class for the pyplotgen program.
    After processing the command line options in the main function,
    an instance of this class is created in __process_args__, passing all option parameters.
    Finally the run function of this class is the actual driver of the pyplotgen program,
    creating the Case objects and calling their plot functions to generate the images and html output.
    """

    def __init__(self, output_folder, clubb_folders=None, replace=False, les=False, cgbest=False, hoc=False,
                 benchmark_only=False, nightly=False, zip=False, thin=False, no_legends=False, ensemble=False,
                 plot_e3sm="", sam_folders=[""], wrf_folders=[""], cam_folders=[""],
                 budget_moments=False, bu_morr=False, diff=None, show_alphabetic_id=False,
                 time_height=False, animation=None):
        """
        This creates an instance of PyPlotGen. Each parameter is a command line parameter passed in from the argparser
        below.
        
        :param output_folder: String containing the foldername into which the output files should be put.
        :param clubb_folders: List of foldernames containing CLUBB netcdf files to be plotted.
            Pass a list of folders in with this option.
            These folders must contain the CLUBB nc files in the root directory,
            where each filename is the name of the case to be plotted with the CLUBB-specific suffixes: _zm, _zt
            E.g. name files dycoms2_rfo2_ds_zm.nc and dycoms2_rfo2_ds_zt.nc
            to plot the file's data with that case's setup.
        :param replace: If False, an already existing outout folder will not be overwritten.
        :param les: If True, plot LES dependent_data for comparison.
        :param cgbest: If True, plot Chris Golaz Best Ever dependent_data for comparison.
        :param hoc: If True, plot !HOC 12/17/2015 dependent_data for comparison.
        :param benchmark_only: If True, autoplotting of CLUBB's default output folder is shut off.
            So only benchmark output will be plotted.
        :param nightly: Apply special parameters only relevant when running as part of a nightly test.
            This is currently limited to disabling case output if not all models have data for a given case.
            E.g. this prevents wrf plots from including cases that only have clubb plots and no wrf plots.
            Do not plot this with clubb-only plots, just plot clubb normally for clubb nightly tests.
        :param zip: If True, output dependent_data into a compressed zip file. Not implemented.
        :param thin: If True, plot using thin solid lines.
        :param no_legends: If True, plots will not have legend boxes listing the line types.
        :param ensemble: If True, plot ensemble tuner runs. Not implemented.
        :param plot_e3sm: Plot E3SM dependent_data for comparison.
            This parameter works exactly like the clubb_folders parameter, except E3SM uses only one nc file.
        :param sam_folders: Plot SAM dependent_data for comparison.
            This parameter works exactly like the clubb_folders parameter, except SAM uses only one nc file.
        :param wrf_folders: Plot WRF dependent_data for comparison.
            This parameter works exactly like the clubb_folders parameter, except WRF uses the following nc files:
            _zm, _zt, and _sfc
        :param cam_folders: Plot CAM dependent_data for comparison.
            This parameter works exactly like the clubb_folders parameter, except CAM uses only one nc file. (CHECK!)
        :param budget_moments: If True, plot all defined budgets of moments.
        :param bu_morr: For morrison microphysics: If True, break microphysical source terms into component processes.
            Not implemented.
        :param diff: Plot the difference between two clubb folders. (MORE DESCRIPTION)
        :param show_alphabetic_id: If True, add an identifying character to the top right of a panel.
        :param time_height: If True, plot time-height (contourf) plots instead of profile-like plots
        :param animation: If True, create time animations instead of time-averaged plots
            (works with profile and budget plots) (Not yet implemented).
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
        self.nightly = nightly
        self.time_height = time_height
        self.animation = animation

        if os.path.isdir(self.output_folder) and self.replace_images == False:
            self.output_folder = self.output_folder + '_generated_on_' + str(datetime.now())

        self.output_folder = clean_path(self.output_folder)

    def run(self):
        """
        Main driver of the pyplotgen program executing the following steps:
        - Download benchmark files if needed
        - Loop through cases listed in Case_definitions.ALL_CASES
        - Create instances of these cases
        - Call plot function on the instances
        - Create output folder
        - Generate html page containing plots
        
        :return: None
        """
        diff_datasets = None
        # Load data used for difference plots
        if self.diff is not None:
            diff_datasets = self.diff_files_data_reader.loadFolder(self.diff)
        all_cases = Case_definitions.ALL_CASES

        # Downloads model output (sam, les, clubb) if it doesn't exist
        if self.__benchmarkFilesNeeded__():
            self.__downloadModelOutputs__()

        # If --replace flag was set, delete old output folder
        if self.replace_images:
            print('###########################################')
            print("\nDeleting old plots")
            subprocess.run(['rm', '-rf', self.output_folder + '/'])
            # TODO: Use for Windows
            # shutil.rmtree(self.output_folder)
        num_cases_plotted = 0
        # Loop through cases listed in Case_definitions.ALL_CASES
        for case_def in all_cases:
            if self.__dataForCaseExists__(case_def):
                num_cases_plotted += 1
                print('###########################################')
                print("plotting ", case_def['name'])
                case_diff_datasets = None
                casename = case_def['name']
                if self.diff is not None:
                    case_diff_datasets = diff_datasets[casename]
                case = Case(case_def, clubb_folders=self.clubb_folders, plot_les=self.les,
                            plot_budgets=self.plot_budgets, sam_folders=self.sam_folders, wrf_folders=self.wrf_folders,
                            diff_datasets=case_diff_datasets, plot_r408=self.cgbest, plot_hoc=self.hoc,
                            e3sm_dirs=self.e3sm_dir, cam_folders=self.cam_folders,
                            time_height=self.time_height, animation=self.animation)
                # Call plot function of case instance
                case.plot(self.output_folder, replace_images=self.replace_images, no_legends=self.no_legends,
                          thin_lines=self.thin, show_alphabetic_id=self.show_alphabetic_id)
                self.cases_plotted.append(casename)
        print('###########################################')
        if num_cases_plotted == 0:
            warn(
                "Warning, no cases were plotted! Please either specify an input folder for a supported model "
                "(e.g. using --sam, --clubb, --e3sm, --wrf, --cam) or make sure the "
                "default clubb output folder contains .nc output. "
                "Please run ./pyplotgen.py -h for more information on parameters.")
        print("\nGenerating webpage for viewing plots ")

        if not os.path.exists(self.output_folder):
            os.mkdir(self.output_folder)

        self.__copySetupFiles__()
        # Generate html pages
        gallery.main(self.output_folder)
        print('###########################################')
        print("Output can be viewed at file://" + self.output_folder + "/index.html with a web browser")

    def __dataForCaseExists__(self, case_def):
        """
        Returns true if there's an nc file for a given case name and any model that should be plotted.

        :param case_def: The case definition object
        :return: True if data exists, False if not
        """
        e3sm_given = len(self.e3sm_dir) != 0
        sam_given = len(self.sam_folders) != 0
        wrf_given = len(self.wrf_folders) != 0
        cam_given = len(self.cam_folders) != 0
        clubb_given = len(self.clubb_folders) != 0

        e3sm_file_defined = case_def['e3sm_file'] is not None
        sam_file_defined = case_def['sam_file'] is not None
        wrf_file_defined = case_def['wrf_file'] is not None
        cam_file_defined = case_def['cam_file'] is not None
        clubb_file_defined = case_def['clubb_file'] is not None

        e3sm_file_exists = self.__caseNcFileExists__(self.e3sm_dir, case_def['e3sm_file'])
        sam_file_exists = self.__caseNcFileExists__(self.sam_folders, case_def['sam_file'])
        wrf_file_exists = self.__caseNcFileExists__(self.wrf_folders, case_def['wrf_file'])
        cam_file_exists = self.__caseNcFileExists__(self.cam_folders, case_def['cam_file'])
        clubb_file_exists = self.__caseNcFileExists__(self.clubb_folders, case_def['clubb_file'])

        e3sm_has_case = e3sm_given and e3sm_file_defined and e3sm_file_exists
        sam_has_case = sam_given and sam_file_defined and sam_file_exists
        wrf_has_case = wrf_given != 0 and wrf_file_defined and wrf_file_exists
        cam_has_case = cam_given and cam_file_defined and cam_file_exists
        clubb_has_case = clubb_given and clubb_file_defined and clubb_file_exists #case_def['name'] in self.clubb_datasets.keys()

        if self.nightly:
            return clubb_has_case and (e3sm_has_case or sam_has_case or cam_has_case or wrf_has_case) \
                   or self.benchmark_only
        else:
            return e3sm_has_case or sam_has_case or cam_has_case or wrf_has_case \
                   or clubb_has_case or self.benchmark_only

    def __caseNcFileExists__(self, list_of_src_folders, rel_filepath):
        """

        :param list_of_src_folders: TODO
        :param rel_filepath: TODO
        :return: TODO
        """
        any_nc_file_found = False
        if rel_filepath is not None and list_of_src_folders is not None:
            for folder in list_of_src_folders:
                if isinstance(rel_filepath, dict):
                    for temp_filename in rel_filepath.values():
                        filename = folder + temp_filename
                        if os.path.exists(filename):
                            any_nc_file_found = True
                else:
                    filename = folder + '/' + rel_filepath
                    if os.path.exists(filename):
                        any_nc_file_found = True
        return any_nc_file_found

    def __copySetupFiles__(self):
        """
        Copies case setup files from the clubb folder(s)
        into the pyplotgen output folders.

        :return: None
        """
        print("Looking for case_setup.txt files")
        for folder in self.clubb_folders:
            setup_file_search_pattern = folder + '/*_setup.txt'
            folder_basename = os.path.basename(folder)

            # Loop through search results of files in folder matching the given pattern
            for file in glob.glob(setup_file_search_pattern):
                # Split filename from entire path
                file_basename = os.path.basename(file)
                # Removing '_setup.txt' from file_basename gives us the case name
                casename = file_basename[:-10]
                # Build destination folder and file names
                copy_dest_folder = self.output_folder + '/' + casename + '/'
                copy_dest_file = copy_dest_folder + file_basename[:-10] + '_' + folder_basename + file_basename[-10:]

                # If case is actually plotted, create output folder and copy setup file to destination folder
                if casename in self.cases_plotted:
                    if not os.path.exists(copy_dest_folder):
                        os.mkdir(copy_dest_folder)
                    shutil.copy(file, copy_dest_file)
                    print("\tFound setup file " + str(file))

    def __benchmarkFilesNeeded__(self):
        """
        Returns true if the user requested to plot model output
        that needs to be downloaded
        
        :return: True if any of self.les or self.hoc or self.cgbest is True, else False 
        """
        return self.les or self.hoc or self.cgbest

    def __downloadModelOutputs__(self):
        """
        Checks for model output in the folder specified under Case_definitions.BENCHMARK_OUTPUT_ROOT,
        e.g. sam benchmark runs, and if it doesn't exist, downloads it.
        Supports: SAM, COAMPS, CLUBB

        :return: None
        """
        # Ensure benchmark output is available
        print("Checking for model benchmark output...")
        # Check if the folder specified in Case_definitions.BENCHMARK_OUTPUT_ROOT exists
        if not os.path.isdir(Case_definitions.BENCHMARK_OUTPUT_ROOT) and \
                not os.path.islink(Case_definitions.BENCHMARK_OUTPUT_ROOT):
            print("\tDownloading the benchmarks to " + Case_definitions.BENCHMARK_OUTPUT_ROOT)
            subprocess.run(['git', 'clone', 'https://carson.math.uwm.edu/les_and_clubb_benchmark_runs.git',
                            Case_definitions.BENCHMARK_OUTPUT_ROOT])
        else:
            print("Benchmark output found in " + Case_definitions.BENCHMARK_OUTPUT_ROOT)


def __trimTrailingSlash__(args):
    """
    Takes in a list filepaths and removes any trailing /'s

    :param arg: list of string file paths
    :return: list of filepaths with trailing /'s removed
    """
    for i in range(len(args)):
        if args[i][-1] == "/":
            args[i] = args[i][:-1]
    return args


def __processArguments__():
    """
    This method takes arguments in from the command line and feeds them into a PyPlotGen object

    :return: A PyPlotGen object containing the parameters as given from the commandline.
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
                        help="Same as -lgd. Plots LES, Golaz Best Ever, and HOC 2005 dependent_data for comparison.",
                        action="store_true")
    parser.add_argument("-z", "--zip", help="Output dependent_data into a compressed zip file.", action="store_true")
    parser.add_argument("--show-alphabetic-id", help="Add an identifying character to the top right of a panel.",
                        action="store_true")
    parser.add_argument("--thin", help="Plot using thin solid lines.", action="store_true")
    parser.add_argument("--no-legends", help="Plot without legend boxes defining the line types.", action="store_true")
    parser.add_argument("--ensemble", help="Plot ensemble tuner runs", action="store_true")  # TODO is this needed?
    parser.add_argument("-b", "--plot-budgets", help="Plot all defined budgets of moments.",
                        action="store_true")
    parser.add_argument("-t", "--time-height-plots",
                        help="Instead of averaged profiles, create contour plots from 2d data."+
                        " Cannot be used with -m.",
                        action="store_true")
    parser.add_argument("-m", "--movies",
                        help="Instead of averaged profiles, plot animations of time steps. Cannot be used with -t."+
                        " Valid arguments: 'gif' or 'mp4'",
                        action="store", choices=['gif','mp4'])
    parser.add_argument("--bu-morr",
                        help="For morrison microphysics: breaks microphysical source terms into component processes",
                        action="store_true")
    parser.add_argument("--benchmark-only",
                        help="Prevents autoplotting of clubb's default output folder when no input folders are "
                             "specified. This results in only plotting the benchmark output, though this "
                             "output doesn't guarantee all text fields or plots are filled.",
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
    parser.add_argument("--nightly", help="Apply special parameters only relevant when running as part of a nightly "
                                          "test. "
                                          "This is currently limited to disabling case output if not all models have"
                                          "data for a given case. E.g. this prevents wrf plots from including cases "
                                          "that "
                                          "only have clubb plots and no wrf plots. Do not plot this with clubb-only "
                                          "plots, just plot clubb normally for clubb nightly tests.",
                        action="store_true")
    args = parser.parse_args()

    if args.zip:
        print("Zip flag detected, but that feature is not yet implemented")
    if args.ensemble:
        print("Ensemble flag detected, but that feature is not yet implemented")
    if args.bu_morr:
        print("Morrison breakdown flag detected, but that feature is not yet implemented")

    # If the last char in folder path is /, remove it
    args.clubb = __trimTrailingSlash__(args.clubb)
    args.sam = __trimTrailingSlash__(args.sam)
    args.cam = __trimTrailingSlash__(args.cam)
    args.e3sm = __trimTrailingSlash__(args.e3sm)
    args.wrf = __trimTrailingSlash__(args.wrf)

    # If no input is specified at all, use the nc files in the default CLUBB output folder
    no_folders_inputted = len(args.e3sm) == 0 and len(args.sam) == 0 and len(args.clubb) == 0 \
                          and len(args.wrf) == 0 and len(args.cam) == 0
    if no_folders_inputted and not args.benchmark_only:
        args.clubb = ["../../output"]

    # Set flags for special dependent_data
    if args.all_best:
        les = True
        cgbest = True
        hoc = True
    else:
        les = args.les
        cgbest = args.plot_golaz_best
        hoc = args.plot_hoc_2005

    if args.time_height_plots and args.movies is not None:
        raise ValueError('Error: Command line parameter -t and -m cannot be used in conjunction.')

    # Call __init__ function of PyPlotGen class defined above and store an instance of that class in pyplotgen
    pyplotgen = PyPlotGen(args.output, clubb_folders=args.clubb, replace=args.replace, les=les, plot_e3sm=args.e3sm,
                          cgbest=cgbest, cam_folders=args.cam, nightly=args.nightly,
                          hoc=hoc, zip=args.zip, thin=args.thin, sam_folders=args.sam,
                          wrf_folders=args.wrf, benchmark_only=args.benchmark_only,
                          no_legends=args.no_legends, ensemble=args.ensemble, budget_moments=args.plot_budgets,
                          bu_morr=args.bu_morr, diff=args.diff, show_alphabetic_id=args.show_alphabetic_id,
                          time_height=args.time_height_plots, animation=args.movies)
    return pyplotgen


if __name__ == "__main__":
    pyplotgen = __processArguments__()
    pyplotgen.run()
