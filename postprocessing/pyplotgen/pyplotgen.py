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
import multiprocessing
import os
import logging
import shutil
import subprocess
import time
from datetime import datetime
from difflib import SequenceMatcher
from multiprocessing import Pool, Array
from multiprocessing import freeze_support

from fpdf import FPDF

from config import Case_definitions, Style_definitions, Case_definitions_red_height
from python_html_gallery import gallery
from src import Panel
from src.CaseGallerySetup import CaseGallerySetup
from src.DataReader import DataReader
from src.interoperability import clean_path
import src.OutputHandler
from src.OutputHandler import logToFile, logToFileAndConsole
from src.OutputHandler import initializeProgress, writeFinalErrorLog, warnUser

class PyPlotGen:
    """
    Main class for the pyplotgen program.
    After processing the command line options in the main function,
    an instance of this class is created in __process_args__, passing all option parameters.
    Finally the run function of this class is the actual driver of the pyplotgen program,
    creating the CaseGallerySetup objects and calling their plot functions to generate the images and html output.
    """

    def __init__(self, output_folder, clubb_folders=None, replace=False, les=False, cgbest=False, hoc=False,
                 benchmark_only=False, nightly=False, zip=False, thin=False, no_legends=False, ensemble=False,
                 e3sm_folders=[""], sam_folders=[""], wrf_folders=[""], cam_folders=[""], priority_vars=False,
                 plot_budgets=False, bu_morr=False, lumped_buoy_budgets=False, background_rcm=False, diff=None,
                 show_alphabetic_id=False, time_height=False, animation=None, samstyle=False, disable_multithreading=False,
                 pdf=False, pdf_filesize_limit=None, plot_subcolumns=False, image_extension=".png",
                 grid_adapt_plot=False, grid_comparison_plot=False, red_height=False):
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
        :param e3sm_folders: Plot E3SM dependent_data for comparison.
            This parameter works exactly like the clubb_folders parameter, except E3SM uses only one nc file.
        :param sam_folders: Plot SAM dependent_data for comparison.
            This parameter works exactly like the clubb_folders parameter, except SAM uses only one nc file.
        :param wrf_folders: Plot WRF dependent_data for comparison.
            This parameter works exactly like the clubb_folders parameter, except WRF uses the following nc files:
            _zm, _zt, and _sfc
        :param cam_folders: Plot CAM dependent_data for comparison.
            This parameter works exactly like the clubb_folders parameter, except CAM uses only one nc file. (CHECK!)
        :param plot_budgets: If True, plot all defined budgets of moments.
        :param bu_morr: For morrison microphysics: If True, break microphysical source terms into component processes.
            Not implemented.
        :param lumped_buoy_budgets: Lump together wpxp_bp and wpxp_pr3 terms in CLUBB's budgets
        :param background_rcm: Show a height-based "contour" plot of time-averaged rcm behind CLUBB profiles.
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
        self.e3sm_folders = e3sm_folders
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
        self.plot_budgets = plot_budgets
        self.plot_subcolumns = plot_subcolumns
        self.bu_morr = bu_morr
        self.lumped_buoy_budgets = lumped_buoy_budgets
        self.background_rcm = background_rcm
        self.diff = diff
        self.cases_plotted = []
        self.clubb_datasets = None
        self.data_reader = DataReader()
        self.diff_files_data_reader = DataReader()
        self.sam_data_reader = DataReader()
        self.show_alphabetic_id = show_alphabetic_id
        self.output_folder = os.path.abspath(self.output_folder)
        self.benchmark_only = benchmark_only
        self.priority_vars = priority_vars
        self.nightly = nightly
        self.time_height = time_height
        self.animation = animation
        self.sam_style_budgets = samstyle
        self.multithreaded = not disable_multithreading
        self.pdf = pdf
        self.pdf_filesize_limit = pdf_filesize_limit
        self.image_extension = image_extension
        self.grid_adapt_plot = grid_adapt_plot
        self.grid_comparison_plot = grid_comparison_plot
        self.red_height = red_height

        # Append date to output folder name
        if os.path.isdir(self.output_folder) and self.replace_images is False:
            current_date_time = datetime.now()
            rounded_down_datetime = current_date_time.replace(microsecond=0)
            datetime_generated_on = str(rounded_down_datetime)
            self.output_folder = self.output_folder + '_generated_on_' + datetime_generated_on

        self.output_folder = clean_path(self.output_folder)

        # If --replace flag was set, delete old output folder
        if self.replace_images:
            subprocess.run(['rm', '-rf', self.output_folder + '/'])
            # TODO: Use for Windows
            # shutil.rmtree(self.output_folder)

        # create output folder to store error log file
        # and configure logging
        try:
            os.mkdir(self.output_folder)
        except FileExistsError:
            pass # do nothing
        self.errorlog = self.output_folder+"/error_temp.log"
        self.finalerrorlog = self.output_folder+"/error.log"
        logging.basicConfig(filename=self.errorlog, filemode='w', level=logging.INFO,
                    format='%(asctime)s.%(msecs)03d %(message)s', datefmt='%y-%m-%d %H:%M:%S')

        # check that all input folders exist
        all_folders = self.clubb_folders + self.e3sm_folders + self.sam_folders \
                    + self.cam_folders+self.wrf_folders
        for folder in all_folders:
            if os.path.isdir(folder)==False:
                raise RuntimeError("The directory " + folder + " does not exist.")

    def run(self):
        """
        Main driver of the pyplotgen program executing the following steps:
        - Download benchmark files if needed
        - Loop through cases listed in Case_definitions.CASES_TO_PLOT
        - Create instances of these cases
        - Call plot function on the instances
        - Create output folder
        - Generate html page containing plots

        :return: None
        """
        logToFileAndConsole('*******************************************')
        logToFileAndConsole("Welcome to PyPlotGen.")
        logToFileAndConsole('*******************************************')
        logToFileAndConsole('                                           ')
        self.diff_datasets = None
        # Load data used for difference plots
        if self.diff is not None:
            self.diff_datasets = self.diff_files_data_reader.loadFolder(self.diff)
        if not self.red_height:
            all_enabled_cases = Case_definitions.CASES_TO_PLOT
        else:
            all_enabled_cases = Case_definitions_red_height.CASES_TO_PLOT

        # Downloads model output (sam, les, clubb) if it doesn't exist
        if self.__benchmarkFilesNeeded__():
            self.__downloadModelOutputs__()
        self.num_cases_plotted = 0
        # Loop through cases listed in Case_definitions.CASES_TO_PLOT
        # for case_def in all_enabled_cases:
        cases_plotted_bools = []

        # initialize counter and progress display
        total_progress_counter = Array('i',[0,0])
        initializeProgress(self.image_extension, self.animation)

        if self.multithreaded:
            freeze_support()  # Required for multithreading
            n_processors = multiprocessing.cpu_count()
            with Pool(processes=n_processors,initializer=tpc_init,initargs=(total_progress_counter, )) as pool:
                cases_plotted_bools = pool.map(self.__plotCase__, all_enabled_cases)
        else:
            for case_def in all_enabled_cases:
                cases_plotted_bools.append(self.__plotCase__(case_def))
        logToFileAndConsole('')
        logToFileAndConsole('-------------------------------------------')

        self.num_cases_plotted = self.__extractNumCasesPlotted__(cases_plotted_bools)

        if self.num_cases_plotted == 0:
            all_cases_casenames = []
            for case in Case_definitions.CASES_TO_PLOT:
                all_cases_casenames.append(case["name"])

            logToFileAndConsole(
                "Error:  No cases were plotted.  Please confirm your input folder path,\n"
                "(i.e. following --sam, --clubb, --e3sm, --wrf, --cam), and confirm that\n"
                "your netcdf filename matches the expected filename given in\n"
                "config/Case_definitions.py (see the 'Make sure filenames match' section\n"
                "of the pyplotgen README for more info).  Also confirm that the CASES_TO_PLOT\n"
                "variable in config/Case_Defintions.py lists all of the cases you expected\n"
                "to plot.\n"
                "CASES_TO_PLOT = " + str(all_cases_casenames) +
                "\nPlease run ./pyplotgen.py -h for more information on parameters.")
        logToFile("Generating webpage for viewing plots ")

        if not os.path.exists(self.output_folder):
            os.mkdir(self.output_folder)

        self.__copySetupFiles__()
        # Generate html pages
        # Multithreading changes the order the cases are plotted on the webpage, so it has been disabled.
        # The capability is being left here as a demo.
        if self.animation is not None:
            movie_extension = "." + self.animation
            gallery.main(self.output_folder, multithreaded=False, file_extension=movie_extension)
        else:
            gallery.main(self.output_folder, multithreaded=False, file_extension=self.image_extension)
        logToFileAndConsole('-------------------------------------------')
        logToFileAndConsole("Output can be viewed at file://" + self.output_folder + "/index.html with a web browser")

    def __printToPDF__(self):
        """
        If --pdf was specified, this prints a pdf. Otherwise, this does nothing.
        If --pdf and --pdf-filesize-limit were specified, this will loop over the run() method, reducing image size on
        each loop, until a pdf is output less than or equal to the specified target pdf filesize.

        :return: None
        """

        lowest_output_folder_level = str.split(self.output_folder, '/')[-1]
        pdf_output_path_plus_filename = self.output_folder + '/' + lowest_output_folder_level + '.pdf'
        case_descriptions = {}
        case_times = {}
        for case in Case_definitions.ALL_CASES:
            case_descriptions[case['name']] = case['description']
            case_times[case['name']] = [case['start_time'], case['end_time']]
        if self.pdf and self.pdf_filesize_limit is None:
            logToFileAndConsole('-------------------------------------------')
            logToFileAndConsole('Generating PDF file ' + pdf_output_path_plus_filename)
            # config = pdfkit.configuration(wkhtmltopdf='/usr/bin/')
            self.__writePdfToDisk__(pdf_output_path_plus_filename, case_descriptions, case_times)
            # pdfkit.from_file(html_input_filename, pdf_output_path_plus_filename)
            logToFileAndConsole("PDF Output can be viewed at file://" + pdf_output_path_plus_filename + " with a web browser/ pdf viewer")
            logToFileAndConsole('-------------------------------------------')
        if self.pdf_filesize_limit is not None:
            output_dpi = Style_definitions.IMG_OUTPUT_DPI
            pdf_too_large = True
            filesize_impossible = False
            logToFileAndConsole('-------------------------------------------')
            logToFileAndConsole('Generating PDF file ' + pdf_output_path_plus_filename)
            logToFileAndConsole('Searching for minimum viable DPI for output images to print within '
                  + str(self.pdf_filesize_limit) + 'MB')
            attempted_prints = 0
            while pdf_too_large and not filesize_impossible:
                attempted_prints += 1
                bytes_to_mb = 1 / 1000000
                self.__writePdfToDisk__(pdf_output_path_plus_filename, case_descriptions, case_times)
                pdf_filesize = os.path.getsize(pdf_output_path_plus_filename) * bytes_to_mb
                logToFileAndConsole("PDF generated using a DPI of " + str(Style_definitions.IMG_OUTPUT_DPI) +
                                    " with a filesize of " + str(pdf_filesize) + "MB.")
                Style_definitions.IMG_OUTPUT_DPI = output_dpi

                if pdf_filesize < self.pdf_filesize_limit:
                    pdf_too_large = False
                    logToFileAndConsole("PDF output can be found at: file://" + pdf_output_path_plus_filename)
                    logToFileAndConsole("Printing PDF to the target filesize took " + str(attempted_prints) +
                                        " attempts to find the right DPI.")
                else:
                    # output downscaled images to a new folder so that the original quality ones can still be referenced
                    # via the web page
                    if "_downscaled_for_pdf" not in self.output_folder:
                        self.output_folder = self.output_folder + "_downscaled_for_pdf"
                        pdf_output_path_plus_filename = self.output_folder + '/pyplotgen_output.pdf'

                    output_dpi = self.__getDecreasedDpiValue__(output_dpi, pdf_filesize)
                    logToFileAndConsole("Attempted to print but the file was too large (" + str(pdf_filesize) +
                          "MB insead of <" + str(self.pdf_filesize_limit) + "MB). Reducing DPI and trying again.")
                    logToFileAndConsole("Attempting to print pdf with dpi of " + str(output_dpi))

                    # delete 'downscaled' output folder and run again
                    subprocess.run(['rm', '-rf', self.output_folder + '/'])
                    self.run()
                if output_dpi <= 1:
                    logToFileAndConsole("There is no possible dpi that fits within " + str(self.pdf_filesize_limit) + "MB.")
                    logToFileAndConsole("The most recent PDF output attempt can be found at: file://" + pdf_output_path_plus_filename)
                    filesize_impossible = True

    def __writePdfToDisk__(self, pdf_output_filename, case_descriptions, case_times):
        """
        This is a helper function that actually writes the PDF to the disk. This uses the fpdf package to generate the
        pdf.

        :param pdf_output_filename: Name of the file to be created. While it can be a relative name, it's recommended
            to use a full file path.
        :param case_descriptions: A dict of name -> description maps. E.g. {'bomex': "I am the bomex case. Fear me!"}
        :return: None
        """
        pdf = FPDF()
        for foldername in sorted(os.listdir(self.output_folder)):
            if os.path.isdir(foldername):
                pdf.add_page()
                pdf.set_font('Arial', 'B', 18)
                pdf.cell(0, 10, foldername + " minutes " + str(case_times[foldername][0]) + "-" + str(case_times[foldername][1]))
                pdf.ln()
                pdf.set_font('Arial', '', 12)
                pdf.multi_cell(0, 8, case_descriptions[foldername])
                current_date_time = datetime.now()
                rounded_down_datetime = str(current_date_time.replace(microsecond=0))
                pdf.set_font('Arial', '', 10)
                pdf.multi_cell(0, 6, "Generated on: " + rounded_down_datetime)
                loop_counter = 0
                num_imgs_per_row = 3
                for filename in sorted(os.listdir(self.output_folder + "/" + foldername)):
                    filename = self.output_folder + '/' + foldername + '/' + filename

                    if "html" not in filename and "txt" not in filename and os.path.isfile(filename):
                        x_coord = 20 + loop_counter * 60
                        pdf.set_x(x_coord)
                        pdf.image(filename, w=50, h=30)
                        pdf.set_y(pdf.get_y() - 30)
                        if loop_counter == num_imgs_per_row - 1:
                            pdf.ln()
                            loop_counter = 0
                            pdf.set_y(pdf.get_y() + 25)
                        else:
                            loop_counter += 1

        pdf.output(pdf_output_filename, 'F')

    def __getDecreasedDpiValue__(self, previous_dpi, previous_filesize):
        """
        When printing files to PDF, dpi lowering the dpi is useful for fitting the pdf's filesize under a certain value.
        This method returns a suggested dpi value to print the pdf to that's more likely to fit within a given filesize.
        Note, it is not guranteed that this function will return a dpi that DOES fit within a certain file size. This
        function is intended to be used multiple times throughout a loop as it gets closer to an optimal dpi.
        This function takes an incremental methodology to approaching an optimal dpi.

        Since the relationship between dpi and pdf filesize is mostly linear, the algorithm used to calculate the next
        recommended dpi is based on subtracting a percentage of the dpi relative to how far off the filesize was.

        An optional step_multiplier can be modified in the code to increase or decrease the strength of the dpi change.

        This method was inspired by the optimization functions used in artificial intelligence.

        :param previous_dpi: A dpi value that was larger than desired. The dpi value returend will be less than this one
        :param previous_filesize: The filesize of the pdf generated by the previous_dpi
        :return: A recommended dpi value that is closer to an optimal dpi value for the target filesize
        """
        # Values >1 make run faster, but may not find the best dpi, values <1 are slower but are more likely to find the
        # best dpi. It's recommended to use values from 0->1. Recommended value: 0.25.
        # General rule of thumb:
        #   1.0 takes ~2-3 iterations, 0.75 takes ~3 iterations, 0.5 takes ~4 iterations, 0.25 takes ~5 iterations
        # Note: A step size of 0 or close to 0 will result in dpi reductions of 1 per iteration. This is the slowest
        #       option but it's guaranteed to find the best dpi.
        step_multiplier = 0.25
        filesize_ratio = (previous_filesize / self.pdf_filesize_limit)
        dpi_reduction_amount = max(previous_dpi * step_multiplier / filesize_ratio, 1)
        new_dpi = round(previous_dpi - dpi_reduction_amount)
        return new_dpi

    def __extractNumCasesPlotted__(self, plotCaseDataArray):
        """
        The __plotCase__() function is often called via a pooling map.
        As such, the data returned from the pool map is an array of True/False values
        where each True represents a case being plotted and each False represents a case
        not being plotted. This function counts the number of True's and retuns that.

        :param plotCaseDataArray: Array of True/False values
        :return: number of True's in the given list.
        """
        num_true = 0
        for element in plotCaseDataArray:
            if element is True:
                num_true += 1
        return num_true

    def __plotCase__(self, case_def):
        """
        Plots the given case. Pulled out into it's own function to help move towards multithreading
        :param case_def:
        :return: True if case was plotted, False if case was not plotted
        """
        self.case_diff_datasets = None
        casename = case_def['name']
        case_plotted = False
        if self.__dataForCaseExists__(case_def):
            logToFile('-------------------------------------------')
            logToFile("Processing: {}".format(case_def['name'].upper()))
            if self.diff is not None:
                self.case_diff_datasets = self.diff_datasets[casename]
            case_gallery_setup = CaseGallerySetup(case_def, clubb_folders=self.clubb_folders, plot_les=self.les,
                                                  plot_budgets=self.plot_budgets, sam_folders=self.sam_folders,
                                                  wrf_folders=self.wrf_folders, diff_datasets=self.case_diff_datasets,
                                                  plot_r408=self.cgbest, plot_hoc=self.hoc, e3sm_folders=self.e3sm_folders,
                                                  cam_folders=self.cam_folders, time_height=self.time_height,
                                                  animation=self.animation, samstyle=self.sam_style_budgets,
                                                  plot_subcolumns=self.plot_subcolumns, lumped_buoy_budgets=self.lumped_buoy_budgets,
                                                  background_rcm=self.background_rcm, image_extension=self.image_extension,
                                                  total_panels_to_plot=0, priority_vars=self.priority_vars)
            # Call plot function of case instance
            case_gallery_setup.plot(self.output_folder, replace_images=self.replace_images, no_legends=self.no_legends,
                                    thin_lines=self.thin, show_alphabetic_id=self.show_alphabetic_id,
                                    total_progress_counter=total_progress_counter, 
                                    generate_grid_adapt_plot = self.grid_adapt_plot, 
                                    grid_comparison_plot = self.grid_comparison_plot,
                                    read_file_paths=self.clubb_folders)
            self.cases_plotted.append(case_def)
            case_plotted = True

        return case_plotted

    def __dataForCaseExists__(self, case_def):
        """
        Returns true if there's an nc file for a given case name and any model that should be plotted.

        :param case_def: The case definition object
        :return: True if data exists, False if not
        """
        e3sm_given = len(self.e3sm_folders) != 0
        sam_given = len(self.sam_folders) != 0
        wrf_given = len(self.wrf_folders) != 0
        cam_given = len(self.cam_folders) != 0
        clubb_given = len(self.clubb_folders) != 0

        e3sm_file_defined = case_def['e3sm_file'] is not None
        sam_file_defined = case_def['sam_file'] is not None
        wrf_file_defined = case_def['wrf_file'] is not None
        cam_file_defined = case_def['cam_file'] is not None
        clubb_file_defined = case_def['clubb_file'] is not None

        e3sm_file_exists = self.__caseNcFileExists__(self.e3sm_folders, case_def['e3sm_file'])
        sam_file_exists = self.__caseNcFileExists__(self.sam_folders, case_def['sam_file'])
        wrf_file_exists = self.__caseNcFileExists__(self.wrf_folders, case_def['wrf_file'])
        cam_file_exists = self.__caseNcFileExists__(self.cam_folders, case_def['cam_file'])
        clubb_file_exists = self.__caseNcFileExists__(self.clubb_folders, case_def['clubb_file'])

        e3sm_has_case = e3sm_given and e3sm_file_defined and e3sm_file_exists
        sam_has_case = sam_given and sam_file_defined and sam_file_exists
        wrf_has_case = wrf_given != 0 and wrf_file_defined and wrf_file_exists
        cam_has_case = cam_given and cam_file_defined and cam_file_exists
        clubb_has_case = clubb_given and clubb_file_defined and clubb_file_exists  # case_def['name'] in self.clubb_datasets.keys()

        if self.nightly:
            return (clubb_has_case or not clubb_file_defined) \
                   and (e3sm_has_case or sam_has_case or cam_has_case or wrf_has_case) \
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
        logToFile("Looking for case_setup.txt files")
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
                if os.path.exists(copy_dest_folder):
                    shutil.copy(file, copy_dest_file)
                    logToFile("\tFound setup file " + str(file))

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
        logToFileAndConsole("Checking for model benchmark output...")
        # Check if the folder specified in Case_definitions.BENCHMARK_OUTPUT_ROOT exists
        if not os.path.isdir(Case_definitions.BENCHMARK_OUTPUT_ROOT) and \
                not os.path.islink(Case_definitions.BENCHMARK_OUTPUT_ROOT):
            logToFileAndConsole("\tDownloading the benchmarks to {}.".format(Case_definitions.BENCHMARK_OUTPUT_ROOT))
            subprocess.run(['git', 'clone', 'https://carson.math.uwm.edu/les_and_clubb_benchmark_runs.git',
                            Case_definitions.BENCHMARK_OUTPUT_ROOT])
        else:
            logToFileAndConsole("Benchmark output found in {}.".format(Case_definitions.BENCHMARK_OUTPUT_ROOT))
        logToFileAndConsole('-------------------------------------------')


def __convertCasenamesToCaseInstances__(casenames_list):
    """
    Takes in a list of case names and searches through all the cases in Case_definitions.py for matches.
    Will then return a tuple containing the case instances that were found and a list of any names that were not found.

    :param casenames_list: List of case names. Names must match the 'name' parameter of the cases definition in
        Case_definitions.py
    :return: Tuple of (converted_cases, cases_not_found) where converted_cases is a list of case instances and
        cases_not_found is a list of casenames that could not be found.
    """
    converted_cases = []
    cases_not_found = []
    for casename in casenames_list:
        case_found = False
        for case_def in Case_definitions.ALL_CASES:
            if case_def['name'] == casename:
                case_found = True
                converted_cases.append(case_def)
                break
        if not case_found:
            cases_not_found.append(casename)

    return converted_cases, cases_not_found


def __getSimilarCasenames__(invalid_casenames):
    """
    Given a list of casenames, this method returns a list of existing case names that are closest to the list given.

    E.g. if the list given contains ['bomax', 'dycoms_rf01'] this function returns ['bomex, 'dycoms2_rf01']

    :param invalid_casenames: List of casenames that could not be found/likely have typos
    :return:
    """
    valid_casenames = []
    suggested_casenames = []
    for case_def in Case_definitions.ALL_CASES:
        valid_casenames.append(case_def['name'])
    for invalid_casename in invalid_casenames:
        best_probability = 0
        best_name_match = "unknown"
        for valid_casename in valid_casenames:
            test_probability = SequenceMatcher(None, invalid_casename, valid_casename).ratio()
            if test_probability > best_probability:
                best_probability = test_probability
                best_name_match = valid_casename
        suggested_casenames.append(best_name_match)

    return suggested_casenames


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
    parser.add_argument("--grid-adapt-plot", help="Plot grid adaptation plot.", action="store_true")
    parser.add_argument("--red-height", help="Plot with reduced height.", action="store_true")
    parser.add_argument("--grid-comparison-plot", help="Plot grid comparison plot.", action="store_true")
    parser.add_argument("-b", "--plot-budgets", help="Plot all defined budgets of moments.",
                        action="store_true")
    parser.add_argument("--lumped-buoy-budgets", help="Lump together wpxp_bp and wpxp_pr3 terms in CLUBB's budgets.",
                        action="store_true")
    parser.add_argument("--background-rcm", help="Show a height-based 'contour' plot of time-averaged rcm behind CLUBB profiles.",
                        action="store_true")
    parser.add_argument("--plot-subcolumns", help="Plot all defined subcolumns.",
                        action="store_true")
    parser.add_argument("-t", "--time-height-plots",
                        help="Instead of averaged profiles, create contour plots from 2d data." +
                             " Cannot be used with -m.",
                        action="store_true")
    parser.add_argument("-m", "--movies",
                        help="Instead of averaged profiles, plot animations of time steps. Cannot be used with -t, " +
                             "--pdf, --eps, or --svg. FRAMES_PER_SECOND can be adjusted in " +
                             "config/Style_definitions.py.",
                        action="store", nargs='?', const='mp4', choices=['mp4','avi'])
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
    parser.add_argument("--disable-multithreading", help="This forces pyplotgen to run on a single thread. This isn't "
                                                         "recommended, but if you're having issues with multithreading "
                                                         "or want text output to appear sequentially and don't care"
                                                         "about performance, use this option.",
                        action="store_true")
    parser.add_argument("--high-quality", "--hq",
                        help="Outputs higher resolution images. The dpi used for hi resolution images"
                             " can be customized in Style_definitions.py",
                        action="store_true")
    parser.add_argument("--svg",
                        help="Outputs images to a lossless vector-graphics format .svg  instead of a rasterized image "
                             "like png.This argument is not compatible with --eps or --pdf. "
                             "Note that some browsers may have difficulty displaying this format.",
                        action="store_true")
    parser.add_argument("--eps",
                        help="Outputs images to Encapsulated PostScript eps instead of a rasterized image "
                             "like png. This argument is not compatible with --svg or --pdf. "
                             "Note that some browsers may have difficulty displaying this format.",
                        action="store_true")
    parser.add_argument("--pdf",
                        help="In addition to pyplotgen's regular output, this also generates a pdf file from the "
                             "plots.html output. Note that this argument depends on wkhtmltopdf. To install wkhtmltopdf"
                             " please visit this page: "
                             "https://github.com/JazzCore/python-pdfkit/wiki/Installing-wkhtmltopdf",
                        action="store_true")
    parser.add_argument("--pdf-filesize-limit", help="Adjust pdf filesize so that it is no larger than the given size "
                                                     "in MB. Note that this argument only works if --pdf is also "
                                                     "specified",
                        action="store", type=int)
    parser.add_argument("--cases", help="A set of case name(s) to be ran. Cases not listed here will not be ran. The "
                                        "casename specified must match the 'name' parameter of the case's definition "
                                        "Case_definitions.py. E.g. --cases bomex arm wangara",
                        action="store",
                        default=[], nargs='+')
    parser.add_argument("--priority-variables", help="Plot only variables with the 'priority' key.",
                        action="store_true")
    parser.add_argument("--sam-style-budgets", help="Lump together certain CLUBB budget terms so that the relevant "
                                                    "CLUBB budgets look comparable to SAM's budgets.",
                        action="store_true")
    args = parser.parse_args()

    if args.zip:
        logToFileAndConsole("Zip flag detected, but that feature is not yet implemented")
    if args.bu_morr:
        logToFileAndConsole("Morrison breakdown flag detected, but that feature is not yet implemented")

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

    image_extension = ".png"
    if args.movies is not None and (args.svg or args.eps):
       raise RuntimeError("The --movies option currently only works with .png images.  Please remove --eps or --svg "
                          "in order to generate animated plots.")
    if args.svg:
        image_extension = ".svg"
    if args.eps:
        image_extension = ".eps"

    if args.eps and args.svg:
        raise RuntimeError("The --svg and --eps options are not compatible with one another. Please select either --eps "
                           "or --svg but not both.")

    if (args.eps or args.svg) and args.pdf:
        raise RuntimeError("SVG and EPS are not supported alongside the pdf parameter. This is due to a limitation of "
                           "the FPDF engine used. Please remove either the --svg or --eps option (whichever was used) or "
                           "remove the --pdf option. Note that you can create --svg/--eps output if desired, and then "
                           "rerun pyplotgen without that option to produce pdf output.")

    if args.high_quality:
        Style_definitions.IMG_OUTPUT_DPI = Style_definitions.HQ_DPI

    if args.time_height_plots and args.movies is not None:
        raise ValueError('Error: Command line parameter -t and -m cannot be used in conjunction.')

    if args.pdf and args.movies is not None:
        raise ValueError('Error: Command line parameters --pdf and --movies cannot be used in conjunction.')

    if args.time_height_plots and args.plot_budgets:
        raise ValueError('Error: Command line parameters --time-height-plots and -b (--plot-budgets) cannot '
                         'be used in conjunction.')

    if len(args.cases) > 0:
        cases_to_plot, cases_not_found = __convertCasenamesToCaseInstances__(args.cases)
        Case_definitions.CASES_TO_PLOT = cases_to_plot
        suggested_casenames = __getSimilarCasenames__(cases_not_found)
        if len(cases_not_found) > 0:
            raise NameError("Error. The following cases could not be found in Case_definitions.py: " +
                            str(cases_not_found) + ".\n"
                                                   "Perhaps you meant " + str(suggested_casenames) + "?")

    # Call __init__ function of PyPlotGen class defined above and store an instance of that class in pyplotgen
    pyplotgen = PyPlotGen(args.output, clubb_folders=args.clubb, replace=args.replace, les=les, e3sm_folders=args.e3sm,
                          cgbest=cgbest, cam_folders=args.cam, nightly=args.nightly,
                          hoc=hoc, zip=args.zip, thin=args.thin, sam_folders=args.sam,
                          wrf_folders=args.wrf, benchmark_only=args.benchmark_only, priority_vars=args.priority_variables,
                          no_legends=args.no_legends, plot_budgets=args.plot_budgets, bu_morr=args.bu_morr,
                          lumped_buoy_budgets=args.lumped_buoy_budgets, background_rcm=args.background_rcm, diff=args.diff,
                          show_alphabetic_id=args.show_alphabetic_id, time_height=args.time_height_plots, animation=args.movies,
                          samstyle=args.sam_style_budgets, disable_multithreading=args.disable_multithreading, pdf=args.pdf,
                          pdf_filesize_limit=args.pdf_filesize_limit, plot_subcolumns=args.plot_subcolumns,
                          image_extension=image_extension, grid_adapt_plot=args.grid_adapt_plot,
                          grid_comparison_plot=args.grid_comparison_plot, red_height=args.red_height)
    return pyplotgen


# Added to track progress with multithreading
def tpc_init(x):
    global total_progress_counter
    total_progress_counter = x


if __name__ == "__main__":
    pyplotgen = __processArguments__()
    start_time = time.time()
    pyplotgen.run()
    pyplotgen.__printToPDF__()
    total_runtime = round(time.time() - start_time)
    logToFileAndConsole("Pyplotgen ran in {} seconds.".format(total_runtime))
    writeFinalErrorLog(pyplotgen.errorlog,pyplotgen.finalerrorlog)
    print("See error.log in the output folder for detailed info including warnings.\n")
    warnUser(pyplotgen.finalerrorlog)