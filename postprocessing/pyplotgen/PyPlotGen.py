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
from _warnings import warn

from pyplotgen.Case_arm import Case_arm
from pyplotgen.Case_arm_97 import Case_arm_97
from pyplotgen.Case_astex_a209 import Case_astex_a209
from pyplotgen.Case_atex import Case_atex
from pyplotgen.Case_bomex import Case_bomex
from pyplotgen.Case_cgils_s11 import Case_cgils_s11
from pyplotgen.Case_cgils_s12 import Case_cgils_s12
from pyplotgen.Case_cgils_s6 import Case_cgils_s6
from pyplotgen.Case_clex9_nov02 import Case_clex9_nov02
from pyplotgen.Case_clex9_oct14 import Case_clex9_oct14
from pyplotgen.Case_dycoms2_rf01_fixed_sst import Case_dycoms2_rf01_fixed_sst
from pyplotgen.Case_dycoms2_rf02_do import Case_dycoms2_rf02_do
from pyplotgen.Case_dycoms2_rf02_ds import Case_dycoms2_rf02_ds
from pyplotgen.Case_fire import Case_fire
from pyplotgen.Case_gabls2 import Case_gabls2
from pyplotgen.Case_gabls3 import Case_gabls3
from pyplotgen.Case_gabls3_night import Case_gabls3_night
from pyplotgen.Case_jun25_altocu import Case_jun25_altocu
from pyplotgen.Case_lba import Case_lba
from pyplotgen.Case_mc3e import Case_mc3e
from pyplotgen.Case_mpace_a import Case_mpace_a
from pyplotgen.Case_mpace_b import Case_mpace_b
from pyplotgen.Case_mpace_b_silhs import Case_mpace_b_silhs
from pyplotgen.Case_nov11_altocu import Case_nov11_altocu
from pyplotgen.Case_rico import Case_rico
from pyplotgen.Case_twp_ice import Case_twp_ice
from pyplotgen.Case_wangara import Case_wangara
from pyplotgen.DataReader import DataReader
from pyplotgen.Case_dycoms2_rf01 import Case_dycoms2_rf01


class PyPlotGen:
    def __init__(self, input_folder, output_folder, replace=False, les=False, cgbest=False, hoc=False, plotrefs=False,
                 nightly=False, hq_imgs=False, eps=False, zip=False, thin=False, no_legends=False, ensemble=False,
                 budget_moments=False, bu_morr=False, diff=False):
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
        self.sam_data_reader = DataReader()

    def run(self):
        '''
        Runs PyPlotGen
        :return: n/a
        '''
        self.nc_datasets = self.data_reader.loadFolder(self.input_folder)
        # self.sam_datasets = self.sam_data_reader.loadFolder("/home/strike/sam_benchmark_runs")
        # cases_plotted = {}
        for case_key in self.nc_datasets.keys():
            ncdf_files = self.nc_datasets[case_key]
            dataset_filenames = [dataset.filepath() for dataset in self.nc_datasets[case_key].values()]
            print("##################################################")
            print("Plotting case ", case_key)
            print("##################################################")
            print("Case will be built with the following file(s):")
            for file in dataset_filenames:
                print("\t", file)
            print()
            if case_key == 'astex_a209':
                # print("\tSkipping case")
                astex_a209 = Case_astex_a209(ncdf_files, plot_sam=self.les)
                astex_a209.plot(self.output_folder)
            elif case_key == 'arm_97':
                # print("\tSkipping case")
                arm_97 = Case_arm_97(ncdf_files, plot_sam=self.les)
                arm_97.plot(self.output_folder)
            elif case_key == 'arm':
                print("\tSkipping case")
                cgils_arm = Case_arm(ncdf_files, plot_sam=self.les)
                cgils_arm.plot(self.output_folder)
            elif case_key == 'atex':
                print("\tSkipping case")
                atex = Case_atex(ncdf_files, plot_sam=self.les)
                atex.plot(self.output_folder)
            elif case_key == 'bomex':
                print("\tSkipping case")
                bomex = Case_bomex(ncdf_files, plot_sam=self.les)
                bomex.plot(self.output_folder)
            elif case_key == 'cgils_s6':
                # print("\tSkipping case")
                cgils_s6 = Case_cgils_s6(ncdf_files, plot_sam=self.les)
                cgils_s6.plot(self.output_folder)
            elif case_key == 'cgils_s11':
                # print("\tSkipping case")
                cgils_s11 = Case_cgils_s11(ncdf_files, plot_sam=self.les)
                cgils_s11.plot(self.output_folder)
            elif case_key == 'cgils_s12':
                # print("\tSkipping case")
                cgils_s12 = Case_cgils_s12(ncdf_files, plot_sam=self.les)
                cgils_s12.plot(self.output_folder)
            elif case_key == 'clex9_oct14':
                # print("\tSkipping case")
                clex9_oct14 = Case_clex9_oct14(ncdf_files, plot_sam=self.les)
                clex9_oct14.plot(self.output_folder)
            elif case_key == 'clex9_nov02':
                # print("\tSkipping case")
                clex9_nov02 = Case_clex9_nov02(ncdf_files, plot_sam=self.les)
                clex9_nov02.plot(self.output_folder)
            elif case_key == 'lba':
                # print("\tSkipping case")
                lba = Case_lba(ncdf_files, plot_sam=self.les)
                lba.plot(self.output_folder)
            elif case_key == 'dycoms2_rf01':
                # print("\tSkipping case")
                dycoms_rf01 = Case_dycoms2_rf01(ncdf_files, plot_sam=self.les)
                dycoms_rf01.plot(self.output_folder)
            elif case_key == 'dycoms2_rf01_fixed_sst':
                # print("\tSkipping case")
                dycoms_rf01_fixed_sst = Case_dycoms2_rf01_fixed_sst(ncdf_files, plot_sam=self.les)
                dycoms_rf01_fixed_sst.plot(self.output_folder)
            elif case_key == 'dycoms2_rf02_ds':
                # print("\tSkipping case")
                dycoms_rf02_ds = Case_dycoms2_rf02_ds(ncdf_files, plot_sam=self.les)
                dycoms_rf02_ds.plot(self.output_folder)
            elif case_key == 'dycoms2_rf02_do':
                # print("\tSkipping case")
                dycoms_rf02_do = Case_dycoms2_rf02_do(ncdf_files, plot_sam=self.les)
                dycoms_rf02_do.plot(self.output_folder)
            elif case_key == 'fire':
                # print("\tSkipping case")
                fire = Case_fire(ncdf_files, plot_sam=self.les)
                fire.plot(self.output_folder)
            elif case_key == 'gabls2':
                # print("\tSkipping case")
                gabls2 = Case_gabls2(ncdf_files, plot_sam=self.les)
                gabls2.plot(self.output_folder)
            elif case_key == 'gabls3':
                # print("\tSkipping case")
                gabls3 = Case_gabls3(ncdf_files, plot_sam=self.les)
                gabls3.plot(self.output_folder)
            elif case_key == 'gabls3_night':
                # print("\tSkipping case")
                gabls3_night = Case_gabls3_night(ncdf_files, plot_sam=self.les)
                gabls3_night.plot(self.output_folder)
            elif case_key == 'jun25_altocu':
                # print("\tSkipping case")
                jun25_altocu = Case_jun25_altocu(ncdf_files, plot_sam=self.les)
                jun25_altocu.plot(self.output_folder)
            elif case_key == 'mc3e':
                # print("\tSkipping case")
                mc3e = Case_mc3e(ncdf_files, plot_sam=self.les)
                mc3e.plot(self.output_folder)
            elif case_key == 'mpace_a':
                # print("\tSkipping case")
                mpace_a = Case_mpace_a(ncdf_files, plot_sam=self.les)
                mpace_a.plot(self.output_folder)
            elif case_key == 'mpace_b':
                # print("\tSkipping case")
                mpace_b = Case_mpace_b(ncdf_files, plot_sam=self.les)
                mpace_b.plot(self.output_folder)
            elif case_key == 'mpace_b_silhs':
                # print("\tSkipping case")
                mpace_b_silhs = Case_mpace_b_silhs(ncdf_files, plot_sam=self.les)
                mpace_b_silhs.plot(self.output_folder)
            elif case_key == 'nov11_altocu':
                # print("\tSkipping case")
                nov11_altocu = Case_nov11_altocu(ncdf_files, plot_sam=self.les)
                nov11_altocu.plot(self.output_folder)
            elif case_key == 'rico':
                # print("\tSkipping case")
                rico = Case_rico(ncdf_files, plot_sam=self.les)
                rico.plot(self.output_folder)
            elif case_key == 'twp_ice':
                # print("\tSkipping case")
                twp_ice = Case_twp_ice(ncdf_files, plot_sam=self.les)
                twp_ice.plot(self.output_folder)
            elif case_key == 'wangara':
                # print("\tSkipping case")
                wangara = Case_wangara(ncdf_files, plot_sam=self.les)
                wangara.plot(self.output_folder)
            else:
                warn("Case " + case_key + " has not been implemented or is commented out. Skipping. ")

        print("##################################################")
        print("Generating webpage for viewing plots ")
        print("##################################################")
        subprocess.run(['sigal', 'build', '-f', self.output_folder + '/'])  # Use sigal to build html in '_build/'


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
    parser.add_argument("-n", "--nightly", help="Run in nightly mode.", action="store_true")  # TODO update this description with what it means
    parser.add_argument("-q", "--high-quality-imgs", help="Output high quality images (does not auto scale)",
                        action="store_true")  # TODO update if this does auto scale
    parser.add_argument("-e", "--save-eps", help="Does not delete EPS images after conversion.", action="store_true")  # TODO what does this mean?
    parser.add_argument("-z", "--zip", help="Output data into a compressed zip file.", action="store_true")
    parser.add_argument("--thin", help="Plot using thin solid lines.", action="store_true")
    parser.add_argument("--no-legends", help="Plot without legend boxes defining the line types.", action="store_true")
    parser.add_argument("--ensemble", help="Plot ensemble tuner runs", action="store_true")  # TODO is this needed?
    parser.add_argument("--budget-moments", help="Plot all defined budgets of moments", action="store_true")
    parser.add_argument("--bu-morr", help="For morrison microphysics: breaks microphysical source terms into component processes",
                        action="store_true")
    parser.add_argument("--diff", help="Plot the difference between two input folders", action="store_true")
    parser.add_argument("input", help="Input folder containing netcdf output data.", action="store")
    parser.add_argument("output", help="Name of folder to create and store plots into.", action="store")
    args = parser.parse_args()

    if args.replace:
        print("Replace flag detected, but that feature is not yet implemented")
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
                          hoc=args.plot_hoc_2005, plotrefs=args.all_best, nightly=args.nightly,
                          hq_imgs=args.high_quality_imgs, eps=args.save_eps, zip=args.zip, thin=args.thin,
                          no_legends=args.no_legends, ensemble=args.ensemble, budget_moments=args.budget_moments,
                          bu_morr=args.bu_morr, diff=args.diff)
    return pyplotgen


if __name__ == "__main__":
    pyplotgen = process_args()
    pyplotgen.run()
