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
from pyplotgen.Case_dycoms2_rf02_nd import Case_dycoms2_rf02_nd
from pyplotgen.Case_dycoms2_rf02_so import Case_dycoms2_rf02_so
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
    '''

    '''

    def __init__(self, input_folder, output_folder, replace=False, les=False, cgbest=False, hoc=False, plotrefs=False,
                zip=False, thin=False, no_legends=False, ensemble=False,
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
        self.replace_images = replace
        self.les = les
        self.cgbest = cgbest
        self.hoc = hoc
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
        self.sam_data_reader = DataReader()

    def run(self):
        '''
        Runs PyPlotGen
        :return: n/a
        '''
        self.nc_datasets = self.data_reader.loadFolder(self.input_folder)
        cases = [Case_astex_a209, Case_arm, Case_arm_97, Case_arm_97, Case_atex, Case_bomex, Case_cgils_s6, Case_cgils_s11, Case_cgils_s12,
                 Case_clex9_oct14, Case_clex9_nov02, Case_dycoms2_rf01, Case_dycoms2_rf01_fixed_sst, Case_dycoms2_rf02_do, Case_dycoms2_rf02_ds, Case_dycoms2_rf02_nd, Case_dycoms2_rf02_so,
                 Case_fire, Case_gabls2, Case_gabls3, Case_gabls3_night, Case_jun25_altocu, Case_lba, Case_mc3e, Case_mpace_a, Case_mpace_b,
                 Case_mpace_b_silhs, Case_nov11_altocu, Case_rico, Case_twp_ice, Case_wangara]
        # cases = [Case_gabls2]

        for case in cases:
            # try:
            self.run_case(case, self.nc_datasets[case.name])
            # except (KeyError):
                # raise FileNotFoundError("The dataset for " + case.name + " was not found in " + self.output_folder +
                #                         ". Please make sure the dataset exists and uses the nameing pattern " + case.name + "_EXT.nc")

        print('###########################################')
        print("\nGenerating webpage for viewing plots ")
        subprocess.run(['sigal', 'build', '-f', self.output_folder + '/'])  # Use sigal to build html in '_build/'

    def run_case(self, case, ncdf_files):
        '''
        Run a case
        :param case: The case class object ot be ran. E.g. pass in  `Case_astex_a209 class` as in the class name and NOT an instance
        :param ncdf_files: Dictionary of netcdf files to load data from
        :return: none
        '''
        print('###########################################')
        print("plotting ", case.name)
        temp_case = case(ncdf_files, plot_sam=self.les, plot_budgets = self.plot_budgets)
        temp_case.plot(self.output_folder, replace_images=self.replace_images, no_legends = self.no_legends)
        print("done plotting ", case.name)

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
    if args.zip:
        print("Zip flag detected, but that feature is not yet implemented")
    if args.thin:
        print("Thin plotAll lines flag detected, but that feature is not yet implemented")
    if args.no_legends:
        print("No legends flag detected, but that feature is not yet implemented")
    if args.ensemble:
        print("Ensemble flag detected, but that feature is not yet implemented")
    if args.plot_budgets:
        print("Budget moments flag detected, but that feature is not yet implemented")
    if args.bu_morr:
        print("Morrison breakdown flag detected, but that feature is not yet implemented")
    if args.diff:
        print("Diff flag detected, but that feature is not yet implemented")

    pyplotgen = PyPlotGen(args.input, args.output, replace=args.replace, les=args.les, cgbest=args.plot_golaz_best,
                          hoc=args.plot_hoc_2005, plotrefs=args.all_best, zip=args.zip, thin=args.thin,
                          no_legends=args.no_legends, ensemble=args.ensemble, budget_moments=args.plot_budgets,
                          bu_morr=args.bu_morr, diff=args.diff)
    return pyplotgen

if __name__ == "__main__":
    pyplotgen = process_args()
    pyplotgen.run()
