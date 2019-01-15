'''
This file is the main driver for PyPlotGen. It handles the overall
process of retrieving the user's input and dictating the loading
of files, generation of plots, and then the export of plots, though
these processes are mostly carried out by other classes/files.

Author: Nicolas Strike
Date: Jan 2019
'''
import argparse

def main():
    '''
    This is the main method that starts the program
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






# run main if file is run by interpreter
if __name__ == "__main__":
    main()