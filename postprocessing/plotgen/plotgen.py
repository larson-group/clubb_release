#!/usr/bin/python
###########################################################################
# Plotgen v4.0a
#
# Plotgen4 is written to be identical from a user standpoint to plotgen3, 
# elif arg == not requiring the use of Matlab.
# Documentation for plotgen3 is available here until plotgen4 docs are complete:
# http://larson-group.com/twiki/bin/view.pl/Documentation/CarsonDoc/Plotgen3
###########################################################################
#TODO add imports
import sys, os

#TODO document these default variables even though they're kinda self-explanatory.
dirPrefix = "/path/to/script/" # TODO
VERSION	= 4.0
DPI = 300
QUALITY = 100
caseCount = 0
nightly = False
diffRun = False
plotgenMode = "plotgen"
dataFileType = "grads"
plotBudgets = False
plotMorrBudgets = False
inputDirs = []
output = ""
outputTemp = ""
outputAsMaff = False
displayLegend = 1
thinLines = False
ensembleTuner = False
linestyles = ["--", "-", "-.", "-"]
lineColors = [
    "[ 1.000, 0.498, 0.000 ]",    # orange
    "[ 0.216, 0.494, 0.722 ]",    # blue
    "[ 0.596, 0.306, 0.639 ]",    # purple
    "[ 0.471, 0.471, 0.471 ]",    # gray
    "[ 0.969, 0.506, 0.749 ]",    # pink
    "[ 0.302, 0.686, 0.290 ]",    # green
    "[ 0.651, 0.337, 0.157 ]",    # brown
    "[ 0.894, 0.102, 0.110 ]",    # red
    "[ 1.000, 1.000, 0.200 ]"     # yellow
]
lineWidths = [4.5, 3, 2.5, 1.5]
lineWidthsBudget = 3
lineStyleCounter = 0
lineColorCounter = 0
lineWidthCounter = 0
outputIndex = ""
plotAll = False
plotLes = False
plotBest = False
plotDec = False
navigationPage = "navigation.html"
indexPage = "index.html"
plotcount = 0
consoleOutput = dirPrefix + "/console_output.py"
##permissions on file copies needs TODO (possibly?)
casesExecuted = []

###############################################################################
#  Main function
###############################################################################
def main():
##TODO do stuff
    exit(0)

###############################################################################
# Executes the PlotCreator.
# Arguments:
#   executeMatlab(String matlabArgs)
#     - matlabArgs: The arguments to be passed to the plot creator
###############################################################################
def executePlot():
#TODO maybe objectify this?  Each plot perhaps will be its own object
    exit(1)

###############################################################################
# Changes to the next line style, width, and color.
# Arguments:
#   None.
###############################################################################
def incrementLineTypes():
#TODO find where to change put this, since it seems like really slow code as is
# Also, change it over to 
    exit(1)

###############################################################################
# Does necessary cleanup code.
# Arguments:
#   None.
###############################################################################
def cleanup():
# TODO clean up all temp files and other cruft.
    exit(1)

###############################################################################
# Checks all input directories to see if one of them contains the current case.
# Will return true if at least on input folder contains data, otherwise false.
# Arguments:
#   dataExists(Hash CASE)
#     - CASE: The case file to check
###############################################################################
def dataExists():
# TODO take in case file, pull list of plots, iterate through all and find datafile names
# check input folders for data files and if the right data file is there for the case
# If so, return 1 (which then runs the case plots elsewhere) else return 0.
# there may be a quicker way for this to work, eventually.
    exit(1)

###############################################################################
#  Casefile path generator function
###############################################################################
def getCasePath():
    casePath = dirPrefix + "cases"
    if plotgenMode == "plotgen":
        return casePath + "/clubb"
    elif plotgenMode == "splotgen":
        return casePath + "/sam_clubb"
    elif plotgenMode == "wrfgen":
        return casePath + "/wrf"
    elif plotgenMode == "camgen":
        return casePath + "/cam"
    elif plotgenMode == "gfdlgen":
        return casePath + "/gfdl"
    else: sys.exit("Plotgen run mode unknown, please use plotgen.py -h for help")
# This error should not be possible unless code was changed.

def HELP_MESSAGE():
    print("Usage: plotgen [OPTION]... INPUT... OUTPUT\n")
    print("  -c\tPlot CLUBB cases [DEFAULT] (equiv to plotgen)\n")
    print("  -s\tPlot SAM_CLUBB cases (equiv to splotgen)\n")
    print("  -w\tPlot WRF_CLUBB cases\n")
    print("  -cam\tPlot CAM cases\n")
    print("  -gfdl\tPlot GFDL cases\n")
    print("  -r\tIf the output folder already exists, replace the contents\n")
    print("  -l\tPlot LES data for comparison.\n")
    print("  -b\tPlot HOC Best Ever data for comparison.\n")
    print("  -d\tPlot HOC 12/17/2005 data for comparison.\n")
    print("  -a\tSame as -lbd. Plots LES, Best Ever, and 12/17/2005 " + 
          "data for comparison.\n")
    print("  -n\tRuns in nightly mode.\n")
    print("  -q\tOutputs high quality images (does not auto scale).\n")
    print("  -e\tDoes not delete EPS images after conversion.\n")
    print("  -m\tOutputs plots compressed inside a .maff directory.\n")
    print("  -g\tUses GrADS data files. [DEFAULT]\n")
    print("  -t\tUses NetCDF data files.\n")
    print("  -thin\tUses thin solid lines\n")
    print("  -dfrnce\tPerforms a 'difference' plot of two output folders\n")
    print("  -nolegend\tPlot without legends\n")
    print("  -ensemble\tUsed for plotting ensemble tuner runs\n")
    print("  -bu\tUsed to plot standard budget plots\n")
    print("  -bumorr\tUsed to plot Morrison budget plots\n")
    print("  -h\tPrints this help message.\n")
    print("Each option must be seperate, eg -r -a not -ra\n")
    exit(0)

###############################################################################
#  Parse arguments
###############################################################################
def parseArgs(argList):
    numArgs = 0
    if len(argList) == 0:
        helpMessage()
    for arg in argList:
        if arg[0] == "-":
            numArgs += 1
            arg = arg[1:]  ##TODO consider allowing -rbd as well as -r -b -d(add for loop, here, and change all the if arg to if arg[0])
            if arg == 'm':
                if keepEps:
                    sys.exit("Argument conflict: Please do not simultaneously" + 
                             "choose to save .eps files (-e) and create a .maff" +
                             "file (-m)\n")
                outputAsMaff = True
            elif arg == 'a':
                plotLes = True
                plotBest = True
                plotDec = True
            elif arg == 'r': overwrite = True
            elif arg == 'b': plotBest = True
            elif arg == 'd': plotDec = True
            elif arg == 'l': plotLes = True
            elif arg == 'n': nightly = True
            elif arg == 'q': highQuality = True
            elif arg == 'e': keepEps = True
            elif arg == 'g': dataFileType = "grads"
            elif arg == 't': dataFileType = "netcdf"
            elif arg == 'thin': thinLines = True
            elif arg == 'dfrnce': diffRun = True
            elif arg == 'ensemble': ensembleTuner = True
            elif arg == 'nolegend': displayLegend = False
            elif arg == 's': plotgenMode = "splotgen"
            elif arg == 'c': plotgenMode = "plotgen"
            elif arg == 'w': plotgenMode = "wrfgen"
            elif arg == 'cam': plotgenMode = "camgen"
            elif arg == 'gfdl': plotgenMode = "gfdlgen"
            elif arg == 'bu': plotBudgets = True
            elif arg == 'bumorr': plotMorrBudgets = True
            elif arg == 'h': HELP_MESSAGE()
            else: HELP_MESSAGE()

    #TODO make sure all exception handling is good, might need to make folders into absolute paths with os.path.join
    if len(argList)-numArgs < 2:
        HELP_MESSAGE()
    inputDirs = []
    for dir in argList[numArgs:len(argList)-1]:
        if os.path.isdir(dir):
            inputDirs.append(dir) 
        else:
            print("Input folder " + dir + " does not exist.")
            sys.exit(1)
    outputDir = argList[len(argList)-1]
    if os.path.isdir(outputDir):
        if overwrite:
            try:
                import shutil
                shutil.rmtree(outputDir)
                os.mkdir(outputDir)
            except Exception, e:
                print(e)
                exit(1)
        else:
            print("Output folder alread exists.  Please use -r to allow overwrite")
            sys.exit(1)
    else:
        try:
            os.mkdir(outputDir)
        except Exception, e:
            print e
            exit(1)

###############################################################################
#  Print the version number of Plotgen
###############################################################################
def VERSION_MESSAGE():
    print('Plotgen version $VERSION, Copyright (C) 2013 Larson Group.\n')

parseArgs(sys.argv[1:])
main()




