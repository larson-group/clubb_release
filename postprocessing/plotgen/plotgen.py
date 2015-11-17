#!/usr/bin/python
###########################################################################
# Plotgen v4.0a
#
# Plotgen4 is written to be identical from a user standpoint to plotgen3, 
# except not requiring the use of Matlab.
# Documentation for plotgen3 is available here until plotgen4 docs are complete:
# http://larson-group.com/twiki/bin/view.pl/Documentation/CarsonDoc/Plotgen3
###########################################################################
#TODO add imports
import sys, os, OutputWriter, importlib

#TODO document these default variables even though they're kinda self-explanatory.

# Location of Script
dirPrefix = os.path.dirname(os.path.realpath(__file__))

# Plotgen Version
VERSION	= 4.0

# Run in Nightly mode: 		Default False
nightly = False

# Difference run: 		Default False
diffRun = False

# Plotgen mode dictates where to find .case files.
#   plotgen (Default)
#   splotgen
#   wrfgen
#   camgen
#   gfdlgen
plotgenMode = "plotgen"

# Data file type:  grads or netcdf
dataFileType = "grads"

# Budget plotting mode:		Default false
plotBudgets = False

# Morrison budget plot mode:	Default false
plotMorrBudgets = False

# Output as .maff file:		Default false
outputAsMaff = False

# Display legend on plots:	Default true
displayLegend = True

# Allow output folder to exist: Default false
overwrite = False

# Ensemble run flag:		Default false
ensembleTuner = False
thinLines = False

# Line style, etc, arrays
lineStyles = ["-", "--", "-.", "-"]
lineColors = [
    "[ 0.000, 0.000, 1.000 ]",    # blue
    "[ 0.000, 1.000, 0.000 ]",    # green
    "[ 1.000, 0.000, 1.000 ]",    # magenta
    "[ 0.471, 0.471, 0.471 ]",    # gray
    "[ 0.969, 0.506, 0.749 ]",    # pink
    "[ 1.000, 0.498, 0.000 ]",    # orange
    "[ 0.651, 0.337, 0.157 ]",    # brown
    "[ 0.894, 0.102, 0.110 ]",    # red
    "[ 1.000, 1.000, 0.200 ]"     # yellow
]
lineWidths = [4.5, 3, 2.5, 1.5]
lineStyleCounter = 0
lineColorCounter = 0
lineWidthCounter = 0

# Set width for budget cases.
lineWidthsBudget = 3

outputIndex = ""

# Plot benchmark lines:		Default False
plotAll = False
plotLes = False
plotBest = False
plotDec = False

# HTML page names
navigationPage = "navigation.html"
indexPage = "index.html"

plotcount = 0
consoleOutput = dirPrefix + "/console_output.py"
inputDirs = []
output = ""
outputTemp = ""
casesExecuted = []
DPI = 300
QUALITY = 100
caseCount = 0

###############################################################################
#  Main function
###############################################################################
def main():
##TODO do stuff
    sys.exit(0)

#CASE.CASE

###############################################################################
#  Runs through all cases in case folder.
###############################################################################
def runCases():
    imgNumber = 0
    casePath = getCasePath()
    #This portion will check if we are in nightly mode.  If not, we'll ignore
    #any case files with _nightly.
    if nightly:
        cases = [f for f in os.listdir(casePath) if f.endswith(".case")]
    else:
        cases = [f for f in os.listdir(casePath) 
                 if f.endswith(".case") and "_nightly" not in f]
    for filename in cases:
        runCase = true
        CASE = imp.load_source('CASE',casePath + filename)  ## Python doesn't like import filename
        if CASE.CASE['name'] in casesExecuted: runCase = false

        if runCase and dataExists() and CASE.CASE['enabled'] != false:
            casesExecuted.append(CASE.CASE['name'])
            if (CASE.CASE['type'] == "standard") or
               (CASE.CASE['type'] == "budget" and plotBudgets) or
               (CASE.CASE['type'] == "morrbudget" and plotMorrBudgets):
                OutputWriter.writeCaseTitle(outputIndex,case['headerText'])
                OutputWriter.writeNavePageCase(navigationPage, 
                                               CASE.CASE['name'], CASE.CASE['headerText'])
                # Output additional text from casefile
                if nightly:
                    if CASE.CASE['nightlyOutput']['subText']: 
                        OutputWriter.writeSubHeader(outputIndex, CASE.CASE['nightlyOutput']['subText'])
                    if CASE.CASE['nightlyOutput']['subHtml']:
                        OutputWriter.writeSubHtml(outputIndex, CASE.CASE['nightlyOutput']['subHtml'])
                else:
                    if CASE.CASE['type'] == "morrbudget":  ##This was in the original code
                        OutputWriter.writeMorrBudgetSubHeader(
                            outputIndex, CASE.CASE['additionalOutput']['subText'])
                    if CASE.CASE['additionalOutput']['subText']: 
                        OutputWriter.writeSubHeader(outputIndex, CASE.CASE['additionalOutput']['subText'])
                    if CASE.CASE['additionalOutput']['subHtml']:
                        OutputWriter.writeSubHtml(outputIndex, CASE.CASE['additionalOutput']['subHtml'])
            

                

###############################################################################
# Executes the PlotCreator.
# Arguments:
#   executeMatlab(String matlabArgs)
#     - matlabArgs: The arguments to be passed to the plot creator
###############################################################################
def executePlot():
#TODO maybe objectify this?  Each plot perhaps will be its own object
    sys.exit(1)

###############################################################################
# Changes to the next line style, width, and color.
# Arguments:
#   None.
###############################################################################
def incrementLineTypes():
    lineColorCounter = ((lineColorCounter + 1) % len(lineColors))
    lineStyleCounter = ((lineStyleCounter + 1) % len(lineStyles))
    lineWidthCounter = ((lineWidthCounter + 1) % len(lineWidths))
### This should work fine, if not, below will work when testing.
#    if lineColorCounter >= len(lineColors): lineColorCounter = 0
#    else: lineColorCounter=lineColorCounter+1
#    if lineStyleCounter >= len(lineStyles): lineStyleCounter = 0
#    else: lineStyleCounter=lineStyleCounter+1
#    if lineWidthCounter >= len(lineWidths): lineWidthCounter = 0
#    else: lineWidthCounter=lineWidthCounter+1

###############################################################################
# Does necessary cleanup code.
# Arguments:
#   None.
###############################################################################
def cleanup():
# TODO clean up all temp files and other cruft.
    sys.exit(1)

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
# To do this, I'll need to translate the case files over.
    return 0

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

###############################################################################
#  User is confused.  Possibly drunk.  Print a help message.
###############################################################################
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
    sys.exit(0)

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
            arg = arg[1:]  
            for curarg in arg:
                if   curarg == 'm': outputAsMaff = True
                elif curarg == 'e': keepEps = True
                elif curarg == 'a':
                    plotLes = True
                    plotBest = True
                    plotDec = True
                elif curarg == 'r': overwrite = True
                elif curarg == 'b': plotBest = True
                elif curarg == 'd': plotDec = True
                elif curarg == 'l': plotLes = True
                elif curarg == 'n': nightly = True
                elif curarg == 'q': highQuality = True
                elif curarg == 'g': dataFileType = "grads"
                elif curarg == 't': dataFileType = "netcdf"
                elif curarg == 'thin': thinLines = True
                elif curarg == 'dfrnce': diffRun = True
                elif curarg == 'ensemble': ensembleTuner = True
                elif curarg == 'nolegend': displayLegend = False
                elif curarg == 's': plotgenMode = "splotgen"
                elif curarg == 'c': plotgenMode = "plotgen"
                elif curarg == 'w': plotgenMode = "wrfgen"
                elif curarg == 'cam': plotgenMode = "camgen"
                elif curarg == 'gfdl': plotgenMode = "gfdlgen"
                elif curarg == 'bu': plotBudgets = True
                elif curarg == 'bumorr': plotMorrBudgets = True
                elif curarg == 'h': HELP_MESSAGE()
                else: HELP_MESSAGE()

    if keepEps and outputAsMaff:
        sys.exit("Argument conflict: Please do not simultaneously" + 
                 "choose to save .eps files (-e) and create a .maff" +
                 "file (-m)\n")
    if len(argList)-numArgs < 2:
        HELP_MESSAGE()
    inputDirs = []
    for dir in argList[numArgs:len(argList)-1]:
        if os.path.isdir(dir):
            inputDirs.append(os.path.join(dir)) 
        else:
            print("Input folder " + dir + " does not exist.")
            sys.exit(1)
    outputDir = argList[len(argList)-1]
    if os.path.isdir(outputDir):
        outputDir = os.path.join(outputDir)
        if overwrite:
            try:
                import shutil
                shutil.rmtree(outputDir)
                os.mkdir(outputDir)
            except Exception, e:
                print(e)
                sys.exit(1)
        else:
            print("Output folder alread exists.  Please use -r to allow overwrite")
            sys.exit(1)
    else:
        try:
            os.mkdir(outputDir)
        except Exception, e:
            print e
            sys.exit(1)

###############################################################################
#  Print the version number of Plotgen
###############################################################################
def VERSION_MESSAGE():
    print('Plotgen version $VERSION, Copyright (C) 2013 Larson Group.\n')


# Now that 
main()
