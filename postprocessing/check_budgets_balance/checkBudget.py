#!/usr/bin/env python
# Author: Cavyn VonDeylen
# Date: September 2010
# Larson-Group UWM

import sys  # Handles command line arguments
import re   # Regular expressions
import readBinaryData  # Reads GrADS .dat files
from pupynere import NetCDFFile # Read NetCDF .cdf / .nc files
from numpy import * # External library for handling large data sets
from time import strftime

# NOTE: This script contains some kludges to hide budget errors we don't have time to fix, marked with TODO

# Modify this to point to the directory containing the output files
FILEPATH = "../../output/"

# Set this to false to skip the completeness tests
COMPLETENESS_TEST = True

# Scale for calculating budget balance tolerance since we cannot easily access the
# model timestep (dt_main). Completeness tests use the frequency of statistical
# output which is obtained from the output files
TIME_SCALE_DENOMINATOR = 60 # Seconds

# How strict the test should be. Anything above this number will be considered
# a failure.
TEST_LENIENCY = 100 # Percent Error

# Any case name containing these strings will be skipped. Usually done due to time constraints
SKIP_LIST = ["nov11_altocu", \
             "twp_ice", \
             "cloud_feedback", \
             "arm_97", \
             "gabls3_night"]
             
# If more failures than this are reported, the test will exit to avoid a huge failure log.
MAX_FAILURES = 1000
numFails = 0 # Global variable indicating the current number of errors found.

#--------------------------------------------------------------------------------------------------
def checkGradsBudgets(fileName, iteration):
    """
    Checks the balance of CLUBB budgets. This is done by assuming
    the budget term (henceforth: leftHandValue)
    is the sum of the budget term components(RightHandValue).
    For example, vm_bt = vm_ma + vm_gf + vm_cf + vm_ta + vm_f + vm_sdmp.
    The script also assumes that in the .ctl file, the components of a budget
    term are listed immediately following the budget variable.
    
    Requires: ReadBinaryData.py for reading binary Grads ".dat" files
    
    Input: fileName: Name of a GrADS .ctl/.dat file pair. Program automatically applies extensions
           iteration: The iteration to look at when balancing the budgets. 0 for all iterations
    """

    testSuccess = True
    
    # Automatically append proper file extensions
    if fileName.find(".") != -1:
        fileName = fileName[:fileName.find(".")]
    
    ctlFileName = fileName + ".ctl"

    # Extract text from file
    try:
        ctlFile = open(FILEPATH + ctlFileName, "r")
        text = ctlFile.read()
    except IOError:
        sys.stderr.write("Cannot find file " + FILEPATH + ctlFileName + "\n")
        sys.exit(1)
        
    # Find number of z levels (ZDEF)
    try:
        numLevels = re.search("ZDEF\s+[0-9]+", text).group()
        numLevels = int( re.search("[0-9]+", numLevels).group() )
    except AttributeError:
        sys.stderr.write("Error parsing number of z levels from *.ctl file\n")
        sys.exit(1)
        
    # Find the number of iterations (TDEF)
    try:
        numIterations = re.search("TDEF\s+[0-9]+", text).group()
        numIterations = int( re.search("[0-9]+", numIterations).group() )
    except AttributeError:
        sys.stderr.write("Error parsing number of iterations from *.ctl file\n")
        sys.exit(1)
        
    # Find the timestep and convert to seconds
    try:
        timestep = re.search("TDEF.+", text).group()
        timestep = int( re.search("[0-9]+mn", timestep).group()[:-2] ) * 60
    except AttributeError:
        sys.stderr.write("Error parsing timestep from *.ctl file\n")
        sys.exit(1)
    
    # Find number of variables (VARS)
    try:
        numVarsIndx = re.search("VARS\s+[0-9]+", text)
        numVars = int( re.search("[0-9]+", numVarsIndx.group()).group() )
    except AttributeError:
        sys.stderr.write("Error parsing number of variables from *.ctl file\n")
        sys.exit(1)
    
    # Check that user entered a valid iteration
    if iteration > numIterations:
        sys.stderr.write(ctlFileName + " only indicates " + str(numIterations) + " iterations\n")
        return True
        
    # Check the budgets for errors. Loop through all timesteps if iteration <= 0, otherwise
    # just test the specified timestep
    if iteration <= 0:
        for t in range(1,numIterations+1):
            testSuccess = findGradsErrorsAtTimestep(t, ctlFile, fileName, numVarsIndx, numVars, \
                                                   numLevels, timestep, numIterations, testSuccess)
    else:
        testSuccess = findGradsErrorsAtTimestep(iteration, ctlFile, fileName, numVarsIndx, numVars, \
                                           numLevels, timestep, numIterations, testSuccess)

    ctlFile.close()
    
    return testSuccess
    
#--------------------------------------------------------------------------------------------------
def checkNetcdfBudgets(fileName, iteration):
    """
    Checks the balance of CLUBB budgets. This is done by assuming
    the budget term (henceforth: leftHandValue)
    is the sum of the budget term components(RightHandValue).
    For example, vm_bt = vm_ma + vm_gf + vm_cf + vm_ta + vm_f + vm_sdmp.
    
    Requires: pupynere for reading NetCDF binary files
              NumPy to use pupynere
              
    Input: fileName: Name of a NetCDF file.
           iteration: The iteration to look at when balancing the budgets. 0 for all iterations
    """

    testSuccess = True

    # Open File
    try:
        ncFile = NetCDFFile(FILEPATH + fileName, 'r')
    except IOError:
        sys.stderr.write("\nCannot find file " + FILEPATH + fileName + "\n\n")
        sys.exit(1)

    # Find number of z levels
    try:
        numLevels = ncFile.dimensions['altitude']
    except StandardError:
        sys.stderr.write("Error parsing number of z levels\n")
        sys.exit(1)
    
    # Find the number of iterations
    try:
        # Pick a random variable and see how many timesteps it has
        varName = ncFile.variables.keys()[0]
        numIterations = ncFile.variables[varName].shape[0]
    except StandardError:
        sys.stderr.write("Error parsing number of iterations\n")
        sys.exit(1)
        
    # Find the timestep
    try:
        if "seconds" in ncFile.variables['time'].units.lower():
            timestep = (ncFile.variables['time'].getValue(numIterations-1) - \
                        ncFile.variables['time'].getValue(numIterations-2))
        elif "minutes" in ncFile.variables['time'].units.lower():
           timestep = 60 * (ncFile.variables['time'].getValue(numIterations-1) - \
                            ncFile.variables['time'].getValue(numIterations-2))
        else:
            sys.stderr.write("Invalid units for timestep\n")
            sys.exit(1)
    except StandardError:
        sys.stderr.write("Error parsing timestep\n")
        sys.exit(1)

    # Find number of variables (VARS)
    try:
        numVars = ncFile.variables.keys()
        # List includes the 4 dimensions. We don't want these
        numVars = len(numVars) - 4
    except StandardError:
        sys.stderr.write("Error parsing number of variables\n")
        sys.exit(1)
    
    # Check that user entered a valid iteration
    if iteration > numIterations:
        sys.stderr.write(fileName + " only indicates " + str(numIterations) + " iterations\n")
        sys.exit(1)
    
    #Prep for finding budget variables
    varList =  ncFile.variables.keys()
    varList.sort()
    
    # Check the budgets for errors. Loop through all timesteps if iteration <= 0, otherwise
    # just test the specified timestep
    if iteration <= 0:
        for t in range(1,numIterations+1):
            testSuccess = findNetcdfErrorsAtTimestep(t, ncFile, numVars, varList, numLevels, \
                                                     timestep, numIterations, testSuccess)
    else:
        testSuccess = findNetcdfErrorsAtTimestep(iteration, ncFile, numVars, varList, numLevels, \
                                                 timestep, numIterations, testSuccess)

    ncFile.close()
    
    return testSuccess

#--------------------------------------------------------------------------------------------------
def findGradsErrorsAtTimestep(iteration, ctlFile, fileName, numVarsIndx, numVars, \
                         numLevels, timestep, numIterations, testSuccess):
    """
    Looks through a .ctl file for a budget variable (a variable ending in _bt). Once found,
    it finds all variables that are components of the budget variable (variables with the same
    prefix, e.g. thlm_ma and thlm_sdmp for the budget variable thlm_bt) and adds these together.
    If the difference is greater than a given tolerance the failure is output to std output.
    
    Input: iteration: The value of t
           ctlFile: File object representing .ctl file
           fileName: Name of .ctl and .dat file (no extension)
           numVarsIndx: Index in ctlFile right before variable list
           numVars: number of variables in .ctl file
           numLevels: number of z levels
           timestep: time between model iterations
           numIterations: number of iterations
           testSuccess: Whether the test is succeeding or failing
    """
    
    findingBudget = False
    
    datFileName = fileName + ".dat"
        
    # Count lines to budget term
    ctlFile.seek(numVarsIndx.end()) # Move file pointer to start of variable declarations
    varNum = 0 # Used for finding variables in GrADS .dat files
    rightHandValue = [0]*numLevels # RHS of budget equation. Each z level is treated separately
    
    for line in ctlFile:
            
        # Once a budget term is found, sum each following line that contains the same 
        # variable prefix. Then determine whether the difference between this sum and
        # the budget term is within a passable tolerance
        if findingBudget == True:
        
            # See if line contains same variable prefix as budget term
            componentTerm = re.match("\w+?_", line).group()
            if componentTerm == termName:
            
                # Vars in the budget have descriptions that include "budget:"
                if line.find("budget:") != -1:
                    
                    componentValue = array(readBinaryData.readGradsData \
                        (FILEPATH + datFileName, numLevels, iteration, iteration, varNum, numVars))
                        
                    # Add componentValue to rightHandValue,
                    # gradually summing up all the component variables
                    rightHandValue = add(rightHandValue, componentValue)
        
            # Find error between budget term and sum of component terms
            else:
                errorDifference = leftHandValue - rightHandValue
                allowedTolerance = calcTolerance(termUnits, TIME_SCALE_DENOMINATOR, termName[0:-1])

                # Ignore zt(1) since it is below ground. Still use zm(1) however
                if fileName.find("_zt") != -1:
                    zLevel = 1
                    errorDifference = errorDifference[1:]
                else:
                    zLevel = 0
                
                testSuccess = dispError(leftHandValue, rightHandValue, errorDifference, allowedTolerance, iteration, zLevel, termName[:-1], \
                    termUnits, testSuccess)
                    
                findingBudget = False
                
        # Find budget term of the form [variablePrefix]_bt
        if findingBudget == False and re.match("\w+_bt", line) != None:
            termName = re.match("\w+_bt", line).group()[0:-2]
            termUnits = re.search("[[].+[]]", line).group()[1:-1]
            leftHandValue = array(readBinaryData.readGradsData \
                (FILEPATH + datFileName, numLevels, iteration, iteration, varNum, numVars))

            # Can't do completeness check when iteration is 1
            if iteration != 1 and COMPLETENESS_TEST == True:
                if termName[:-1] == "rtm" or termName[:-1] == "thlm": #TODO Ignore completeness test failures except for rtm and thlm. See ticket 153
                    # Check that the budget is consistent with previous and next time iterations
                    testSuccess = checkGradsCompleteness(fileName, numLevels, iteration, numVars, \
                                  termName[:-1], termUnits, timestep, leftHandValue, numVarsIndx.end(), testSuccess)

                
            rightHandValue = [0]*numLevels    # Clear old data
            findingBudget = True
             
        varNum += 1
        
    return testSuccess
    
#--------------------------------------------------------------------------------------------------
def findNetcdfErrorsAtTimestep(iteration, ncFile, numVars, varList, numLevels, timestep, numIterations, testSuccess):
    """
    Looks through a list of variables in alphabetical order for budget variables (ending in _bt)
    
    Input: iteration: The value of t
           ncFile: File object representing NetCDF file
           numVars: number of variables
           varList: A list of all the variables in alphabetical order
           numLevels: number of z levels
           numIterations: number of iterations
           testSuccess: Whether the test is succeeding or failing
    """
    
    rightHandValue = array( [0]*numLevels ) # RHS of budget equation. Each z level is treated separately
    
    # Find a budget variable (_bt)
    for varName in varList:
        try:
            # If varName is not budget var, exception will be thrown and next var checked
            budgetVarName = re.match("\w+_bt", varName).group()
            
            budgetVar = ncFile.variables[varName]
            leftHandValue = budgetVar.getValue(iteration-1) # First index of list is 0
            
            # Can't do completeness check when iteration is 1
            if iteration != 1 and COMPLETENESS_TEST == True:
                if budgetVarName[:-3] == "rtm" or budgetVarName[:-3] == "thlm": #TODO Ignore completeness test failures except for rtm and thlm. See ticket 153
                    # Check that the budget is consistent with previous and next time iterations
                    testSuccess = checkNetcdfCompleteness(ncFile, numLevels, iteration, numVars, \
                                  budgetVarName[:-3], varList, budgetVar.units, timestep, leftHandValue, testSuccess)
            
            # Find components of the budget variable
            for variableName in varList:
                try:
                    varPrefix = re.match("\w+_", variableName).group()
                    if varPrefix == budgetVarName[:-2] and variableName[-2:] != "bt":
                            var = ncFile.variables[variableName]
                            
                            # Vars in the budget have descriptions that include eg. "thlm budget:"
                            if var.long_name.find("budget:") != -1:
                                componentValue = array( var.getValue(iteration-1) )
                                rightHandValue = rightHandValue + componentValue
                        
                except AttributeError:
                    pass
                    
            # All component terms have been summed up
            errorDifference = [a - b for a,b in zip(leftHandValue, rightHandValue)]
            allowedTolerance = calcTolerance(budgetVar.units, TIME_SCALE_DENOMINATOR, budgetVarName[0:-3])

            # Ignore zt(1) since it is below ground. Still use zm(1) however
            if ncFile.filename.find("_zt") != -1:
                zLevel = 1
                errorDifference = errorDifference[1:]
            else:
                zLevel = 0

            testSuccess = dispError(leftHandValue, rightHandValue, errorDifference, allowedTolerance, iteration, zLevel, budgetVarName[:-3], \
                budgetVar.units, testSuccess)
                
            # Clear old data
            rightHandValue = [0]*numLevels 
           
        except AttributeError:
            pass
        
    return testSuccess

#--------------------------------------------------------------------------------------------------
def checkGradsCompleteness(fileName, numLevels, iteration, numVars, \
                           termName, termUnits, timestep, leftHandValue, numVarsIndx, testSuccess):
    """
    Check that the budgets are complete with the equation:
    rtm_bt = (rtm_after - rtm_before) / timestep
    This cannot run when t = 1 as rtm_before does not exist in this case.
    Input: datFileName: Needed to read Grads .dat file
           numLevels:       ""
           iteration:       ""
           numVars:         ""
           termName: Name of variable
           termUnits: Units of variable, used for deciding error tolerance
           timestep: time in model that passes before variables are updated
           leftHandValue: value of budget term
           numVarsIndx: Index of ctlFile where variables start
           testSuccess: Whether the test is succeeding or failing
    """
    
    ctlFileName = fileName + ".ctl"
    datFileName = fileName + ".dat"
    
    ctlFile = open(FILEPATH + ctlFileName, "r")
    
    # Find variable number of the real, non-budget variable
    ctlFile.seek(numVarsIndx)
    varNum = 0
    for line in ctlFile:
        if re.match(termName+"[^_]", line) != None:
            break
        varNum += 1

    varAfter = array(readBinaryData.readGradsData \
        (FILEPATH + datFileName, numLevels, iteration, iteration, varNum, numVars))
    varBefore = array(readBinaryData.readGradsData \
        (FILEPATH + datFileName, numLevels, iteration-1, iteration-1, varNum, numVars))
    
    rightHandValue = (varAfter - varBefore) / timestep
    
    errorDifference = leftHandValue - rightHandValue
    allowedTolerance = calcTolerance(termUnits, timestep, termName)
    
    # Specify that this is the completeness test when printing to stdout
    termName = termName + " completeness test"
    
    zLevel = 0
    testSuccess = dispError(leftHandValue, rightHandValue, errorDifference, allowedTolerance, iteration, zLevel, termName, termUnits, testSuccess)
                    
    return testSuccess
    
#--------------------------------------------------------------------------------------------------
def checkNetcdfCompleteness(ncFile, numLevels, iteration, numVars, \
                              termName, varList, termUnits, timestep, leftHandValue, testSuccess):
    """
    Check that the budgets are complete with the equation (using rtm as an example):
    rtm_bt = (rtm_after - rtm_before) / timestep
    This cannot run when t = 1 as rtm_before does not exist in this case.
    Input: ncFile: The NetCDF data file
           numLevels: Number of z levels in data file
           iteration: Value of t
           numVars: Number of variables in data file
           termName: Name of variable
           varList: List of variables in data file
           termUnits: Units of variable, used for deciding error tolerance
           timestep: time in model that passes before variables are updated
           leftHandValue: value of budget term
           testSuccess: Whether the test is succeeding or failing
    """
    for varName in varList:
        # Look through list of variables for main stats variable related to budget variable
        if termName == varName:
        
            statVar = ncFile.variables[varName]
            
            varAfter = array( statVar.getValue(iteration-1) )  # First index in python list is 0
            varBefore = array( statVar.getValue(iteration-2) ) # while iteration starts at t = 1
            
            rightHandValue = (varAfter - varBefore) / timestep
            
            errorDifference = leftHandValue - rightHandValue
            allowedTolerance = calcTolerance(termUnits, timestep, termName)
            
            # Specify that this is the completeness test when printing to stdout
            varName = varName + " completeness test"
            
            zLevel = 0
            testSuccess = dispError(leftHandValue, rightHandValue, errorDifference, allowedTolerance, iteration, zLevel, varName, termUnits, testSuccess)
                    
    return testSuccess
    
#--------------------------------------------------------------------------------------------------
def dispError(leftHandValue, rightHandValue, errorDifference, allowedTolerance, iteration, zLevel, termName, termUnits, testSuccess):
    """
    Find and display any budget failures to stdout
    Input: leftHandValue: list of floats representing the left side of equation being examined
           rightHandValue: list of floats representing the right side of equation being examined
           errorDifference: list of differences between two lists being compared.
           allowedTolerance: Tolerance the errorDifference may be off before considered a fail
           zLevel: Starting z level. For zm grid this should be 0, for zt it should be 1.
           termName: Name of variable that is being tested
           termUnits: Units of term being tested
           testSuccess: Whether the test is succeeding or failing
    """
    global numFails
    
    for value in errorDifference:
        zLevel += 1
        if(zLevel > 1): #TODO: Hides errors at z level 1. See ticket 360 for more info.
            # Check to make sure tolerance is not smaller than a values precision
            valuePrecision = calcPrecision(leftHandValue[zLevel-1], rightHandValue[zLevel-1])
            allowedTolerance = max( abs(allowedTolerance), abs(valuePrecision) )
            
            percentError = calcPercentError(leftHandValue, rightHandValue, allowedTolerance)
            
            # Display failure message for each error difference thats greater than the tolerance
            if abs(percentError[zLevel-1]) >= TEST_LENIENCY: # [zLevel-1] because array starts at 0
                testSuccess = False
                numFails += 1
                print " ".join([ termName, "fails at t=", \
                    str(iteration), "and z=", str(zLevel), "with a difference of", "%e" % value, \
                    termUnits, "and error", "%.9f" % percentError[zLevel-1], "%" ])
                
                if numFails >= MAX_FAILURES:
                    print "Too many failures: exiting test. (Change MAX_FAILURES variable to view more)"
                    sys.exit(1)
                
    return testSuccess
    
#--------------------------------------------------------------------------------------------------
def calcPercentError(experimental, accepted, denominator):
    """
    Finds the percent error using the given inputs
    Input: experimental: Your results in list form
           accepted: The correct values in list form
           denominator: Integer same as accepted unless accepted is extremely small (i.e. 0)
    """
    return (experimental - accepted) / denominator * 100
    
#--------------------------------------------------------------------------------------------------
def calcTolerance(termUnits, timestep, termName):
    """
    Finds the tolerance above which we consider the budgets do not balance. Which one is chosen
    is dependant on a terms units.
    Input: termUnits: String containting the units of the variable in question.
           timestep: Integer of real time in seconds between model iterations
           termName: Name of variable being looked at. Used for finding tolerance for special cases
    """
    
    # These were copied from constants_clubb.F90
    w_tol = 2e-2 # m/s
    thl_tol = 1e-2 # K
    rt_tol = 1e-8 # kg/kg
    Nr_tol = 1e-5 # num/kg
    
    # Special Cases for tolerance
    if termName == "Ncm":
        tol = 2e2
    elif termName == "Nrm":
        tol = 5e3
    # Get error tolerance based on variable units
    elif termUnits == "(num/kg)/s":                                    # num/kg
        tol = Nr_tol
    elif termUnits == "(kg/kg)/s" or termUnits == "kg kg^{-1} s^{-1}": # kg/kg
        tol = rt_tol
    elif termUnits == "(kg^2)/(kg^2 s)": 					        		 	 	 	 	 # kg^2/kg^2
        tol = rt_tol * rt_tol * 100 # Multiply by 100 because otherwise it's too small for tests to pass
    elif termUnits == "m s^{-2}":                                      # m/s
        tol = w_tol
    elif termUnits == "m^2/s^3":                                       # m^2/s^2
        tol = w_tol * w_tol
    elif termUnits == "m^{3} s^{-4}":                                  # m^3/s^3
        tol = w_tol * w_tol * w_tol
    elif termUnits == "K s^{-1}":                                      # K
        tol = thl_tol
    elif termUnits == "(K^2)/s":                                       # K^2
        tol = thl_tol * thl_tol
    elif termUnits == "(m kg)/(s^2 kg)":                               # m/s kg/kg
        tol = w_tol * rt_tol
    elif termUnits == "(kg K)/(kg s)":                                 # K kg/kg
        tol = thl_tol * rt_tol
    elif termUnits == "(m K)/s^2":                                     # K m/s
        tol = thl_tol * w_tol
    else:
        sys.stderr.write("Error parsing units: " + termUnits + "\nCheck this script's calcTolerance method")
        sys.exit(1)
    
    return tol / timestep
    
#--------------------------------------------------------------------------------------------------
def calcPrecision(value1, value2):
    """
    Finds the smallest precision of a given list of numbers and
    returns the smallest representable number with the same precision.
    
    Input: value1: Floating point number
           value2: A second floating point number
    """
    retVal = 0
    
    # Force use of tolerance instead of precision for allowed tolerance if at least one value is 0.
    # Otherwise we get errors.
    if value1 == 0:
        value1 = 1e-40
    if value2 == 0:
        value2 = 1e-40
    
    # Convert the number in scientific notation to a string,
    # then substring the last three digits to obtain the precision.
    valueStr1 = '%e' %(value1)
    valuePrecision1 = int(valueStr1[-3:])
        
    valueStr2 = '%e' %(value2)
    valuePrecision2 = int(valueStr2[-3:])
    
    retVal = min(valuePrecision1, valuePrecision2)
    
    return 10**(retVal)
    
#--------------------------------------------------------------------------------------------------
# Allows this module to be run as a script
if __name__ == "__main__":

    import sys
    import os
    
    # If wrong arguments were given, print a helpful message
    if len(sys.argv) != 3:
        print "Program Help:"
        print "Arguments must be: filename iteration"
        print ""
        print "Filename can be either a NetCDF file (.nc or .cdf) or grads files (.ctl/.dat)"
        print "Set filename to 'all' for all files in the directory specified by the FILEPATH variable"
        print "Set iteration to '0' for all iterations"
        print ""
        print "This script verifies that the budget variable equals the sum of the budget terms"
        print "e.g. vm_bt = vm_ma + vm_gf + vm_cf +...+ vm_sf"
        print "The completeness test verifies that the budgets are consistent"
        print "e.g. rtm_bt(t=5) = ( rtm(t=6) - rtm(t=5) ) / timestep"
        print ""
        sys.exit(1)
        
    testSuccess = True
    skipCase = False
    
    # Check all files in FILEPATH
    if sys.argv[1] == "all":
        testableFiles = []
        
        # Make a list of all files in FILEPATH
        files = os.listdir(FILEPATH)
        
        for dataFile in files:
            for case in SKIP_LIST:
                if dataFile.find(case) > -1:
                    skipCase = True
            
            # Only keep files ending with zt or zm not including .dat or rad files
            if (dataFile.find("zt.") > -1 or dataFile.find("zm.") > -1) and dataFile.find(".dat") == -1 and dataFile.find("rad") == -1 and skipCase == False:
                testableFiles.append(dataFile)
                
            skipCase = False

        testableFiles.sort()
        
        # Make sure data exists
        if len(testableFiles) == 0:
            print "Unable to find testable data"
            testSuccess = False
        
        # Test all remaining files
        for dataFile in testableFiles:
            print "".join(["\n", strftime("%H:%M:%S"), " - Testing ", dataFile])
            
            # Check if file is NetCDF, otherwise assume GrADS
            if dataFile.find(".nc") != -1 or dataFile.find(".cdf") != -1:
                if checkNetcdfBudgets( dataFile, int(sys.argv[2]) ) == False:
                    testSuccess = False
            else:
                if checkGradsBudgets( dataFile, int(sys.argv[2]) ) == False:
                    testSuccess = False
                    
    # Only check 1 file
    else:
        print "\nTesting " +  sys.argv[1]
        
        # Check if file is NetCDF, otherwise assume GrADS
        if sys.argv[1].find(".nc") != -1 or sys.argv[1].find(".cdf") != -1:
            if checkNetcdfBudgets( sys.argv[1], int(sys.argv[2]) ) == False:
                testSuccess = False
        else:
            if checkGradsBudgets( sys.argv[1], int(sys.argv[2]) ) == False:
                testSuccess = False
            
    # Print resolution to the screen
    if testSuccess == True:
        print "Budgets successfully balance!"
    else:
        print "Budgets fail to balance"
