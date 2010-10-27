#! /usr/bin/python
# Author: Cavyn VonDeylen
# Date: September 2010
# Larson-Group UWM

import sys  # Handles command line arguments
import re   # Regular expressions
import readBinaryData  # Reads GrADS .dat files

# Modify this to point to the directory containing the output files
FILEPATH = "../../output/"

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
        
    # Find the timestep
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
            print "Checking budgets when t = " + str(t)
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
    from pupynere import NetCDFFile

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
        timestep = 60 * (ncFile.variables['time'].getValue(numIterations-1) - \
                        ncFile.variables['time'].getValue(numIterations-2))
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
            print "Checking budgets when t = " + str(t)
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
    varNum = 0
    rightHandValue = [0]*numLevels # RHS of budget equation. Each z level is treated separately
    
    for line in ctlFile:
            
        # Once a budget term is found, sum each following line that contains the same 
        # variable prefix. Then determine whether the difference between this sum and
        # the budget term is within a passable tolerance
        if findingBudget == True:
        
            # See if line contains same variable prefix as budget term
            componentTerm = re.match("\w+?_", line).group()
            if componentTerm == termName:
            
                # Vars in the budget have descriptions that include eg. "thlm budget:"
                if line.find("budget:") != -1:
                    
                    componentValue = readBinaryData.readGradsData \
                        (FILEPATH + datFileName, numLevels, iteration, iteration, varNum, numVars)
                        
                    # Add componentValue to rightHandValue,
                    # gradually summing up all the component variables
                    rightHandValue = [a + b for a,b in zip(rightHandValue, componentValue)]
        
            # Find error between budget term and sum of component terms
            else:
                errorDifference = [a - b for a,b in zip(leftHandValue, rightHandValue)]
                allowedTolerance = calcTolerance(termUnits, timestep, termName[0:-1])
                percentError = calcPercentError(leftHandValue, rightHandValue, allowedTolerance)

                # Ignore zt(1) since it is below ground. Still use zm(1) however
                if fileName.find("_zt") != -1:
                    zLevel = 1
                    errorDifference = errorDifference[1:]
                else:
                    zLevel = 0
                
                testSuccess = dispError(leftHandValue, rightHandValue, errorDifference, allowedTolerance, zLevel, termName[:-1], \
                    termUnits, percentError, testSuccess)
                    
                findingBudget = False
                
        # Find budget term of the form [variablePrefix]_bt
        if re.match("\w+_bt", line) != None and findingBudget == False:
            termName = re.match("\w+_bt", line).group()[0:-2]
            termUnits = re.search("[[].+[]]", line).group()[1:-1]
            leftHandValue = readBinaryData.readGradsData \
                (FILEPATH + datFileName, numLevels, iteration, iteration, varNum, numVars)

            # Can't do completeness check when iteration is 1
            if iteration != 1:
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
    
    rightHandValue = [0]*numLevels # RHS of budget equation. Each z level is treated separately
    
    # Find a budget variable (_bt)
    for varName in varList:
        try:
            # If varName is not budget var, exception will be thrown and next var checked
            budgetVarName = re.match("\w+_bt", varName).group()
            
            budgetVar = ncFile.variables[varName]
            leftHandValue = budgetVar.getValue(iteration-1) # First index of list is 0
            
            # Can't do completeness check when iteration is 1
            if iteration != 1:
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
                                componentValue = var.getValue(iteration-1)
                                rightHandValue = [a + b for a,b in zip(rightHandValue, componentValue)]
                        
                except AttributeError:
                    pass
                    
            # All component terms have been summed up
            errorDifference = [a - b for a,b in zip(leftHandValue, rightHandValue)]
            allowedTolerance = calcTolerance(budgetVar.units, timestep, budgetVarName[0:-3])
            percentError = calcPercentError(leftHandValue, rightHandValue, allowedTolerance)

            # Ignore zt(1) since it is below ground. Still use zm(1) however
            if ncFile.filename.find("_zt") != -1:
                zLevel = 1
                errorDifference = errorDifference[1:]
            else:
                zLevel = 0

            testSuccess = dispError(leftHandValue, rightHandValue, errorDifference, allowedTolerance, zLevel, budgetVarName[:-3], \
                budgetVar.units, percentError, testSuccess)
                
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

    varAfter = readBinaryData.readGradsData \
        (FILEPATH + datFileName, numLevels, iteration, iteration, varNum, numVars)
    varBefore = readBinaryData.readGradsData \
        (FILEPATH + datFileName, numLevels, iteration-1, iteration-1, varNum, numVars)
    
    rightHandValue = [(a - b) / timestep for a,b in zip(varAfter, varBefore)]
    
    errorDifference = [a - b for a,b in zip(leftHandValue, rightHandValue)]
    allowedTolerance = calcTolerance(termUnits, timestep, termName)
    percentError = calcPercentError(leftHandValue, rightHandValue, allowedTolerance)
    
    # Specify that this is the completeness test when printing to stdout
    termName = termName + " completeness test"
    
    zLevel = 0
    testSuccess = dispError(leftHandValue, rightHandValue, errorDifference, allowedTolerance, zLevel, termName, termUnits, percentError, testSuccess)
                    
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
            
            varAfter = statVar.getValue(iteration-1)  # First index in python list is 0
            varBefore = statVar.getValue(iteration-2) # while iteration starts at t = 1
            
            rightHandValue = [(a - b) / timestep for a,b in zip(varAfter, varBefore)]
            
            errorDifference = [a - b for a,b in zip(leftHandValue, rightHandValue)]
            allowedTolerance = calcTolerance(termUnits, timestep, termName)
            percentError = calcPercentError(leftHandValue, rightHandValue, allowedTolerance)
            
            # Specify that this is the completeness test when printing to stdout
            varName = varName + " completeness test"
            
            zLevel = 0
            testSuccess = dispError(leftHandValue, rightHandValue, errorDifference, allowedTolerance, zLevel, varName, termUnits, percentError, testSuccess)
                    
    return testSuccess
    
#--------------------------------------------------------------------------------------------------
def dispError(leftHandValue, rightHandValue, errorDifference, allowedTolerance, zLevel, termName, termUnits, percentError, testSuccess):
    """
    Find and display any budget failures to stdout
    Input: leftHandValue: list of floats representing the left side of equation being examined
           rightHandValue: list of floats representing the right side of equation being examined
           errorDifference: list of differences between two lists being compared.
           allowedTolerance: Tolerance the errorDifference may be off before considered a fail
           zLevel: Starting z level. For zm grid this should be 0, for zt it should be 1.
           termName: Name of variable that is being tested
           termUnits: Units of term being tested
           percentError: List of percent errors between two lists being compared
           testSuccess: Whether the test is succeeding or failing
    """
        
    for value in errorDifference:
        zLevel += 1
        
        # Check to make sure tolerance is not smaller than a values precision
        valuePrecision = calcPrecision(leftHandValue[zLevel-1], rightHandValue[zLevel-1])
        allowedTolerance = max( abs(allowedTolerance), abs(valuePrecision) )
        
        # Display failure message for each error difference thats greater than the tolerance
        if abs(value) > abs(allowedTolerance):
            testSuccess = False
            print termName + " fails at z-level " + str(zLevel) + \
                " with a difference of " + "%e" % value + " " + termUnits \
                + " and error " + "%.9f" % percentError[zLevel-1] + " %"
                
    return testSuccess
    
#--------------------------------------------------------------------------------------------------
def calcPercentError(experimental, accepted, denominator):
    """
    Finds the percent error using the given inputs
    Input: experimental: Your results in list form
           accepted: The correct values in list form
           denominator: Integer same as accepted unless accepted is extremely small (i.e. 0)
    """
    return [((a - b) / denominator) * 100 for a,b in zip(experimental, accepted)]
    
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
    Nr_tol = 1e-10 # num/kg
    
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
        tol = rt_tol * rt_tol
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
    
    # Check all files in FILEPATH
    if sys.argv[1] == "all":
        testableFiles = []
        
        # Make a list of all files in FILEPATH
        files = os.listdir(FILEPATH)
        
        # Only keep files ending with zt or zm not including .dat or rad files
        for dataFile in files:
            if (dataFile.find("zt.") > -1 or dataFile.find("zm.") > -1) and dataFile.find(".dat") == -1 and dataFile.find("rad") == -1:
                testableFiles.append(dataFile)

        testableFiles.sort()
        
        # Test all remaining files
        for dataFile in testableFiles:
            print "\n----------------------------------------------------\nTesting " + dataFile
            
            # Check if file is NetCDF, otherwise assume GrADS
            if dataFile.find(".nc") != -1 or dataFile.find(".cdf") != -1:
                if checkNetcdfBudgets( dataFile, int(sys.argv[2]) ) == False:
                    testSuccess = False
            else:
                if checkGradsBudgets( dataFile, int(sys.argv[2]) ) == False:
                    testSuccess = False
                    
    # Only check 1 file
    else:
    
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
