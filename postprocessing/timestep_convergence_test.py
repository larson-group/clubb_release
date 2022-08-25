#!usr/bin/env python3
# To use this script, edit caseName, ztVariable, sfcVariable, and dirRoot. Then type
# python timestep_convergence_test.py
# You will have to run the timestep_convergence_test.bash to create the input files for this script. Make sure clubb is compiled when you do this. After that, move all output files from the timestep_convergence_test.bash in clubb/output into the the location specified in the dirRoot var in this script. 
import pdb

def main():

    from netCDF4 import Dataset
    import numpy as np
    import matplotlib
    import matplotlib.pyplot as plt
    import os
    import sys

    # Set directory to the location of this script
    # Typically the script is run from the CLUBB postprocessing directory
    #os.chdir(os.path.dirname(sys.argv[0]))  # used to work with Canopy
    os.chdir(os.path.dirname(os.path.realpath(__file__))) # works with /usr/bin/python

    #data = Dataset('../output/bomex/bomex_300_sfc.nc', "r")
    #pdb.set_trace() 

    #dir = '../output/'
    caseName = 'bomex'

    ztVariable = 'T_in_K'
    
    sfcVariable = 'lwp'

    #listOfFilenames = findTimestepFiles(dir,caseName)
    dirRoot = '../output/bomex_restarted'
    ztFilenameSuffix = 'zt.nc'
    listOfZtFilenames = constructListOfFilenames(dirRoot,caseName,ztFilenameSuffix)
    print(listOfZtFilenames)

    numFiles = len(listOfZtFilenames)

    sfcFilenameSuffix = 'sfc.nc'
    listOfSfcFilenames = constructListOfFilenames(dirRoot,caseName,sfcFilenameSuffix)

    timestepArrayMinus1 = np.array([0.2, 0.3, 0.5, 1, 2, 3, 5, 10, 20, 30, 60, 100, 150, 300])

    # Find number of vertical levels in netcdf file
    print(listOfZtFilenames[0])
    data = Dataset( listOfZtFilenames[0], "r" )
    altitudes = data.variables['altitude']
    numAltitudes = len(altitudes)
    
    # Find number of time of outputs in each simulation
    data = Dataset( listOfSfcFilenames[0], "r" )    
    times = data.variables['time'][:]
    #variableArray[:,idx] = data.variables[variable][numTimes-1,:,0,0]    

    ztVariableArray = extractZtVariableFromNetcdfFiles(ztVariable, listOfZtFilenames, numFiles, numAltitudes)

    ztRmseArray = computeRmseArray(ztVariableArray, numFiles)

    convergenceExponent = np.polyfit( np.log(timestepArrayMinus1), np.log(ztRmseArray), 1 )[0]
    convergenceExponent = np.array([convergenceExponent])

    print(listOfZtFilenames)

    print(ztVariableArray)

    print(ztRmseArray)

    sfcVariableArray = extractSfcVariableFromNetcdfFiles(sfcVariable, listOfSfcFilenames, numFiles, 1)

    fig = plt.figure(figsize=(8,6))
    matplotlib.rcParams.update({'font.size': 18})
    ax1 = fig.add_subplot(111)
    #plt.title(caseName+' shallow cumulus case')
    plt.title(caseName+' case')
    ax1.set_xlabel('Time step [s]')
    #ax1.set_ylabel('RMSE of Temperature  [K]')
    ax1.set_ylabel('RMSE of ' + ztVariable)
    clubb_errors, = ax1.loglog(timestepArrayMinus1, ztRmseArray,'r.', markersize=14)
    conv_1, = ax1.loglog(timestepArrayMinus1, (ztRmseArray[0]/timestepArrayMinus1[0])*timestepArrayMinus1,'k')
    np.set_printoptions(precision=3)
    ax1.text(0.95, 0.05, 'CLUBB convergence\n exponent\n ='+np.array_str(convergenceExponent), \
        transform=ax1.transAxes, fontsize=16, verticalalignment='bottom', horizontalalignment='right')
    plt.legend([clubb_errors, conv_1], ['CLUBB convergence', 'Convergence rate of 1'], loc=0, prop={'size': 16})
    plt.show()
#
#    pdb.set_trace()

    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111)
    plt.title(caseName+' case')
    ax2.set_xlabel('Time since beginning of simulation [s]')
    ax2.set_ylabel(sfcVariable)
    line = ax2.plot(times, sfcVariableArray)
    #plt.legend(handles=[line])
    plt.show()


    print("Finished!")



def constructListOfFilenames(dirRoot,caseName,filenameSuffix):

    listOfFilenames = [dirRoot+'/'+caseName+'_0p1_'+filenameSuffix, \
                       dirRoot+'/'+caseName+'_0p2_'+filenameSuffix, \
                       dirRoot+'/'+caseName+'_0p3_'+filenameSuffix, \
                       dirRoot+'/'+caseName+'_0p5_'+filenameSuffix, \
                       dirRoot+'/'+caseName+'_1_'+filenameSuffix, \
                       dirRoot+'/'+caseName+'_2_'+filenameSuffix, \
                       dirRoot+'/'+caseName+'_3_'+filenameSuffix, \
                       dirRoot+'/'+caseName+'_5_'+filenameSuffix, \
                       dirRoot+'/'+caseName+'_10_'+filenameSuffix, \
                       dirRoot+'/'+caseName+'_20_'+filenameSuffix, \
                       dirRoot+'/'+caseName+'_30_'+filenameSuffix, \
                       dirRoot+'/'+caseName+'_60_'+filenameSuffix, \
                       dirRoot+'/'+caseName+'_100_'+filenameSuffix, \
                       dirRoot+'/'+caseName+'_150_'+filenameSuffix, \
                       dirRoot+'/'+caseName+'_300_'+filenameSuffix]    

    return listOfFilenames

def computeRmseArray(variableArray, numColumns):

    import numpy as np

    rmseArray = np.zeros( numColumns-1 )

    for idx in range(1,numColumns):

           rmseArray[idx-1] = np.sqrt(( (variableArray[:,idx] - variableArray[:,0])**2 ).mean()) 

    #pdb.set_trace()

    return rmseArray

def extractZtVariableFromNetcdfFiles(variable, listOfFilenames, numFiles, numAltitudes):
    """ Pull one profile per netcdf file and glue them together in an array."""

    from netCDF4 import Dataset
    import numpy as np

    # Find number of time steps in each simulation
    data = Dataset( listOfFilenames[0], "r" )    
    times = data.variables['time']
    numTimes = len(times)    
        
    # Initialize array to hold values of variable
    variableArray = np.zeros( (numAltitudes, numFiles) )
    
    for idx, fileName in enumerate(listOfFilenames):
            data = Dataset(fileName, "r")
            # Extract a profile at the last time step
            variableArray[:,idx] = data.variables[variable][numTimes-1,:,0,0]            
            #pdb.set_trace() 
                                 
    return variableArray

def extractSfcVariableFromNetcdfFiles(variable, listOfFilenames, numFiles, numAltitudes):
    """ Pull one profile per netcdf file and glue them together in an array."""

    from netCDF4 import Dataset
    import numpy as np

    # Find number of time steps in each simulation
    data = Dataset( listOfFilenames[0], "r" )    
    times = data.variables['time']
    numTimes = len(times)    
        
    # Initialize array to hold values of variable
    sfcVariableArray = np.zeros( (numTimes, numFiles) )
    
    for idx, fileName in enumerate(listOfFilenames):
            data = Dataset(fileName, "r")
            # Extract a profile at the last time step
            #pdb.set_trace()
            sfcVariableArray[:,idx] = data.variables[variable][0:numTimes,0,0,0]            
            #pdb.set_trace() 
                                 
    return sfcVariableArray


def findTimestepFiles(dir,caseName):
    """Create a list of netcdf files, with one filename per time step value."""

    import glob

    listOfFilenames = glob.glob(dir + caseName + "*" + "_zt.nc")

    return listOfFilenames

# Standard boilerplate to call the main() function to begin
# the program.
if __name__ == '__main__':
    main()


