
import pdb

def main():

    from netCDF4 import Dataset
    import numpy as np
    import matplotlib.pyplot as plt
    import os
    import sys

    # Set directory to the location of this script
    os.chdir(os.path.dirname(sys.argv[0]))

    #data = Dataset('../output/bomex_300_zt.nc', "r")

    #dir = '../output/'
    caseName = 'dycoms2_rf01'

    variable = 'T_in_K'

    #listOfFilenames = findTimestepFiles(dir,caseName)
    
    listOfFilenames = ['../output/'+caseName+'_0p1_zt.nc', \
                       '../output/'+caseName+'_0p2_zt.nc', \
                       '../output/'+caseName+'_0p3_zt.nc', \
                       '../output/'+caseName+'_0p5_zt.nc', \
                       '../output/'+caseName+'_1_zt.nc', \
                       '../output/'+caseName+'_2_zt.nc', \
                       '../output/'+caseName+'_3_zt.nc', \
                       '../output/'+caseName+'_5_zt.nc', \
                       '../output/'+caseName+'_10_zt.nc', \
                       '../output/'+caseName+'_20_zt.nc', \
                       '../output/'+caseName+'_30_zt.nc', \
                       '../output/'+caseName+'_60_zt.nc', \
                       '../output/'+caseName+'_100_zt.nc', \
                       '../output/'+caseName+'_150_zt.nc', \
                       '../output/'+caseName+'_300_zt.nc']

    timestepArrayMinus1 = np.array([0.2, 0.3, 0.5, 1, 2, 3, 5, 10, 20, 30, 60, 100, 150, 300])

    numFiles = len(listOfFilenames)

    variableArray = extractVariableFromNetcdfFiles(variable, listOfFilenames, numFiles)

    rmseArray = computeRmseArray(variableArray, numFiles)

    print listOfFilenames

    print variableArray

    print rmseArray

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    plt.title(caseName+' case')
    ax1.set_xlabel('Time step [s]')
    ax1.set_ylabel('RMSE of ' + variable)
    clubb_errors, = ax1.loglog(timestepArrayMinus1, rmseArray,'r.')
    conv_1, = ax1.loglog(timestepArrayMinus1, (rmseArray[0]/timestepArrayMinus1[0])*timestepArrayMinus1,'k')
    plt.legend([clubb_errors, conv_1], ['CLUBB convergence', 'Convergence rate of 1'])
    plt.show()

    print np.polyfit( np.log(timestepArrayMinus1), np.log(rmseArray), 1 )

    print("Finished!")

    #pdb.set_trace()

def computeRmseArray(variableArray, numColumns):

    import numpy as np

    rmseArray = np.zeros( numColumns-1 )

    for idx in range(1,numColumns):

           rmseArray[idx-1] = np.sqrt(( (variableArray[:,idx] - variableArray[:,0])**2 ).mean()) 

    #pdb.set_trace()

    return rmseArray

def extractVariableFromNetcdfFiles(variable, listOfFilenames, numFiles):
    """ Pull one profile per netcdf file and glue them together in an array."""

    from netCDF4 import Dataset
    import numpy as np
    
    # Find dimensions of variables in netcdf file
    data = Dataset( listOfFilenames[0], "r" )
    altitudes = data.variables['altitude']
    numAltitudes = len(altitudes)
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

def findTimestepFiles(dir,caseName):
    """Create a list of netcdf files, with one filename per time step value."""

    import glob

    listOfFilenames = glob.glob(dir + caseName + "*" + "_zt.nc")

    return listOfFilenames

# Standard boilerplate to call the main() function to begin
# the program.
if __name__ == '__main__':
    main()


"""
data2 = Dataset('/home/wedowski/Downloads/gq/flags_true/camclubb711_L30_T600.cam.h0.0097-06-18-84585.nc', "r")
data3 = Dataset('/home/wedowski/Downloads/gq/flags_true/camclubb713_L30_T600.cam.h0.0097-06-18-84585.nc', "r")

# Grid cell altitudes
lev = data.variables['lev'][:]
nlev = len(data.dimensions['lev'])
times = data.variables['time']
ntimes = len(times)

analytic = data.variables["PRCOAN"][0:ntimes,:,0,0]
#analytic = data.variables["PRAOAN"][0:ntimes,:,0,0]
gq = data.variables["PRCOGQ"][0:ntimes,:,0,0]
silhs_ub = data.variables["PRC_UB"][0:ntimes,:,0,0]
silhs = data.variables["PRCO"][0:ntimes,:,0,0]
#silhs = data.variables["PRAO"][0:ntimes,:,0,0]
#silhs_ub = data.variables["PRA_UB"][0:ntimes,:,0,0]


analytic_2 = data2.variables["PRCOAN"][0:ntimes,:,0,0]
gq_2 = data2.variables["PRCOGQ"][0:ntimes,:,0,0]
silhs_ub_2 = data2.variables["PRC_UB"][0:ntimes,:,0,0]
silhs_2 = data2.variables["PRCO"][0:ntimes,:,0,0]

#analytic_2 = data2.variables["PRAOAN"][0:ntimes,:,0,0]
#gq_2 = data2.variables["PRCOGQ"][0:ntimes,:,0,0]
#silhs_ub_2 = data2.variables["PRA_UB"][0:ntimes,:,0,0]
#silhs_2 = data2.variables["PRAO"][0:ntimes,:,0,0]

analytic_3 = data3.variables["PRCOAN"][0:ntimes,:,0,0]
gq_3 = data3.variables["PRCOGQ"][0:ntimes,:,0,0]
silhs_ub_3 = data3.variables["PRC_UB"][0:ntimes,:,0,0]
silhs_3 = data3.variables["PRCO"][0:ntimes,:,0,0]

#analytic_2 = data2.variables["PRAOAN"][0:ntimes,:,0,0]
#gq_2 = data2.variables["PRCOGQ"][0:ntimes,:,0,0]
#silhs_ub_2 = data2.variables["PRA_UB"][0:ntimes,:,0,0]
#silhs_2 = data2.variables["PRAO"][0:ntimes,:,0,0]

blue_label  =  '100'
black_label = '1000'
green_label = '2000'

fig, axarr = plt.subplots(1,4,sharex=True, sharey=True)
fig.suptitle("Autoconversion rate", fontsize=20)
plt.subplots_adjust(left=0.05,right=0.97,wspace=0.1)
#plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)

axarr[0].set_xlabel('Analytic (kg/kg/s)', fontsize=20)
axarr[0].set_ylabel('SILHS (kg/kg/s)', fontsize=20)
axarr[0].plot(0,0,'b.',label=blue_label)
axarr[0].plot(0,0,'k.',label=black_label)
axarr[0].plot(0,0,'g.',label=green_label)
axarr[0].plot(analytic[:][:],gq[:][:],'b.')
axarr[0].plot(analytic_2[:][:],gq_2[:][:],'k.')
axarr[0].plot(analytic_3[:][:],gq_3[:][:],'g.')
axarr[0].legend(numpoints=1, loc='upper left')
axarr[0].ticklabel_format(style='sci', axis='x', scilimits=(-3,3))
axarr[0].ticklabel_format(style='sci', axis='y', scilimits=(-3,3))
axarr[0].plot(axarr[0].get_xlim(),axarr[0].get_xlim(),'r')

axarr[1].set_xlabel('Analytic (kg/kg/s)', fontsize=20)
#axarr[1].set_ylabel('SILHS (kg/kg/s)', fontsize=20)
axarr[1].plot(0,0,'b.',label=blue_label)
axarr[1].plot(analytic[:][:],gq[:][:],'b.')
axarr[1].legend(numpoints=1, loc='upper left')
axarr[1].ticklabel_format(style='sci', axis='x', scilimits=(-3,3))
axarr[1].ticklabel_format(style='sci', axis='y', scilimits=(-3,3))
axarr[1].plot(axarr[1].get_xlim(),axarr[1].get_xlim(),'r')

axarr[2].set_xlabel('Analytic (kg/kg/s)', fontsize=20)
#axarr[2].set_ylabel('SILHS (kg/kg/s)', fontsize=20)
axarr[2].plot(0,0,'k.',label=black_label)
axarr[2].plot(analytic_2[:][:],gq_2[:][:],'k.')
axarr[2].legend(numpoints=1, loc='upper left')
axarr[2].ticklabel_format(style='sci', axis='x', scilimits=(-3,3))
axarr[2].ticklabel_format(style='sci', axis='y', scilimits=(-3,3))
axarr[2].plot(axarr[2].get_xlim(),axarr[2].get_xlim(),'r')

axarr[3].set_xlabel('Analytic (kg/kg/s)', fontsize=20)
#axarr[3].set_ylabel('SILHS (kg/kg/s)', fontsize=20)
axarr[3].plot(0,0,'g.',label=green_label)
axarr[3].plot(analytic_3[:][:],gq_3[:][:],'g.')
axarr[3].legend(numpoints=1, loc='upper left')
axarr[3].ticklabel_format(style='sci', axis='x', scilimits=(-3,3))
axarr[3].ticklabel_format(style='sci', axis='y', scilimits=(-3,3))
axarr[3].plot(axarr[3].get_xlim(),axarr[3].get_xlim(),'r')


plt.show()

data.close()
data2.close()
data3.close()
"""
