from help import plot_budgets as pb
from help import OutputWriter as ow
from cases import rico_budgets_case as cf
from cases import general_budget_variables as bv
import numpy as np
import os
import sys
from netCDF4 import Dataset
import logging

#-------------------------------------------------------------------------------
#   L O G G E R
#-------------------------------------------------------------------------------
FORMAT='%(asctime)s:%(levelname)s:%(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger('plotgen')
#logger.setLevel(logging.DEBUG)

#-------------------------------------------------------------------------------
#    F U N C T I O N S
#-------------------------------------------------------------------------------
def isFunction(value):
    logger.info('isFunction')
    isFunc = False
    if '+' in value:
        isFunc = True
    elif '-' in value:
        isFunc = True
    elif '*' in value:
        isFunc = True
    elif '/' in value:
        isFunc = True
    else:
        isFunc = False
    return isFunc

#-------------------------------------------------------------------------------
#    M A I N
#-------------------------------------------------------------------------------
if __name__ == "__main__":
    logger.info('plotgen.py')

    if not os.path.exists(cf.out_dir):
        os.makedirs(cf.out_dir)
    
    if not os.path.exists(cf.sam_file):
        logger.error('The .nc file does not exist.')
        sys.exit("The .nc file does not exist.")
    
    nc = Dataset(cf.sam_file, "r")
    
    logger.info('Read SAM profiles')
    
    # Grap cell altitudes
    
    # get the specific levels
    level = pb.get_budgets_from_nc(nc, 'z', 1.,1,1)
    idx_z0 = (np.abs(level[:] - cf.startHeight)).argmin()
    idx_z1 = (np.abs(level[:] - cf.endHeight)).argmin() +1
    level = level[idx_z0:idx_z1]
    
    # get the specific time interval
    time = pb.get_budgets_from_nc(nc, 'time',1.,1,1)
    idx_t0 = (np.abs(time[:] - cf.startTime)).argmin()
    idx_t1 = (np.abs(time[:] - cf.endTime)).argmin()

    n = len(level)
    t = len(time)
    
    imageName = cf.out_dir + 'jpg/'
    imageNames = []
    
    # grap the data
    for i in range(len(bv.lines)):
        budget = bv.lines[i]
        plot = bv.plotNames[i]
        functions = []
        func_names = []
        budgets_data = []
        logger.debug("plot: %s", str(bv.sortPlots[i]))
        for j in range(len(budget)):
            if not isFunction(budget[j][2]):
            # grap data of each variable that is not a function
                logger.info("Grap data of: %s", budget[j][0])
                value = pb.mean_profiles(pb.get_budgets_from_nc(nc, budget[j][0], budget[j][3], n, t), idx_t0, idx_t1, idx_z0, idx_z1)
                if np.any(value < -100000):
                # if there are no values for the variable
                    value = np.zeros(n)
                    logger.warning("Could not find the variable %s of %s", budget[j][0], bv.sortPlots[i])
                budgets_data.append([budget[j][0], budget[j][1], budget[j][2], value])
            else:
            # save a function for an evaluation
                functions.append(budget[j][2])
                func_names.append(budget[j][0])
        for k in range(len(functions)):
        # evaluate all functions
            function = functions[k]
            logger.debug(func_names[k])
            logger.debug(function)
            for l in range(len(budgets_data)):
                logger.info("Calculate %s", budgets_data[l][0])
                function = function.replace(budgets_data[l][0], "budgets_data["+str(l)+"][3]")
            if function != "":
                logger.debug(function)
                res = eval(function)
                logger.debug(res)
            budgets_data.append([func_names[k], True, func_names[k], res])
        # plot the budget
        name = cf.plot_case_name + bv.sortPlots[i] + '.jpg'
        imageNames.append(name)
        pb.plot_budgets(budgets_data, level, plot[1], cf.yLabel, plot[0], imageName + name, linewidth = cf.lineWidth, color = cf.color)
        
    # write html page
    logger.info("Write HTML page")
    index = cf.out_dir + 'index.html'
    mode = 'Splotgen'    
    ow.writeNavPage(cf.out_dir, cf.headerText)
    ow.writePlotsPage(cf.out_dir, cf.headerText, mode, imageNames)
    ow.writeIndex(index, mode)