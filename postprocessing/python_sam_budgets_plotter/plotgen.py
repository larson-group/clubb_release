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
    #logger.debug(str(value) + " is function? " + str(isFunc))
    return isFunc

def evalFunction(function, variable_names, variables, place):
    logger.info('evalFunction')
    value = 0
    if function != "":        
        for j in range(len(variable_names)):
            function = function.replace(variable_names[j], "variables["+str(j)+"][" + str(place) + "]")
        logger.debug(function)
        value = eval(function)
    return value

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
    
    # grap the data
    budgets_data = []
    for plot in bv.plotNames:
        logger.debug("plot: %s", str(plot[1]))
        functions = []
        func_names = []
        for i in range(len(bv.lines)):
            if bv.lines[i][0] == plot[0] and not isFunction(bv.lines[i][3]):
                logger.info("Grap data of: %s", bv.lines[i][1])
                value = pb.mean_profiles(pb.get_budgets_from_nc(nc, bv.lines[i][1], bv.lines[i][4], n, t), idx_t0, idx_t1, idx_z0, idx_z1)
                if np.any(value < -100000):
                    value = np.zeros(n)
                budgets_data.append([bv.lines[i][0], bv.lines[i][1], bv.lines[i][2], bv.lines[i][3], value])
            elif bv.lines[i][0] == plot[0] and isFunction(bv.lines[i][3]):
                functions.append(bv.lines[i][3])
                func_names.append(bv.lines[i][1])
        for k in range(len(functions)):
            function = functions[k]
            logger.debug(func_names[k])
            logger.debug(function)
            for j in range(len(budgets_data)):
                if budgets_data[j][0] == plot[0]:
                    logger.info("Calculate %s", budgets_data[j][1])
                    function = function.replace(budgets_data[j][1], "budgets_data["+str(j)+"][4]")
            if function != "":
                logger.debug(function)
                res = eval(function)
                logger.debug(res)
            budgets_data.append([plot[0], func_names[k], True, func_names[k], res])
    
    pb.plot_many_budgets(budgets_data, bv.plotNames, level, cf.yLabel, cf.out_dir + 'jpg/' + cf.plot_case_name, cf.lineWidth)
    
    logger.info("Write HTML page")
    index = cf.out_dir + 'index.html'
    mode = 'Splotgen'    
    ow.writeNavPage(cf.out_dir, cf.headerText)
    ow.writePlotsPage(cf.out_dir, cf.headerText, mode)
    ow.writeIndex(index, mode)