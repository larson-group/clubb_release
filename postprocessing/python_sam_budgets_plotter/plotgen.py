"""
-------------------------------------------------------------------------------
   G E N E R A L   I N F O R M A T I O N
-------------------------------------------------------------------------------
This is the python version of the plotgen perl script used to visualize output
of the SAM and SAM_CLUBB simulations. With this program it is possible to plot
specific budgets for specific cases of SAM. They look like the plotgen output.

Before running plotgen, make sure that the information in the imported cases
files is correct as they are used as input for plotgen.py. They must be located
in a subfolder named 'cases' in the same directory as plotgen.py. The same
is true for the imported help files.

To run this program from the command line, first change the working directory to
the directory containing the plotgen.py file, then type the following line into
the console using any python2 interpreter, and finally hit <Enter>:
python2 plotgen.py

TODO:   - Change setup lists in variable setup files to dictionary with meaningful string keys:
id, title, xlabel, text, textpos, sam lines, clubb lines
"""

#-------------------------------------------------------------------------------
from help import plot_budgets as pb
from help import OutputWriter as ow

# budget plots
#from cases import lba_budgets as cf
#from cases import rico_budgets_case as cf
#from cases import dycoms2_rf02_ds_budgets_case as cf
#from cases import general_budget_variables as bv

# corr and covars plots
#from cases import lba_corrs_covars as cf
#from cases import dycoms2_rf02_corrs_covars as cf
#from cases import general_corrs_covars_variables as bv

# standalone
#from cases import lba_standalone as cf
#from cases import dycoms2_rf02_standalone as cf
#from cases import rico_standalone as cf
#from cases import general_standalone_variables as bv

import numpy as np
import os
import sys
from netCDF4 import Dataset
import logging
import importlib as ilib
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import pyplot as plt
from PyPDF2 import PdfFileWriter as pdfw, PdfFileReader as pdfr
from datetime import datetime as dt

#-------------------------------------------------------------------------------
#   L O G G E R
#-------------------------------------------------------------------------------
FORMAT='%(asctime)s:%(levelname)s:%(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger('plotgen')
logger.setLevel(logging.INFO)
#logger.setLevel(logging.DEBUG)
#logger.setLevel(logging.CRITICAL)

#-------------------------------------------------------------------------------
#    C H O I C E    O F    C A S E S
#-------------------------------------------------------------------------------
ntrials = 3 # maximum number of input trials before closing the program

case_files = {
    "BOMEX" : "bomex_case",
    "DYCOMS_RF01" : "dycoms2_rf01_case",
    "DYCOMS_RF02" : "dycoms2_rf02_case",
    "RICO" : "rico_case",
    "LBA" : "lba_case",
    } # list of cases, append as needed

#case_files = {
    #"BOMEX" : "bomex_budgets_case",
    #"DYCOMS_RF01" : "dycoms2_rf01_budgets_case",
    #"DYCOMS_RF02" : "dycoms2_rf02_ds_budgets_case",
    #"RICO" : "rico_domke_budgets_case",
    #"LBA" : "lba_domke_budgets_case",
    #} # list of cases, append as needed
case_list = case_files.keys()
case_module_names = case_files.values()

variables_template = 'general_{}_variables'
type_files = {
    "budget" : "Budgets",
    "standalone" : "Standalone",
    "corrs_covars" : "Correlations and Covariances (not implemented)",
    } # list of cases, append as needed
type_list = type_files.keys()
type_names = type_files.values()

date_file_format = '%Y%m%d'

#-------------------------------------------------------------------------------
#    F U N C T I O N S
#-------------------------------------------------------------------------------
def isFunction(value):
    logger.info('__isFunction__')
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

def makeDirectory(pathToFile):
    if not os.path.exists(pathToFile):
        os.makedirs(pathToFile)

#-------------------------------------------------------------------------------
#    M A I N
#-------------------------------------------------------------------------------
def main():
    n = 1
    module = None # Variable for case name
    out_type = None
    now = dt.now()
    
    ## Load variables:
    while(n<=ntrials):
        logger.info("Please specify output variables:")
        # List case module names:
        for i,c in enumerate(type_names):
            logger.info("%d: %s",i+1,c)
            # Try and convert input to list index:
        try:
            user_input = int(raw_input("Specify number -> "))-1
        except ValueError as ve:
            n = n+1
            logger.warning("Given input was not valid. Specify an integer between 1 and %d to select a case.\nAttempt %d/%d:",len(type_list),n,ntrials)
            continue
        out_type = type_list[user_input]
        logger.info(out_type)
        module = variables_template.format(out_type)
        break
    if (module is None):
        logger.error("Too many invalid attempts. Exiting...")
        sys.exit()
    # Dynamically import the specified case module
    logger.debug("Loading variables package: "+module)
    bv = ilib.import_module("..{}".format(module),"cases.subpkg")
    
    # If input parameters were given:
    module = None
    n = 1
    if len(sys.argv)>1:
        # Try and get case module name from list:
        try:
            module = case_list[int(sys.argv[1])]
            module = case_files[module]
        except:
            logger.warning("Invalid input: {}.".format(sys.argv[1]))
            module = None
    else:
        logger.info("No input case was given.")
    if module is None: # No case module specified yet
        while(n<=ntrials):
            logger.info("Please specify:")
            # List case module names:
            for n,c in enumerate(case_list):
                logger.info("%d: %s",n+1,c)
            # Try and convert input to list index:
            try:
                user_input = int(raw_input("Specify number -> "))-1
            except ValueError as ve:
                n = n+1
                logger.warning("Given input was not valid. Specify an integer between 1 and %d to select a case.\nAttempt %d/%d:",len(case_list),n,ntrials)
                continue
            module = case_list[user_input]
            logger.info(module)
            module = case_files[module]
            break
        if (module is None):
            logger.error("Too many invalid attempts. Exiting...")
            sys.exit()
    # Dynamically import the specified case module
    logger.debug("Loading case file: "+module)
    cf = ilib.import_module("..{}".format(module),"cases.subpkg")
    
    ## Input missing values into strings
    # output directory, identified by creation date
    logger.debug("Output directory template: "+cf.out_dir)
    out_dir = cf.out_dir.format(date=now.strftime(date_file_format))
    logger.debug("Generate output directory name: "+out_dir)
    # path to output pdf, identified by creation date
    logger.debug("Output pdf template: "+cf.pdf_out)
    pdf_out = os.path.join(out_dir, cf.pdf_out.format(type=out_type, date=now.strftime(date_file_format))+'.pdf')
    logger.debug("Generate output pdf name: "+pdf_out)
    # name for .jpg plot files, identified by creation date and output type, plot name still missing
    logger.debug("Output jpg template: "+cf.plot_case_name)
    plot_case_name = cf.plot_case_name.format(type=out_type, date=now.strftime(date_file_format), plot='{plot}')
    logger.debug("Generate output jpg name template: "+plot_case_name)

    logger.info('plotgen.py')
    logger.info('Plotting output of case %s into %s', cf.case, out_dir)
    
    
    # Create output directories
    makeDirectory(out_dir)
    
    # Create pdf output page
    pdf = PdfPages(pdf_out)
    
    # Create subdirectory for output images
    imageName = out_dir + 'jpg/'
    imageNames = []
    makeDirectory(imageName)
    
    # Exit if no nc file found
    if not os.path.exists(cf.sam_file):
        logger.error('The .nc file does not exist.')
        sys.exit("The .nc file does not exist.")
    # open nc file
    try:
        nc = Dataset(cf.sam_file, "r")
    except IOError as err:
        logger.error(err)
        logger.error("In order to run plotgen.py, for the chosen case, a valid path to a .nc file must be specified in the file cases/<case_name>_budgets_case.py in the variable 'sam_file'.")
    
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
    for i in range(len(bv.lines)):
        budget = bv.lines[i]
        # default title and units for plots
        title = bv.plotNames[i][0]
        units = bv.plotNames[i][1]
        # TODO
        # Structure of functions:
        # Index | Description
        #______________________________________________
        #   0   | 
        #   1   | 
        #   2   | 
        #   3   | 
        #   4   | 
        functions = []
        func_names = []
        # Structure of budgets_data:
        # Index | Description
        #______________________________________________
        #   0   | plot label
        #   1   | switch: show plot
        #   2   | varname / expression
        #   3   | plot axis ()
        #   4   | value of variable
        budgets_data = []
        # collect and separate data variables and function variables
        # calculate mean of data over time -> array containing means for all height levels z0<=h<=z1
        for j in range(len(budget)):
            if not isFunction(budget[j][2]):
            # grap data of each variable that is not a function
                logger.info("Grap data of: %s", budget[j][0])
                value = pb.mean_profiles(pb.get_budgets_from_nc(nc, budget[j][2], budget[j][3], n, t), idx_t0, idx_t1, idx_z0, idx_z1)
                if np.any(np.isnan(value)) or np.any(value <= -1000):
                # if there are no values for the variable
                    value = np.zeros(n)
                    logger.warning("Could not find the variable %s of %s", budget[j][0], bv.sortPlots[i])
                budgets_data.append([budget[j][0], budget[j][1], budget[j][2], budget[j][4], value])
            else:
            # save a function for an evaluation
                #functions.append([budget[j][2]])
                functions.append([budget[j][0], budget[j][1], budget[j][2], budget[j][4]])
                #func_names.append(budget[j][0])
        # calculate values of functions with data values
        for k in range(len(functions)):
        # evaluate all functions
            function = functions[k][2]
            func_name = functions[k][0]
            logger.info("Calculate %s", func_name)
            logger.debug(function)
            for l in range(len(budgets_data)):
                function = function.replace(budgets_data[l][0], "budgets_data["+str(l)+"][4]")
            if function != "":
                res = eval(function)
                #logger.debug(res)
                budgets_data.append([func_name, functions[k][1], func_name, functions[k][3], res])
        # plot the budget
        name = plot_case_name.format(plot=bv.sortPlots[i]) + '.jpg'
        imageNames.append(name)
        # get title and units from nc-file
        #tmp = pb.get_long_name(nc, bv.sortPlots[i])
        #if tmp != "nm":
            #title = tmp
        tmp = pb.get_units(nc, bv.sortPlots[i])
        if tmp != "nm":
            units = tmp
            
        # Create figure with subplots for pdf
        pdf_fig, pdf_ax = plt.subplots(figsize=(6,10))
        #figheight = int(len(bv.lines)/2.+.5)*10
        #pdf_fig, pdf_axes = plt.subplots(int(len(bv.lines)/2.+.5),2,figsize=(15,figheight))
        # If 2nd x axis is needed, create, and put all pdf_ax objects into list
        logger.debug("Switch list for 2nd axis: %s",str([b[3] for b in budgets_data]))
        if any([b[3]==1 for b in budgets_data]):
            logger.debug("Twiny")
            pdf_ax = [pdf_ax, pdf_ax.twiny()]
        else:
            pdf_ax = [pdf_ax]
        pb.plot_budgets(budgets_data, level, units, cf.yLabel, title, imageName + name, linewidth = cf.lineWidth, color = cf.color,pdf=pdf_ax)
        #pdf_fig.tight_layout()
        pdf.savefig(pdf_fig)
        plt.close(pdf_fig)
        
    # write html page
    logger.info("Write HTML page")
    index = out_dir + 'index.html'
    mode = 'Splotgen'
    logger.debug("HTML header template: "+cf.headerText)
    headerText = cf.headerText.format(type=type_files[out_type])
    logger.debug("HTML header: "+headerText)
    ow.writeNavPage(out_dir, headerText)
    ow.writePlotsPage(out_dir, headerText, mode, imageNames)
    ow.writeIndex(index, mode)
    pdf.close()
    ## Add bookmarks to pdf output:
    logger.info("Add bookmarks to pdf output")
    writer = pdfw() # create pdf output class
    with open(pdf_out, 'rb') as inpdf: # open input file
        reader = pdfr(inpdf) # open pdf input class
        writer.appendPagesFromReader(reader) # copy all pages to output
        for i in range(writer.getNumPages()):
            writer.addBookmark(bv.sortPlots[i], i, parent=None) # add bookmark
        writer.setPageMode("/UseOutlines") # make pdf open bookmarks overview
        with open(os.path.join(out_dir,'tmp.pdf'),'wb') as tmp_pdf: # open output pdf file
            writer.write(tmp_pdf) # write pdf content to output file
    # file streams closed
    os.rename(os.path.join(out_dir,'tmp.pdf'), pdf_out) # rename newly created pdf

# run main if file is run by interpreter
if __name__ == "__main__":
    main()