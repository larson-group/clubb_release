#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
bla
TODO: Clean up logger output
"""

#-------------------------------------------------------------------------------
#    IMPORTS
#-------------------------------------------------------------------------------
import sys
import logging
import importlib as ilib
from help.plot_defs import *
from help.plotgen import plot_dict

#-------------------------------------------------------------------------------
#   L O G G E R
#-------------------------------------------------------------------------------
formatter = logging.Formatter('%(asctime)s:%(levelname)s:%(message)s')
fileHandler = logging.FileHandler('plotgen.log', mode='w')
streamHandler = logging.StreamHandler()
fileHandler.setFormatter(formatter)
streamHandler.setFormatter(formatter)
logger = logging.getLogger('plotgen')
logger.addHandler(fileHandler)
logger.addHandler(streamHandler)
#logger.setLevel(logging.INFO)
logger.setLevel(logging.DEBUG)
#logger.setLevel(logging.CRITICAL)


def main():
    
    logger.info('Python plotgen')
    logger.info('--------------')
    
    type_module = None # Variable for plot type name
    case = None # Variable for case name
    
    ## TODO: Work with command line parameters:
    #if len(sys.argv)>1:
    ## Try and get case module name from list:
    #try:
    #module = case_names[int(sys.argv[1])]
    #module = case_files[module]
    #except:
    #logger.warning("Invalid input: {}.".format(sys.argv[1]))
    #module = None
    #else:
    #logger.info("No input case was given.")
    
    ## Load variables:
    for n in range(ntrials):
        logger.info("Receive user input plot type:")
        # List case module names:
        for i,key in enumerate(type_dict):
            logger.info("%d: %s",i+1,key)
        # Try and convert input to list index:
        try:
            # Read number
            user_input = int(raw_input("Specify number -> "))-1
            # Generate module name
            type_module = type_name_template.format(type_modules[user_input])
        except (ValueError, IndexError) as err:
            logger.warning("Given input was not valid. Specify an integer between 1 and %d to select a plot type.\nAttempt %d/%d:", len(type_dict), n, ntrials)
            continue
        break
    if (type_module is None):
        logger.error("Too many invalid attempts. Exiting...")
        sys.exit()
    # Dynamically import the specified type module
    logger.debug("Loading variables module: "+type_module)
    plots = ilib.import_module("..{}".format(type_module),"cases.subpkg")
    
    for n in range(ntrials):
        logger.info("Receive user input case:")
        # List case module names:
        for i,key in enumerate(case_dict):
            logger.info("%d: %s",i+1,key)
        # Try and convert input to list index:
        try:
            user_input = int(raw_input("Specify number -> "))-1
            case = case_modules[user_input]
        except (ValueError, IndexError) as err:
            logger.warning("Given input was not valid. Specify an integer between 1 and %d to select a case.\nAttempt %d/%d:", len(case_dict), n, ntrials)
            continue
        break
    if (case is None):
        logger.error("Too many invalid attempts. Exiting...")
        sys.exit()
    # Dynamically import the specified case module
    logger.debug("Loading case setup module: "+case)
    cf = ilib.import_module("..{}".format(case),"cases.subpkg")
    
    # Call subprogram:
    for plotter in plot_dict:
        if plotter in type_module:
            plot_dict[plotter](plots, cf)

if __name__ == "__main__":
    main()