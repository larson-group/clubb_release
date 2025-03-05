#!/usr/bin/env/ python3
# -*- coding: utf-8 -*-
"""
Create Module Table (create_module_table.py)

Author: schemenauero
Date: 7/02/14

When ran in a folder which also has module lists, this script returns
an html table of which modules are used in which list
"""

import argparse
import sys
import datetime
from api_commitment_test import findFiles, arrayOfFileNames
from tabulate import tabulate

def getClubbCoreModules(clubbCoreLocation):
    retVals = []
    for mod in arrayOfFileNames(findFiles(clubbCoreLocation)):
        if ".svn-base" not in mod:
            retVals.append(mod)
    return retVals


def main(clubbCoreDir, clubbStandaloneTxtFile, clubbCoreTxtFile, samTxtFile, camTxtFile, wrfTxtFile):
    """
    Main, see header for description.

    :param clubbCoreDir: The CLUBB_core directory
    :param clubbStandaloneTxtFile: Text file containing output of api_commitment_test on CLUBB folder
    :param clubbCoreTxtFile: Text file containing output of api_commitment_test on CLUBB_core folder
    :param samTxtFile: Text file containing output of api_commitment_test on SAM folder
    :param camTxtFile: Text file containing output of api_commitment_test on CAM folder
    :param wrfTxtFile: Text file containing output of api_commitment_test on WRF folder
    :return: None
    """
    # Create data lists
    clubb_standalone = []
    clubb_core = []
    sam = []
    wrf = []
    cam = []

    # Read in modules
    with open(clubbStandaloneTxtFile, 'r') as file:
        clubb_standalone = file.readlines()

    with open(clubbCoreTxtFile, 'r') as file:
        clubb_core = file.readlines()

    with open(samTxtFile, 'r') as file:
        sam = file.readlines()

    with open(camTxtFile, 'r') as file:
        cam = file.readlines()

    with open(wrfTxtFile, 'r') as file:
        wrf = file.readlines()

    # Make a master set and list
    moduleSet = set()
    modules = []

    # Add the modules in clubb core
    for mod in getClubbCoreModules(sys.argv[1]):
        moduleSet.add(str(mod) + "\n")

    modules = list(moduleSet)

    # Sort the master list alphabetically
    modules.sort(key=lambda string: string.lower())

    # Setup some table data
    currentUsers = 0
    usedByClubb = False
    usedBySam = False
    usedByWrf = False
    usedByCam = False
    usedByClubbCore = False
    # Using "sortable" http://www.kryogenix.org/code/browser/sorttable/#symbolsbeforesorting
    print(datetime.datetime.now())
    table = []
    table.append(["Module Name", "# of Users", "CLUBB_core", "CLUBB_standalone", "SAM", "WRF", "CAM"])


    # Create the table
    for module in modules:
        currentUsers = 0
        usedByClubbStandalone = False
        usedBySam = False
        usedByWrf = False
        usedByCam = False
        usedByClubbCore = False

        if module in clubb_standalone:
            usedByClubbStandalone = True
            currentUsers += 1

        if module in sam:
            usedBySam = True
            currentUsers += 1

        if module in wrf:
            usedByWrf = True
            currentUsers += 1

        if module in cam:
            usedByCam = True
            currentUsers += 1

        if module in clubb_core:
            usedByClubbCore = True
            currentUsers += 1

        table_line = [str(module), str(currentUsers), str(usedByClubbCore), str(usedByClubbStandalone), str(usedBySam), str(usedByWrf), str(usedByCam)]
        table.append(table_line)

    print(tabulate(table, headers='firstrow', tablefmt='grid'))

if __name__ == "__main__":
    # parse the command line arguments
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('clubbCoreDir', help="CLUBB_core directory")
    parser.add_argument('clubbStandaloneTxtFile', help="Text file output for CLUBB api_commitment_test")
    parser.add_argument('clubbCoreTxtFile', help="Text file output for CLUBB_core api_commitment_test")
    parser.add_argument('samTxtFile', help="Text file output for SAM api_commitment_test")
    parser.add_argument('camTxtFile', help="Text file output for CAM api_commitment_test")
    parser.add_argument('wrfTxtFile', help="Text file output for WRF api_commitment_test")

    args = parser.parse_args(sys.argv[1:])

    # call the main function
    main(args.clubbCoreDir, args.clubbStandaloneTxtFile, args.clubbCoreTxtFile, args.samTxtFile, args.camTxtFile, args.wrfTxtFile)
