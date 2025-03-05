#!/usr/bin/env/ python3
# -*- coding: utf-8 -*-
"""
API Commitment Tests (api_commitment_test.py)

Author: schemenauero
Date: 6/24/14

Edited by: Matt McFadden (mcfadd)
Date: 6/7/19

Updated to Python 3 and adapted for Jenkins by: Steffen Domke and ChatGPT (OpenAI)
Date: 2/27/25

This script first creates a string array of every file's name (sans .F90) in the source folder of CLUBB.
Next, it recursively searches through given folders for files with lines including
"use anyStringInThatArray" (except "use clubb_api_module").
If found, the script exits with a nonzero exit code and prints the offending file.

Usage:
  python api_commitment_test.py [-h] [-cpu] clubbSrcDir searchDir [--exclude-dir EXCLUDEDIR [EXCLUDEDIR ...]]

  For cpu readable output, use the -cpu flag
"""
import os
import sys
import argparse
import re

# Search pattern for use statements
# We want to find all the occurrences of 'use' followed by a CLUBB non-API module name
# Check that the line is not a comment. i.e. no '!' before 'use' in the same line
# '{}' will in turn be later replaced by each of CLUBB's module names
# Check that the character after the module name is one of ',', '!', any whitespace, or a line break
bannedNamePattern = r'\n[^!\n]*use {}[,\s\n!]'

# Whitelist for file extensions that we want to check
fortranFileExtensions = ('.f', '.f90')

def findFiles(dir, exclude=[]):
    """
    Returns a list of every file in a directory.

    This function is adapted from Stack Overflow:
    http://stackoverflow.com/questions/19932130/python-iterate-through-folders-then-subfolders-and-print-filenames-with-path-t

    :param dir: The directory to find the files in
    :return: A list containing all of the files in dir
    """
    result = []
    # Find all subdirectories within the host model folder that do not have any of the excluded foldernames in their path
    # e.g. if we want to exclude all of ['CLUBB', 'SILHS'], then 'SAM/UTIL/SRC' is fine, but 'SAM/SRC/CLUBB' is not
    subdirs = [x[0] for x in os.walk(dir) if not any([e in x[0]+'/' for e in exclude])]
    # Iterate through subdirectories
    for subdir in subdirs:
        # Get all files within that subfolder
        files = next(os.walk(subdir))[2]
        # If it is an actual file, combine foldername and filename and append to result
        if len(files) > 0:
            for file in files:
                result.append(subdir + "/" + file)
    return result


def arrayOfFileNames(files):
    """
    Returns a list of file names (strings) without any ".F90" on the end.

    :param files: A list of file names
    :return: A list of file names without ".F90"
    """
    result = []
    for filename in files:
        # It's not clubb_api_module is it?
        if not ("clubb_api_module.F90" in filename):
            # Ends with .F90? -> Remove that extension from the name
            if filename.lower().endswith('.f90'):
                filename = filename.rsplit('.',maxsplit=1)[0]
            # Remove directory tree
            filename = filename.rsplit('/', maxsplit=1)[1]
            # Add to return array
            result.append(filename)
    return result


def checkForStrings(bannedNames, file):
    """
    Returns a list of the strings found in the given file.

    :param bannedNames: The list of strings to search for in the file
    :param file: The file to search for the strings in
    :return: A list of the strings found in the file
    """
    foundNames = []
    with open(file, 'r') as f:
        try:
            # Read file content line by line
            content = f.read().lower()
            # Iterate through all CLUBB modules except the API module
            for moduleName in bannedNames:
                # Search for use pattern with moduleName inserted
                res = re.search(bannedNamePattern.format(moduleName), content)
                # If found, append to results
                if res:
                    print("Found: Module ", moduleName, ' in file ', file)
                    foundNames.append("use " + moduleName)
        except Exception as e:
            print('Read error for file ', file)
            print(str(e))
            return foundNames

    return foundNames

def printErrors(errorDictionary, excludedDirectories):
    """
    Prints the given files as errors.

    :param errorDictionary: A list of directories to report an error
    :param excludedDirectories: A list of directories to exclude
    :return: None
    """
    for file in errorDictionary:
        if not any(e in file for e in excludedDirectories):
            print("An API error occurred in ", file, errorDictionary[file])


def printCpuOutput(errorDictionary, excludedDirectories):
    """
    Prints any errors as the respective modules only.

    :param errorDictionary: A list of directories to report an error
    :param excludedDirectories: A list of directories to exclude
    :return: None
    """
    modulesUsedSet = set()

    # Check for excluded directories
    for file in errorDictionary:
        if not any(e in file for e in excludedDirectories):
            for error in errorDictionary[file]:
                modulesUsedSet.add(file+ ": " + error.split(" ")[1])

    # Add set to a list to be sorted
    modulesUsedList = sorted(modulesUsedSet)

    # print list
    for finalModule in modulesUsedList:
        print(finalModule)


def main(clubbSrcDirectory, searchObjs, excludedDirectories, cpuMode):
    """
    Main, see header for description.

    :param clubbSrcDirectory: The clubb source directory
    :param searchObjs: The directories and files to search in
    :param excludedDirectories: The directories to exclude
    :param cpuMode: If true prints output in cpu mode
    :return: None
    """
    for i in range(len(excludedDirectories)):
        excludedDirectories[i] = '/' + excludedDirectories[i] + '/'

    excludedDirectories.append(".git")

    # Find all of the files in the clubb src directory.
    files = findFiles(clubbSrcDirectory)

    # Put them in a string array.
    bannedNames = arrayOfFileNames(files)

    # Find all the files in the search directory.
    searchFiles = []
    for o in searchObjs:
        if os.path.exists(o):
            if os.path.isdir(o):
                searchFiles.extend(findFiles(o, excludedDirectories))
            elif os.path.isfile(o):
                searchFiles.append(o)
            else:
                print('Error. Given path ', o, 'is neither a folder nor a file.' )
        else:
            print('Error. Given path ', o, ' does not exist.')

    # Do any of them contain an offending line?
    errorDictionary = {}  # All of the files in the host model and their modules which call CLUBB

    for file in searchFiles:
        # Only check fortran source code files
        if file.lower().endswith(fortranFileExtensions):
            foundFiles = checkForStrings(bannedNames, file)
            if len(foundFiles) > 0:
                errorDictionary[file] = foundFiles

    # Print out the errors, if any.
    if len(errorDictionary) > 0:
        if cpuMode:
            printCpuOutput(errorDictionary, excludedDirectories)
        else:
            printErrors(errorDictionary, excludedDirectories)
        print("Test failed.")
        sys.exit(1)
    else:
        print("Test was successful.")
        sys.exit(0)


if __name__ == "__main__":
    # parse the command line arguments
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('clubbSrcDir', help="clubb source directory", nargs=1)
    parser.add_argument('searchObjs', help="directory or file to search in", nargs='*')
    parser.add_argument('-cpu', help="print cpu readable output", dest='cpuMode', required=False, action='store_true')
    parser.add_argument('--exclude-dir', help="list of directories to exclude", dest='excludeDir', required=False,
                        nargs='+', default=[])

    args = parser.parse_args(sys.argv[1:])

    # call the main function
    main(args.clubbSrcDir[0], args.searchObjs, args.excludeDir, args.cpuMode)
