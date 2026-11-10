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
If found, the script exits with status 1 and prints the offending file.
Missing, empty, or unreadable Fortran inputs exit with status 2 instead of passing.

Usage:
  python api_commitment_test.py [-h] [-cpu] clubbSrcDir searchDir [-exclude_dir EXCLUDEDIR [EXCLUDEDIR ...]]

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

def input_error(message):
    """Report an incomplete check on stderr and exit with status 2."""
    print("Error: " + message, file=sys.stderr)
    sys.exit(2)


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
    subdirs = [x[0] for x in os.walk(dir, onerror=lambda error: input_error(str(error))) if not any([e in x[0]+'/' for e in exclude])]
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
        content = f.read().lower()
    for moduleName in bannedNames:
        res = re.search(bannedNamePattern.format(moduleName), content)
        if res:
            print("Found: Module ", moduleName, ' in file ', file)
            foundNames.append("use " + moduleName)

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

    if not os.path.isdir(clubbSrcDirectory):
        input_error("CLUBB source directory does not exist: " + clubbSrcDirectory)
    if not searchObjs:
        input_error("Specify at least one host source path to check.")

    # A missing module list would make every import appear allowed.
    files = [path for path in findFiles(clubbSrcDirectory)
             if path.lower().endswith(fortranFileExtensions)]
    if not files:
        input_error("No Fortran source files found in " + clubbSrcDirectory)

    # Put them in a string array.
    bannedNames = arrayOfFileNames(files)

    # Find all the files in the search directory.
    searchFiles = []
    for o in searchObjs:
        if os.path.isdir(o):
            paths = findFiles(o, excludedDirectories)
        elif os.path.isfile(o):
            paths = [o]
        else:
            input_error("Host source path does not exist or is not a file/directory: " + o)
        paths = [path for path in paths if path.lower().endswith(fortranFileExtensions)]
        if not paths:
            input_error("No Fortran source files to check in " + o)
        searchFiles.extend(paths)

    # Do any of them contain an offending line?
    errorDictionary = {}  # All of the files in the host model and their modules which call CLUBB

    for file in searchFiles:
        # Only check fortran source code files
        if file.lower().endswith(fortranFileExtensions):
            try:
                foundFiles = checkForStrings(bannedNames, file)
            except (OSError, UnicodeError) as error:
                input_error("Cannot read " + file + ": " + str(error))
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
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter, add_help=False, allow_abbrev=False)
    parser.add_argument("-h", "-help", action="help", help="Show this help and exit.")
    parser.add_argument('clubbSrcDir', help="clubb source directory", nargs=1)
    parser.add_argument('searchObjs', help="directory or file to search in", nargs='*')
    parser.add_argument('-cpu', help="print cpu readable output", dest='cpuMode', required=False, action='store_true')
    parser.add_argument('-exclude_dir', help="list of directories to exclude", dest='excludeDir', required=False,
                        nargs='+', default=[])

    args = parser.parse_args(sys.argv[1:])

    # call the main function
    main(args.clubbSrcDir[0], args.searchObjs, args.excludeDir, args.cpuMode)
