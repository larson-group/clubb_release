"""
API Commitment Tests (api_commitment_test.py)

Author: schemenauero
Date: 6/24/14

Edited by: Matt McFadden (mcfadd)
Date: 6/7/19

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


def findFiles(dir):
    """
    Returns a list of every file in a directory.

    This function is adapted from Stack Overflow:
    http://stackoverflow.com/questions/19932130/python-iterate-through-folders-then-subfolders-and-print-filenames-with-path-t

    :param dir: The directory to find the files in
    :return: A list containing all of the files in dir
    """
    result = []
    subdirs = [x[0] for x in os.walk(dir)]
    for subdir in subdirs:
        files = os.walk(subdir).next()[2]
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
            # Ends with .F90? Remove that extension only
            if filename.endswith(".F90"):
                filename = filename.replace(".F90", "")
            elif filename.endswith(".f90"):
                filename = filename.replace(".f90", "")
            # Remove directory tree
            directories = filename.split("/")
            filename = directories[len(directories) - 1]
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
    for bannedName in bannedNames:
        openFile = open(file, "r")
        if any([
            "use " + bannedName + "," in openFile.read(),
            "use " + bannedName + " " in openFile.read(),
            "use " + bannedName + "\n" in openFile.read(),
            "use " + bannedName + "\t" in openFile.read(),
            "use " + bannedName + "!" in openFile.read(),
            "USE " + bannedName + "," in openFile.read(),
            "USE " + bannedName + " " in openFile.read(),
            "USE " + bannedName + "\n" in openFile.read(),
            "USE " + bannedName + "\t" in openFile.read(),
            "USE " + bannedName + "!" in openFile.read(),
        ]):
            foundNames.append("use " + bannedName)

        openFile.close()

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
            print "An API error occurred in ", file, errorDictionary[file]


def printCpuOutput(errorDictionary, excludedDirectories):
    """
    Prints any errors as the respective modules only.

    :param errorDictionary: A list of directories to report an error
    :param excludedDirectories: A list of directories to exclude
    :return: None
    """
    modulesUsedSet = set()
    modulesUsedList = []

    # Check for excluded directories
    for file in errorDictionary:
        if not any(e in file for e in excludedDirectories):
            for error in errorDictionary[file]:
                modulesUsedSet.add(error.split(" ")[1])

    # Add set to a list to be sorted
    for module in modulesUsedSet:
        modulesUsedList.append(module)
    modulesUsedList.sort()

    # print list
    for finalModule in modulesUsedList:
        print finalModule


def main(clubbSrcDirectory, searchDirectory, excludedDirectories, cpuMode):
    """
    Main, see header for description.

    :param clubbSrcDirectory: The clubb source directory
    :param searchDirectory: The directory to search in
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
    searchFiles = findFiles(searchDirectory)

    # Do any of them contain an offending line?
    errorDictionary = {}  # All of the files in the host model and their modules which call CLUBB

    for file in searchFiles:
        foundFiles = checkForStrings(bannedNames, file)
        if len(foundFiles) > 0:
            errorDictionary[file] = foundFiles

    # Print out the errors, if any.
    if len(errorDictionary) > 0:
        if cpuMode:
            printCpuOutput(errorDictionary, excludedDirectories)
        else:
            printErrors(errorDictionary, excludedDirectories)
        sys.exit(1)
    else:
        print "Test was successful."


if __name__ == "__main__":
    # parse the command line arguments
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('clubbSrcDir', help="clubb source directory")
    parser.add_argument('searchDir', help="directory to search in")
    parser.add_argument('-cpu', help="print cpu readable output", dest='cpuMode', required=False, action='store_true')
    parser.add_argument('--exclude-dir', help="list of directories to exclude", dest='excludeDir', required=False,
                        nargs='+', default=[])

    args = parser.parse_args(sys.argv[1:])

    # call the main function
    main(args.clubbSrcDir, args.searchDir, args.excludeDir, args.cpuMode)
