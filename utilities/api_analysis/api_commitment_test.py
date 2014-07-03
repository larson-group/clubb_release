import os, sys, re, fileinput

# ~~~~~~ API Commitment Tests ~~~~~~ schemenauero ~~~~~~ 6/24/14 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# This script first creates a string array of every file's name (sans .F90) in the source folder 
# of CLUBB. Next, it recursively searches through given folders for files with lines including 
# "use anyStringInThatArray" (except "use clubb_api_module"). If found, the script 
# exits with a nonzero exit code and prints the offending file.

#Eg Usage: python thisProgram clubbSrcFolder searchDirectory --exclude-dir=dir1,dir2,etc
# For cpu readable output, use the -cpu flag

# Returns an array of every file in a directory
# This function is adapted from Stack Overflow: http://stackoverflow.com/questions/19932130/python-
# iterate-through-folders-then-subfolders-and-print-filenames-with-path-t
def findFiles(dir):
    r=[]
    subdirs = [x[0] for x in os.walk(dir)]
    for subdir in subdirs:
        files = os.walk(subdir).next()[2]
        if (len(files) > 0):
            for file in files:
                r.append(subdir + "/" + file)
    return r

# Returns a list of file names (strings) without any ".F90" on the end
def arrayOfFileNames(files):
    retFilenames = []
    for filename in files:
        # It's not clubb_api_module is it?
        if not ("clubb_api_module.F90" in filename):
            # Ends with .F90? Remove that extension only
            if (filename.endswith(".F90")):
                filename = filename.replace(".F90", "")
            if (filename.endswith(".f90")):
                filename = filename.replace(".f90", "")
            # Remove directory tree
            directories = filename.split("/")
            filename = directories[len(directories)-1]
            # Add to return array
            retFilenames.append(filename)
    return retFilenames

# Returns array of the strings found in the given file
def checkForStrings(bannedNames, file):
    foundNames = []
    for bannedName in bannedNames:
        openFile = open(file, "r")
        if("use "+bannedName+"," in openFile.read() or # This
        "use "+bannedName+" " in openFile.read() or    # Is
        "use "+bannedName+"\n" in openFile.read() or   # Not
        "use "+bannedName+"\t" in openFile.read() or   # The
        "use "+bannedName+"!" in openFile.read() or    # Best
        "USE "+bannedName+"," in openFile.read() or    # Way
        "USE "+bannedName+" " in openFile.read() or    # To
        "USE "+bannedName+"\n" in openFile.read() or   # Do
        "USE "+bannedName+"\t" in openFile.read() or   # This! 
        "USE "+bannedName+"!" in openFile.read()):
            foundNames.append("use "+bannedName)
        openFile.close()
    return foundNames
    
# Prints the given files as errors
def printErrors(errorDictionary, excludedDirectories):
    for file in errorDictionary:
        if not any(e in file for e in excludedDirectories):
            print "An API error occurred in ", file, errorDictionary[file]

#Prints any errors as the respective modules only
def printCpuOutput(errorDictionary, excludedDirectories):
    modulesUsedSet = set()
    modulesUsedList = []

    #Check for excluded directories
    for file in errorDictionary:
        if not any(e in file for e in excludedDirectories):
            for error in errorDictionary[file]:
                modulesUsedSet.add(error.split(" ")[1])
    
    #Add set to a list to be sorted
    for module in modulesUsedSet:
        modulesUsedList.append(module)
    modulesUsedList.sort()

    #print list
    for finalModule in modulesUsedList:
        print finalModule

# Main, see header for description
def main(clubbSrcDirectory, searchDirectory, excludedDirectories, mode):
    excludedDirectories.append(".svn")
    oldExcludedDirectories = list(excludedDirectories) # have to deep copy the list to not create a for-each loop error
    for directory in oldExcludedDirectories:
        excludedDirectories.remove(directory)
        excludedDirectories.append("/"+directory+"/")

    # Find all of the files in the clubb src directory.
    files = findFiles(clubbSrcDirectory)

    # Put them in a string array.
    bannedNames = arrayOfFileNames(files)

    # Find all the files in the search directory.
    searchFiles = findFiles(searchDirectory)

    # Do any of them contain an offending line?
    errorDictionary = {} # All of the files in the host model and their modules which call CLUBB
    for file in searchFiles:
        foundFiles = checkForStrings(bannedNames, file)
        if (len(foundFiles)>0):
            errorDictionary[file] = foundFiles
    # Print out the errors, if any.
    if (len(errorDictionary) > 0):
        if (mode == "computer"):
            printCpuOutput(errorDictionary, excludedDirectories)
        else:
            printErrors(errorDictionary, excludedDirectories)
        sys.exit(1)
    else :
        print "Test was successful."

# Main Call, handles args
if __name__=="__main__":
    try:
        #help
        if (sys.argv[1] == "-h"):
            print "\nSearches a directory of files for CLUBB references."
            print "Usage: python thisProgram clubbSrcFolder searchDirectory [--exclude-dir=dir1,dir2,etc]\n"
            print "For example, this may be a useful command:"
            print 'python api_commitment_test.py CLUBB/src/CLUBB_core CAM/models/atm/cam/src/ --exclude-dir=clubb"\n'

        #Printing output for computers
        elif (sys.argv[1] == "-cpu"):
            try:
                excludeDirArgs = sys.argv[4].split('=')
                excludedDirectories = excludeDirArgs[len(excludeDirArgs)-1].split(",")
            except Exception:
                excludedDirectories = []
            finally:
                main(sys.argv[2], sys.argv[3], excludedDirectories, "computer")

        #Printing output as a human readbale string
        else:
            try:
                excludeDirArgs = sys.argv[3].split('=')
                excludedDirectories = excludeDirArgs[len(excludeDirArgs)-1].split(",")
            except Exception:
                excludedDirectories = []
            finally:
                main(sys.argv[1], sys.argv[2], excludedDirectories, "")
    except IOError:
        print "Unknown or incorrect arguments. Use -h for example usage."
