# Example call (include trailing slash):
# python check_for_missing_threadprivate.py /path/to/CLUBB_core/ /path/to/SILHS/

"""  A function that checks whether CLUBB-SILHS module variables have been 
declared as threadprivate.  This is necessary to ensure thread safety
in CAM-CLUBB runs that use openmp shared-memory parallelization. """

import pdb  # debugging tools

import re  # regex tools

import os  # for traversing directories

import sys  # for getting cmd line args


# pdb.set_trace()

def extractLine(fileName):
    """  Extract a logical line that declares a real 
    or integer.  The logical line may be continued onto other lines. """
    with open(fileName) as srcFile:
        parts = []  # parts of a line that may have a continuation
        # Are we in a line of code that declares a real or integer variable? Not yet.
        isPartOfDecLine = False
        isPartOfTypeDec = False  # is the current line part of a type declaration? Not yet.
        for line in srcFile:
            line = line.split("!")[0]  # Remove comment from end of line
            if line.isspace() or line == '':
                # The line is only a comment and nothing else; skip
                continue
            if "contains" in line:
                break  # There are no more module variables after "contains"
            if re.search(r"\btype\b", line, re.IGNORECASE) is not None \
                    and re.search(r"\bend\s*\btype\b", line, re.IGNORECASE) is None \
                    and re.search(r"\btype\s*\(", line, re.IGNORECASE) is None:
                # this is the beginning of a type declaration; 
                #   stop search until 'end type'
                isPartOfTypeDec = True
                continue
            if re.search(r"\bend\s*\btype\b", line, re.IGNORECASE) is not None:
                # this is the end of a type declaration; 
                #   resume search for declarations
                isPartOfTypeDec = False
                continue
            if re.search(r"\breal\b", line, re.IGNORECASE) is None \
                    and re.search(r"\binteger\b", line, re.IGNORECASE) is None \
                    and re.search(r"\blogical\b", line, re.IGNORECASE) is None \
                    and re.search(r"\bcharacter\b", line, re.IGNORECASE) is None \
                    and re.search(r"\btype\s*\(\b", line, re.IGNORECASE) is None \
                    and not isPartOfDecLine \
                    or isPartOfTypeDec:
                # this line is not part of a variable declaration, so skip it
                continue
            else:
                # this line is part of a declaration
                if "&" in line:  # then line continues
                    splitLine = line.split("&")
                    parts.append(splitLine[0])
                    isPartOfDecLine = True
                else:  # line ends here
                    yield ''.join(parts) + line  # output line
                    parts = []
                    isPartOfDecLine = False
        if parts:
            yield ''.join(parts)  # return any declarations at end


def extractVarNames(decLine):
    """Given a logical line that declares one or more variables,
    return the variable names in a list."""

    if re.search(r"\bparameter\b(?=.*::)", decLine, re.IGNORECASE):
        return []

    if "::" not in decLine:
        return []
    decLine = decLine.split("::", 1)[1]
    decLine = decLine.split("!")[0]

    # Grab identifiers before '=' or ',' or end
    candidates = re.findall(r"\b([A-Za-z_]\w*)\b\s*(?==|,|$)", decLine)

    # Remove intrinsics/keywords
    blacklist = {"reshape", "kind", "dimension", "public", "private"}
    filtered = [c for c in candidates if c.lower() not in blacklist]

    # Remove false positives:
    #  - things ending in '_c' (your numeric suffixes)
    #  - known dimension symbols like d_var_total
    filtered = [c for c in filtered if not c.endswith("_c")]
    filtered = [c for c in filtered if not c.startswith("_")]
    filtered = [c for c in filtered if c not in {"d_var_total"}]

    return filtered


def findDecVars(fileName):
    """ Given a file, return all declared variables in the file. 
    """
    listOfAllVars = []
    for decLine in extractLine(fileName):
        listOfAllVars.extend(extractVarNames(decLine))
    return listOfAllVars


def findThreadPrivates(fileName):
    """ Given a file, return a list of all variables 
    that are declared threadprivate. """
    with open(fileName) as srcFile:
        srcString = srcFile.read()     # Make file into one big string
        srcString = re.split(r"contains\s*\n", srcString)[0]  # Chop off everything after "contains"
        listOfThreadPrivates = re.findall(r"(?si)!\$omp\s*threadprivate\s*(\(.*?\))", srcString)
        listOfThreadPrivates = [x.strip("()") for x in listOfThreadPrivates]  # Eliminate parentheses
        listOfThreadPrivates = [x.split(",") for x in listOfThreadPrivates]  # Separate variables
        listOfThreadPrivates = sum(listOfThreadPrivates, [])  # Convert nested list to simple list
        listOfThreadPrivates = [re.sub(r"(?i)&\s*!\$omp", "", x) for x in
                                listOfThreadPrivates]  # Remove openmp directives
        listOfThreadPrivates = [x.strip() for x in listOfThreadPrivates]  # Strip off whitespace around variable names
    return listOfThreadPrivates


def findF90Files(dir):
    """ Find all Fortran source files (*.F90) in a directory tree. """
    foundFiles = []
    for subdir, _, files in os.walk(dir):
        for file in files:
            if file.endswith(".F90"):
                foundFiles.append(os.path.join(subdir, file))
    return foundFiles


def main():
    """ Test whether all module variables in a set of directories
    are properly listed as threadprivate.
    """

    print("check_for_missing_threadprivate.py has begun.")

    # get command line args where argv[1] is clubb core path and argv[2] is silhs path
    clubbFileNames = findF90Files(sys.argv[1])
    silhsFileNames = findF90Files(sys.argv[2])
    fileNames = clubbFileNames + silhsFileNames

    # outputFile = open("/home/vlarson/Downloads/threadprivate_output.txt","w+")
    failedFiles = []
    passedTest = True
    for fileName in fileNames:
        # print(fileName)
        listOfVars = findDecVars(fileName)
        # print("listOfVars=")
        # print(listOfVars)
        listOfThreadPrivates = findThreadPrivates(fileName)
        # print("listOfThreadPrivates=")
        # print(listOfThreadPrivates)
        varsWoThreadprivates = list(set(listOfVars) - set(listOfThreadPrivates))
        threadprivatesWoVars = list(set(listOfThreadPrivates) - set(listOfVars))
        if varsWoThreadprivates or threadprivatesWoVars:
            passedTest = False
            failedFiles.append(fileName)
            print(f"File {fileName} failed:")

            if varsWoThreadprivates:
                print("  Missing threadprivate for:")
                for var in varsWoThreadprivates:
                    print(f"    {var}")

            if threadprivatesWoVars:
                print("  Threadprivate without declaration for:")
                for var in threadprivatesWoVars:
                    print(f"    {var}")
        # outputFile.write("%s\n" % fileName)
        # outputFile.write("%s\n" % str(varsWoThreadprivates))
        # outputFile.write("%s\n\n" % str(threadprivatesWoVars))
    # outputFile.close()
    if passedTest:
        print("check_for_missing_threadprivate.py passed.")
    else:
        failed_color = '\033[91m'
        end_color = '\033[0m'
        print(failed_color + "check_for_missing_threadprivate.py did not pass." + end_color)
        print(failed_color + "failedFiles = " + str(failedFiles) + end_color)

    return passedTest


# Standard boilerplate to call the main() function to begin
# the program.
if __name__ == '__main__':
    success = main()
    if success:
        sys.exit(0)
    else:
        sys.exit(1)
