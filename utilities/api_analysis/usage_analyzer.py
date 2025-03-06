#!/usr/bin/env/ python3
# -*- coding: utf-8 -*-
"""
Usage Analyzer (usage_analyzer.py)

Author: Orion Schemenauer (schemena@uwm.edu)

The Usage Analyzer compares how three groups of FORTRAN modules use one FORTRAN module
and outputs a html table to standard out. In general, one would want to save the resultant HTML to a file.
This is a passive parser. The code should not be compiled prior to using this tool.

Usage: python usage_analyzer.py path/to/single/module.F90 path/to/group/one/ path/to/group/two/ path/to/group/three/
Example: python usage_analyzer.py ../../src/CLUBB_core/clubb_api_module.F90 ../../../SAM_CLUBB/SRC ../../../WRF_CLUBB/WRF/phys ../../../CAM_CLUBB/models/atm/cam/src > output.html
"""

# Importing regex for parsing, system functions for arguments and operating system functions for file walking
import re, sys, os
# Importing strftime for the system's time
from time import strftime
from tabulate import tabulate
import datetime

# Settings
fontSize = 12  # Controls the general font of the output table
fontSizeReduced = 7  # Controls the font of the parent modules in the table

# Whitelist for file extensions that we want to check
fortranFileExtensions = ('.f', '.f90')

# Represents a module in FORTRAN and contains lists of objects which represent the various
# subroutines, functions, variables, and derived types which are in that module.
# For now, programs are treated as Modules. In the future this could be split up or
# this class could be renamed to Container - or something similar.
class Module:
    # Module Constructor
    # Initializes Variables
    def __init__(self, name='', subroutines=[], functions=[], variables=[], derivedTypes=[], uses=[], startLine=-1, endLine=-1):
        self.name = name                    # The name of the module (module [name])
        self.subroutines = subroutines      # Every subroutine definition found in the module
        self.functions = functions          # Every function definition found in the module
        self.variables = variables          # Every variable found in the module (including those in subroutines)
        self.derivedTypes = derivedTypes    # Every derived type defined (not referenced, that's a variable) in the module
        self.uses = uses                    # every use statment in the module ([0.. ] use someMod, only : someVarsOrFunctions)
        self.startLine = startLine          # The line number containing "module [name]"
        self.endLine = endLine              # The line number with "end module [name]"

    # Some explanation here for those unfamiliar with the following three instance functions:
    # In Python (like many other OOP languages), objects that are defined (rather than the ones like "String"
    # which are built into the language) can still be included in data structures, like sets. Since the 
    # programmatic version of a set follows the rules from sets in general mathematics, no set can contain
    # duplicate objects. In order to enforce this, Python must know how to check if a custom object is equal 
    # to another of the same type. Thus, the two comparator functions were overridden to check equality based
    # on module name. In Python, the set implementation uses a "hash map" (probably a dictionary actually)
    # which must know how to hash custom objects. Therefore, the hash function was overridden.

    # Comparator function overridden for custom object
    def __eq__(self, other):
        if isinstance(other, Module):
            return self.name == other.name
        else:
            return False

    # Negative comparator function overridden for custom object
    def __ne__(self, other):
        return not self.__eq__(other)

    # Hash function overridden for custom object set implementation
    def __hash__(self):
        return hash(self.name)

    # Represent function overridden. Defines how the object represents itself when using "print element".
    def __repr_(self):
        return self.name

# Represents a variable, function, subroutine, or derived type in FORTRAN
# It would have been unnecessary code to make all four sseparateclasses
# in a dynamically typed programming language like Python. If this were a 
# statically typed language, like c++, each element should have its own class.
class Element:
    # (For variables, startLine == endLine, regardless of & continuations)
    # Element Constructor
    # Initializes Variables
    def __init__(self, name='', uses=[], startLine=-1, endLine=-1):
        self.name = name                # The name of the element (subroutine [name])
        self.uses = uses                # The use statements found between the startLine and endLine of the element
        self.startLine = startLine      # The line number where the element begins in the FORTRAN code
        self.endLine = endLine          # The line number where the element ends

    # The following three functions are explained in the Module definition above.
    # Comparator function overridden for custom object
    def __eq__(self, other):
        if isinstance(other, Element):
            return self.name == other.name
        elif isinstance(other, str):
            return self.name == other
        else:
            return False

    # Negative comparator function overridden for custom object
    def __ne__(self, other):
        return not self.__eq__(other)

    # Hash function overridden for custom object set implementation
    def __hash__(self):
        return hash(self.name)

    # Represent function overridden. Defines how the object represents itself when using "print element".
    def __repr_(self):
        return self.name


# Represents a use statement in FORTRAN (use [name], only : [elements])
class Use:
    # Use Constructor
    # Initializes Variables
    def __init__(self, name='', elements=[], parentElement=''):
        self.name = name                        # The name of the use statment (use [name], only : [elements])
        self.elements = elements                # The "used" elements. It is worthwhile to mention here that this is a list of
                                                # Strings, not of Elements. In other words, elementInstance.name == element[n].
        self.parentElement = parentElement      # The name of the function, subroutine, or module which this use is contained in
                                                # (In that order of precedence.)

    # The following three functions are explained in the Module definition above.
    # Comparator function overridden for custom object
    def __eq__(self, other):
        if isinstance(other, Use):
            return self.name == other.name
        elif isinstance(other, str):
            return self.name == other
        else:
            return False

    # Negative comparator function overridden for custom object
    def __ne__(self, other):
        return not self.__eq__(other)

    # Hash function overridden for custom object set implementation
    def __hash__(self):
        return hash(self.name)

    # Represent function overridden. Defines how the object represents itself when using "print element".
    def __repr_(self):
        return self.name


# Given a filename string, returns an array of modules found in that file.
# Eg. parseModulesInFile("clubb_api_module.f90")
# This function is the workhorse of the program. It is responsible for parsing ALL the information
# out of the last three arguments. (Not the first, api, argument).
#
# When the program is ran, although all of the information is retrieved from each file in every
# folder given, only the use statements are compared when all is said and done. However, it is
# important to parse every element and module in every file so that each Use knows where
# it came from.
def parseModulesInFile(filename):
    # Take in the entire file at once (small enough text files to not matter)
    lines = []
    i = 1
    with open(filename, "r") as f:
        try:
            for line in f:
                lines.append(line.strip().lower())
                i = i+1
        except Exception as e:
            print("Read error in ", filename, ", line ", i)
            print("Read error in ", filename)
            print(str(e))

    # Setup some data to be used as the function parses
    retModules = list() # The modules found in the current file. Return Value.

    # Each data type that is currently being parsed is stored temporarily.
    # This is done to eliminate the need of accessor methods in the objects, which could result
    # in nasty "This object is only half complete!" bugs.
    # I will list the corresponding instance variable each temp variable below relates to.
    moduleName = ""                 # module.name
    moduleStart = -1                # module.startLine
    moduleSubroutines = list()      # module.subroutines
    moduleFunctions = list()        # module.functions
    moduleVariables = list()        # module.variables
    moduleUses = list()             # module.uses
    moduleDerivedTypes = list()     # module.derivedTypes

    subroutineName = ""             # element.name
    subroutineUses = list()         # element.uses
    subroutineStart = -1            # element.startLine

    functionName = ""               # element.name
    functionUses = list()           # element.uses
    functionStart = -1              # element.startLine

    # Derived Types can't "use" in FORTRAN
    derivedTypeName = ""            # element.name
    derivedTypeStart = -1           # element.startLine

    # Variables are found and created on the spot

    # Go through every line in the current file
    for n in range(len(lines)):

        lines[n] = lines[n].split("!")[0]  # remove any comments

        # Parsing Use Statements
        # For anyone not familiar with regex, the "^\s*use" means
        # "match any number of spaces at the start of the line
        # that come before the word 'use'"
        if re.match(r"^\s*use", lines[n], re.IGNORECASE):                   # Found a Use statement!
            if functionStart != -1:                                         # If the use is in a function,
                functionUses.append(parseUse(n, lines, functionName))       # add it to the current function's use list
            elif subroutineStart != -1:                                     # If the use is in a subroutine,
                subroutineUses.append(parseUse(n, lines, subroutineName))   # add it to the current subroutine's use list
            else:                                                           # If it's not in a function or subroutine, it must be in a module
                moduleUses.append(parseUse(n, lines, moduleName))           # Add it the current module's use list

        # Parsing Variables
        # This regex is "match any number of spaces at the start of the line that come before a variable type, a space, and not preceding 'function'."
        if re.match(r"^\s*((integer)|(character)|(logical)|(real)|(type(?=(\(|\s*\())))(?!\sfunction)", lines[n], re.IGNORECASE): # Found a Variable!
            if not lines[n-1].split("!")[0].strip().endswith("&"):                              # If the line before doesn't end with a &
                for variable in getContinuousLine(n, lines).split("::")[-1].split(","):         # merge the &'s and split on ","
                    moduleVariables.append(Element(variable.split("=")[0].strip(), None, n, n)) # Take left of "=" and add the variable to the module's variable list

        # Parsing Derived Types
        if re.match(r"^\s*type(?!(\(|\s*\())", lines[n], re.IGNORECASE):  # Oh Look, a Derived Type is declared!
            derivedTypeStart = n                                          # Get the startLine number
            derivedTypeName = lines[n].split(" ")[-1].strip()             # The name is everything on the line after the first space

        if re.match(r"^\s*end\stype\s",lines[n], re.IGNORECASE):          # The derived type is done being declared.
            moduleDerivedTypes.append(Element(derivedTypeName, None, derivedTypeStart, n))  # Add the derived types to the module's derivedType list
            derivedTypeName = ""  # Reset the derived type temp vars
            derivedTypeStart = -1

        # Parsing Functions
        if re.match(r"^\s*((?!end)\b.*)function\s", lines[n], re.IGNORECASE):  # Function being declared
            functionStart = n                                                  # Get the line number
            functionName = lines[n].partition("function")[-1].split("(")[0].split("&")[0].strip()  # the name is everything after "function" and before "("

        if re.match(r"^\s*end\sfunction", lines[n], re.IGNORECASE):  # End Function
            moduleFunctions.append(Element(functionName, functionUses, functionStart, n))  # Add the function to the module's function list
            functionName = ""  # Reset the derived type temp vars
            functionUses = list()
            functionStart = -1

        # Parsing Subroutines, almost identical to functions
        if re.match(r"^\s*subroutine\s", lines[n], re.IGNORECASE):  # Subroutine Start
            subroutineStart = n
            subroutineName = lines[n].partition("subroutine")[-1].split("(")[0].split("&")[0].strip()

        if re.match(r"^\s*end\ssubroutine", lines[n], re.IGNORECASE):  # Subroutine End
            moduleSubroutines.append(Element(subroutineName, subroutineUses, subroutineStart, n))
            subroutineName = ""
            subroutineUses = list()
            subroutineStart = -1

        # Parsing Modules
        if re.match(r"^\s*(module|program)\s", lines[n], re.IGNORECASE):  # Module Start
            moduleStart = n                                               # Get the line number
            moduleName = re.split("module|program", lines[n])[-1].strip()  # The name is what follows "module"

        if re.match(r"^\s*end\s(module|program)\s", lines[n], re.IGNORECASE):  # Module End
            retModules.append(Module(moduleName, moduleSubroutines, moduleFunctions,
                                     moduleVariables, moduleDerivedTypes, moduleUses, moduleStart,
                                     n))            # Add the module to the list of modules being returned
            moduleName = ""                         # Reset the module temp vars
            moduleStart = -1
            moduleSubroutines = list()
            moduleFunctions = list()
            moduleVariables = list()
            moduleUses = list()
            moduleDerivedTypes = list()

    return retModules


# Parses an API (reads the public vars and subroutines)
# Handles the first argument given - the single .f90 file.
def parseApiModule(filename):
    # Take in the entire file at once (small enough text files to not matter)
    lines = [line.strip().lower() for line in open(filename)]
    # Some temp data
    modulePublics = list()
    for n in range(len(lines)):                                                       # Go through every line
        if re.match(r"^\s*public", lines[n], re.IGNORECASE):                          # until a public declaration is found
            currentLine = getContinuousLine(n, lines)                                 # Fold away any &'d lines to one long line
            currentLine = currentLine.replace("public", "").replace(" ", "").strip()  # Remove public and external spaces from the line
            for subroutine in currentLine.split(","):                                 # Make an array around the commas
                modulePublics.append(Element(subroutine, None, n, n))                 # Add each resultant element to the moudule's subroutine list
    # Here, the only thing in the API module is some "subroutines" which are really
    # every variable, function, subroutine, and derived type in the module
    return Module(None, modulePublics, None, None, None, None, 1, 100000)


# Returns the use statement on the indicated line
# for example, the code
#
# 1: module module2 
# 2: use module1, only : var1, var2, &
# 3:                   var3, var4
# 4: end module2
#
# when parsed with parseUse(2, linesInTheCurrentFile, "module2")
# would return a Use object like:
# name = "module1"
# elements = ["var1","var2","var3","var4"]
# parentElement = module2
def parseUse(n, lines, parentElement):
    useName = ""  # The name of the use found on the given line number
    usedElements = list()  # The elements in the use (use [useName], only [usedElements])
    lines[n] = getContinuousLine(n, lines).replace(" ", "")  # Turn any &'s in the line to one long line
    useName = lines[n].split(",only:")[0][3:]  # The name of the use is everything before ",only:"
    for element in lines[n].split(",only:")[-1].split(","):  # Go through each element after ",only:"
        usedElements.append(element.split("=>")[-1])  # Add each element, removing any pointers, to the elements list
    return Use(useName, usedElements, parentElement)  # Return a new Use from the parsed data


# Returns a list of every use in a module and its children elements
# This function iterates through a module's subroutines and functions and returns a 
# list of every use it finds, as well as those in the top level of the module itself.
def getTotalUsesInModule(module):
    retUses = list()  # The list of Uses being returned
    retUses.extend(module.uses)  # Add every use in module.uses
    for subroutine in module.subroutines:  # Go through each subroutine
        retUses.extend(subroutine.uses)  # Add each use in the subroutine
    for function in module.functions:  # Go through each function
        retUses.extend(function.uses)  # Add each use in the function
    return retUses  # Return the uses found


# Reformats lines split using ampersands to one line
# Since FORTRAN can split lines using &, they must be returned to one single long line
# Eg:       This is a really long line that &
#           has been split up in FORTRAN.
# Becomes:  This is a really long line that has been split up in FORTRAN.
# This function is recursive
def getContinuousLine(n, lines):
    if re.match(r"^\s*!", lines[n], re.IGNORECASE) or \
            lines[n].startswith("#") or \
            lines[n].strip() == "":  # ignore comments, empty lines, and #
        lines[n] = "&"
    lines[n] = lines[n].split("!")[0].strip()  # Split off any comments
    if lines[n].endswith("&"):  # If the line ends with a &
        lines[n] = lines[n].partition("&")[0]  # Remove the &
        lines[n] += getContinuousLine(n + 1, lines)  # Add the result of this function on the next line
    return lines[n]  # Return the reformatted line


# Returns an array of every file in a directory
# This function is adapted from Stack Overflow: http://stackoverflow.com/questions/19932130/python-
# iterate-through-folders-then-subfolders-and-print-filenames-with-path-t
def findFiles(dir):
    foundFiles = []  # The returned full files locations and names
    subdirs = [x[0] for x in
               os.walk(dir)]  # Walk through each directory (including the given) and find subdirectories recursively
    for subdir in subdirs:  # Go through each subdirectory (including this one)
        currentFiles = next(os.walk(subdir))[2]  # Add the found files to a list
        if (len(currentFiles) > 0):  # If files were found
            for file in currentFiles:  # go through each file
                foundFiles.append(subdir + "/" + file)  # Add its location and the filename to the returned files list
    return foundFiles  # Return the found files


# Formats all this fancy data into a HTML table
def makeTable(apiModule, samModules, wrfModules, camModules):
    table = []
    table.append(["API Element", "SAM", "WRF", "CAM"])

    # Every element (and use!) in the API
    apiElements = list()
    apiElements.extend(apiModule.subroutines)

    # add host models to a list
    hostModels = list()
    hostModels.append(samModules)
    hostModels.append(wrfModules)
    hostModels.append(camModules)

    for apiElement in apiElements:
       temp_line = []
       if not isinstance(apiElement, str):  # Add it to the table's leftmost column
            temp_line.append(apiElement.name)
       else:
            temp_line.append(apiElement)
       for hostModel in hostModels:  # Go through each hostModel
            elementsFound = set()  # Hold all of the elements found
            for module in hostModel:  # Go through each module
                for use in module.uses:  # Go through each use in each module's top level uses (module.uses)
                    if apiElement in use.elements and (
                            use.name == "clubb_api_module" or use.name == "silhs_api_module"):  # If there's a match
                        elementsFound.add(use.parentElement)  # Add the element to the list of found elements
                for function in module.functions:  # Go through each function in the module
                    for use in function.uses:  # Go though each use in the function
                        if apiElement in use.elements and (
                                use.name == "clubb_api_module" or use.name == "silhs_api_module"):  # If there's a match
                            elementsFound.add(use.parentElement + " in " + module.name)  # Add the element to the list of found elements
                for subroutine in module.subroutines:  # Same as functions
                    for use in subroutine.uses:
                        if apiElement in use.elements and (
                                use.name == "clubb_api_module" or use.name == "silhs_api_module"):
                            elementsFound.add(use.parentElement + " in " + module.name)
            if len(elementsFound) == 0:  # If no elements have been found
                temp_line.append(" ")
            else:  # otherwise
                temp_line.extend(elementsFound)
       table.append(temp_line)
    return table  # Return it


# Prints the -h menu
def printHelpMenu():
    print()
    print("Cross references how one module (like an API) is used by three other groups of modules (Like SAM, WRF, and CAM).")
    print("Example Usage: python thisProgram.py path/to/api.F90 path/to/host1 path/to/host2 path/to/host3 > output.html")
    print()
    print("This program can also display information about a file.")
    print("Example Usage: python thisProgram.py path/to/file.F90")
    print()


# Runs the program
def main():
    # If only one file is given as an arg, print out what is found in it (debugging purposes)
    if len(sys.argv) == 2:
        if not sys.argv[1].endswith(fortranFileExtensions):
            print('Error! Listed file is not a fortran source file.')
            sys.exit(1)
        baseModule = parseModulesInFile(sys.argv[1]).pop()
        print(baseModule.name, baseModule.startLine, baseModule.endLine)
        print("subroutines: ")
        for subroutine in baseModule.subroutines:
            print(subroutine.name, subroutine.startLine, subroutine.endLine)
            for use in subroutine.uses:
                print("    ", use.name, "only: ", ", ".join(use.elements))
        print("functions: ")
        for function in baseModule.functions:
            print(function.name)
            for use in function.uses:
                print("    ", use.name, "only: ", ", ".join(use.elements))
        print("variables: ", ([str(variable.name) for variable in baseModule.variables]))
        print("derivedTypes: ", ([str(derivedType.name) for derivedType in baseModule.derivedTypes]))
        print("uses:")
        for use in baseModule.uses:
            print(use.name, "only: ", ", ".join(use.elements))
    # If the expected arguments are given
    elif len(sys.argv) == 5:
        samModules = list()
        wrfModules = list()
        camModules = list()
        apiModule = parseApiModule(sys.argv[1])  # Parse the API
        for file in findFiles(sys.argv[2]):  # Parse SAM
            if file.lower().endswith(fortranFileExtensions):
                samModules.extend(parseModulesInFile(file))
        for file in findFiles(sys.argv[3]):  # parse WRF
            if file.lower().endswith(fortranFileExtensions):
                wrfModules.extend(parseModulesInFile(file))
        for file in findFiles(sys.argv[4]):  # parse CAM
            if file.lower().endswith(fortranFileExtensions):
                camModules.extend(parseModulesInFile(file))
        print(datetime.datetime.now())
        print(tabulate(makeTable(apiModule, samModules, wrfModules, camModules), headers='firstrow', tablefmt='grid')) # print the table
    else:
        raise IOError("An error occurred processing the given arguments.")


# Python's way of having multiple main methods in a set of files
# Also added the -h check here
if __name__ == "__main__":
    try:
        if sys.argv[1] == "-h":
            printHelpMenu()
        else:
            main()
    except IOError:
        print("An error occurred processing the given arguments.")
