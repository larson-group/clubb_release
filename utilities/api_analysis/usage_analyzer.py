import re, sys, os
from time import strftime

fontSize = 12
fontSizeReduced = 7

# Represents a module in FORTRAN and contains lists of elements which represent the various
# subroutines, functions, variables, and derived types which are in that module
class Module:
    name = ""
    subroutines = list() # Element list
    functions = list()
    variables = list()
    derivedTypes = list()
    uses = list() # Use list
    
    startLine = -1
    endLine = -1
    def __init__(self, name, subroutines, functions, variables, derivedTypes, uses, startLine, endLine):
        self.name         = name
        self.subroutines  = subroutines
        self.functions    = functions
        self.variables    = variables
        self.derivedTypes = derivedTypes
        self.uses         = uses
        self.startLine    = startLine
        self.endLine      = endLine
    def __eq__(self, other):
        return self.name == other.name
    def __ne__(self, other):
        return not self.__eq__(other)
    def __hash__(self):
        return hash(self.name)

# Represents a variable, function, subroutine, or derived type in FORTRAN
class Element:
    name = "" # Both element and use must have a "name" attribute for comparisons
    uses = list()
    startLine = -1
    endLine = -1
    def __init__(self, name, uses, startLine, endLine):
        self.name      = name
        self.uses      = uses
        self.startLine = startLine
        self.endLine   = endLine
    def __eq__(self, other):
        if isinstance(other, Element):
            return self.name == other.name
        else:
            return self.name == other
    def __ne__(self, other):
        return not self.__eq__(other)
    def __hash__(self):
        return hash(self.name)
    def __str__(self):
        return self.name
    def __repr_(self):
        return self.name
        
class Use:
    name = ""
    elements = list() # Strings
    parentElement = ""
    def __init__(self, name, elements, parentElement):
        self.name = name
        self.elements = elements
        self.parentElement = parentElement
    def __eq__(self, other):
        if isinstance(other, Use):
            return self.name == other.name #and (self.elements.difference_update(other.elements) == None)
        else:
            return self.name == other
    def __ne__(self, other):
        return not self.__eq__(other)
    def __hash__(self):
        return hash(self.name)
    def __str__(self):
        return self.name
    def __repr_(self):
        return self.name


# Given a filename string, returns an array of populated Modules
def parseModulesInFile(filename):
    
    # Take in the entire file at once (small enough text files to not matter)
    lines = [line.strip().lower() for line in open(filename)]
    
    retModules = list() # The modules found in the current file
    
    # Each data type that is currently being parsed is stored temporarily
    moduleName = ""
    moduleStart = -1
    moduleEnd = -1
    moduleSubroutines = list()
    moduleFunctions = list()
    moduleVariables = list()
    moduleUses = list()
    moduleDerivedTypes = list()
    
    subroutineName = ""
    subroutineUses = list()
    subroutineStart = -1
    subroutineEnd = -1
    
    functionName = ""
    functionUses = list()
    functionStart = -1
    functionEnd = -1
    
    derivedTypeName = ""
    derivedTypeStart = -1
    derivedTypeEnd = -1
    
    for n in range(0,len(lines)):
        lines[n] = lines[n].split("!")[0] # Split off any comments
        
        if re.match(r"^\s*use", lines[n], re.IGNORECASE): # Use statement
            if (functionStart != -1 and functionEnd == -1): # Used in a function
                functionUses.append(parseUse(n, lines, functionName))
            elif (subroutineStart != -1 and subroutineEnd == -1): # Used in a subroutine
                subroutineUses.append(parseUse(n, lines, subroutineName))
            else: # Used in a module
                moduleUses.append(parseUse(n, lines, moduleName))
                
        if re.match(r"^\s*((integer)|(character)|(logical)|(real)|(type(?=(\(|\s*\())))(?!\sfunction)", lines[n], re.IGNORECASE): # Variable
            if not lines[n-1].split("!")[0].strip().endswith("&"): # If the line before doesn't end with a &
                for variable in getContinuousLine(n, lines).split("::")[-1].split(","): # merge the &'s and split on ","
                    moduleVariables.append(Element(variable.split("=")[0].strip(), None, n, n)) # Take left of "=" and add
        
        if re.match(r"^\s*type(?!(\(|\s*\())", lines[n], re.IGNORECASE): # Derived Type Start
            derivedTypeStart = n
            derivedTypeName = lines[n].split(" ")[-1].strip()
            
        if re.match(r"^\s*end\stype\s",lines[n], re.IGNORECASE): # Derived Type End
            derivedTypeEnd = n
            moduleDerivedTypes.append(Element(derivedTypeName, None, derivedTypeStart, derivedTypeEnd))
            derivedTypeName = ""
            derivedTypeStart = -1
            derivedTypeEnd = -1
        
        if re.match(r"^\s*((?!end)\b.*)function\s", lines[n], re.IGNORECASE): # Function Start
            functionStart = n
            functionName = lines[n].partition("function")[-1].split("(")[0].split("&")[0].strip()
            
        if re.match(r"^\s*end\sfunction\s", lines[n], re.IGNORECASE): # Function End
            functionEnd = n
            moduleFunctions.append(Element(functionName, functionUses, functionStart, functionEnd))
            functionName = ""
            functionUses = list()
            functionStart = -1
            functionEnd = -1
        
        if re.match(r"^\s*subroutine\s", lines[n], re.IGNORECASE): # Subroutine Start
            subroutineStart = n
            subroutineName = lines[n].partition("subroutine")[-1].split("(")[0].split("&")[0].strip()
            
        if re.match(r"^\s*end\ssubroutine\s", lines[n], re.IGNORECASE): # Subroutine End
            subroutineEnd = n
            moduleSubroutines.append(Element(subroutineName, subroutineUses, subroutineStart, subroutineEnd))
            subroutineName = ""
            subroutineUses = list()
            subroutineStart = -1
            subroutineEnd = -1
        
        if re.match(r"^\s*module\s", lines[n], re.IGNORECASE): # Module Start
            moduleStart = n
            #print lines[n]
            moduleName = lines[n].partition("module")[2].strip()
            
        if re.match(r"^\s*end\smodule\s", lines[n], re.IGNORECASE): # Module End
            moduleEnd = n
            retModules.append(Module(moduleName, moduleSubroutines, moduleFunctions, \
                moduleVariables, moduleDerivedTypes, moduleUses, moduleStart, moduleEnd))
            moduleName = ""
            moduleStart = -1
            moduleEnd = -1
            moduleSubroutines = list()
            moduleFunctions = list()
            moduleVariables = list()
            moduleUses = list()
            moduleDerivedTypes = list()
            
    return retModules

# Parses an API (reads the public vars and subroutines)
def parseApiModule(filename):
    # Take in the entire file at once (small enough text files to not matter)
    lines = [line.strip().lower() for line in open(filename)]
    
    moduleName = filename.split("/")[-1].replace(".F90", "").strip()
    modulePublics = list()
    
    for n in range(0,len(lines)):
        if re.match(r"^\s*public", lines[n], re.IGNORECASE):
            currentLine = getContinuousLine(n, lines)
            currentLine = currentLine.replace("public", "").replace(" ", "").strip()
            for subroutine in currentLine.split(","):
                modulePublics.append(Element(subroutine, None, n, n))
    
    return Module(moduleName, modulePublics, None, None, None, None, 1, 100000)

# Returns the use case on the indicated line
def parseUse(n, lines, parentElement):
    moduleName = ""
    uses = list()
    lines[n] = getContinuousLine(n, lines).replace(" ", "")
    moduleName = lines[n].split(",only:")[0][3:]
    for use in lines[n].split(",only:")[-1].split(","):
        uses.append(use.split("=>")[-1])
    return Use(moduleName, uses, parentElement)

# Returns a list of every use in a module and its children elements    
def getTotalUsesInModule(module):
    retUses = list()
    retUses.extend(module.uses)
    for subroutine in module.subroutines:
        retUses.extend(subroutine.uses)
    for function in module.functions:
        retUses.extend(function.uses)
    return retUses

# Reformats lines split using ampersands to one line
def getContinuousLine(n, lines):
    if re.match(r"^\s*!", lines[n], re.IGNORECASE) or \
        lines[n].startswith("#") or \
        lines[n].strip() == "":              #ignore comments, empty lines, and #
        lines[n]="&"
    lines[n] = lines[n].split("!")[0].strip() # Split off any comments
    if lines[n].endswith("&"):
        lines[n] = lines[n].partition("&")[0]
        lines[n] += getContinuousLine(n+1, lines)
    return lines[n]

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

# Formats all this fancy data into an HTML table
def makeTable(apiModule, samModules, wrfModules, camModules):
    table =  '<!DOCTYPE html><html><script src="sorttable.js"></script>'            # Setup page 
    table += '<style>table,td, th{border:1px solid black;border-collapse: '         # Table CSS
    table += 'collapse;}th {background-color: lightgrey;}</style><head><title>'     # "       "
    table += 'Usage Breakdown: ' + strftime("%Y-%m-%d %H:%M:%S")                    # Page Title
    table += '</title></head><body><table class="sortable">'                        # Setup Table
    table += '<tr><th>API Element</th><th>SAM</th><th>WRF</th><th>CAM</th></tr>'    # Add the Header
    
    # Every element (and use!) in API
    apiElements = list()
    apiElements.extend(apiModule.subroutines)
    
    # add host models to list
    hostModels = list()
    #samModules = cleanupUsesInModules(samModules)
    #wrfModules = cleanupUsesInModules(wrfModules)
    #camModules = cleanupUsesInModules(camModules)
    hostModels.append(samModules)
    hostModels.append(wrfModules)
    hostModels.append(camModules)
    
    # Go through every element in every use in every module in each host model in every element in the api...
    for apiElement in apiElements:
        if not isinstance(apiElement, str):
            table += "<tr><td>" + apiElement.name + "</td>"
        else:
            table += "<tr><td>" + apiElement + "</td>"
       
        for hostModel in hostModels:
            elementsFound = set()
            for module in hostModel: 
                for use in module.uses:
                    if apiElement in use.elements and use.name == apiModule.name:
                        elementsFound.add(use.parentElement)
                for function in module.functions:
                    for use in function.uses:
                        if apiElement in use.elements and use.name == apiModule.name:
                            elementsFound.add(use.parentElement + \
                                "<span style='font-size:"+str(fontSizeReduced)+"pt'> in " + \
                                module.name+"<span style='font-size:"+str(fontSize)+"pt'>")
                for subroutine in module.subroutines:
                    for use in subroutine.uses:
                        if apiElement in use.elements and use.name == apiModule.name:
                            elementsFound.add(use.parentElement + \
                                "<span style='font-size:"+str(fontSizeReduced)+"pt'> in " + \
                                module.name+"<span style='font-size:"+str(fontSize)+"pt'>")
            if (len(elementsFound) == 0):
                table += "<td></td>"
            else:
                table += "<td><p style='font-size:"+str(fontSize)+"pt'>"
                table += ", <br>".join(elementsFound)
                table += "</p></td>"
            
    table += "</table></body></html>"
    return table

# Prints the -h menu
def printHelpMenu():
    print ""
    print "Cross references how one module (like an API) is used by three other groups of modules (Like SAM, WRF, and CAM)."
    print "Example Usage: python thisProgram.py path/to/api.F90 path/to/host1 path/to/host2 path/to/host3 > output.html"
    print ""
    print "This program can also display information about a file."
    print "Example Usage: python thisProgram.py path/to/file.F90"
    print ""

def main():
    if (len(sys.argv) == 2):
        baseModule = parseModulesInFile(sys.argv[1]).pop()
        print baseModule.name
        print "subroutines: "
        for subroutine in baseModule.subroutines:
            print subroutine, subroutine.startLine, subroutine.endLine
            for use in subroutine.uses:
                print "    ", use, "only: ", ", ".join(use.elements)
        print "functions: "
        for function in baseModule.functions:
            print function
            for use in function.uses:
                print "    ", use, "only: ", ", ".join(use.elements)
        print "variables: ", ([str(variable.name) for variable in baseModule.variables])
        print "derivedTypes: ", ([str(derivedType.name) for derivedType in baseModule.derivedTypes])
        print "uses:"
        for use in baseModule.uses:
            print use, "only: ", ", ".join(use.elements)
    elif (len(sys.argv) == 5):
        samModules = list()
        wrfModules = list()
        camModules = list()
        apiModule = parseApiModule(sys.argv[1])
        for file in findFiles(sys.argv[2]):
            samModules.extend(parseModulesInFile(file))
        for file in findFiles(sys.argv[3]):
            wrfModules.extend(parseModulesInFile(file))
        for file in findFiles(sys.argv[4]):
            camModules.extend(parseModulesInFile(file))
        print makeTable(apiModule, samModules, wrfModules, camModules)
    else:
        raise IOError("An error occurred processing the given arguments.")

if __name__ == "__main__":
    try:
        if (sys.argv[1] == "-h"):
            printHelpMenu()
        else:
            main()
    except IOError:
        print "An error occurred processing the given arguments."
