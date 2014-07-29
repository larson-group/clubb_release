import re, sys

# Represents a module in FORTRAN and contains sets of elements which represent the various
# subroutines, functions, variables, and derived types which are in that module
class Module:
    name = ""
    subroutines = set() # Element set
    functions = set()
    variables = set()
    derivedTypes = set()
    uses = set() # Use set
    
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
    name = ""
    uses = set()
    startLine = -1
    endLine = -1
    def __init__(self, name, uses, startLine, endLine):
        self.name      = name
        self.uses      = uses
        self.startLine = startLine
        self.endLine   = endLine
    def __eq__(self, other):
        return self.name == other.name
    def __ne__(self, other):
        return not self.__eq__(other)
    def __hash__(self):
        return hash(self.name)
        
class Use:
    moduleName = ""
    elements = set() # Strings
    def __init__(self, moduleName, elements):
        self.moduleName = moduleName
        self.elements = elements
    def __eq__(self, other):
        return self.moduleName == other.moduleName
    def __ne__(self, other):
        return not self.__eq__(other)
    def __hash__(self):
        return hash(self.moduleName)


# Given a filename string, returns an array of populated Modules
def parseModulesInFile(filename):
    
    # Take in the entire file at once (small enough text files to not matter)
    lines = [line.strip() for line in open(filename)]
    
    retModules = set() # The modules found in the current file
    
    # Each data type that is currently being parsed is stored temporarily
    moduleName = ""
    moduleStart = -1
    moduleEnd = -1
    moduleSubroutines = set()
    moduleFunctions = set()
    moduleVariables = set()
    moduleUses = set()
    moduleDerivedTypes = set()
    
    subroutineName = ""
    subroutineUses = set()
    subroutineStart = -1
    subroutineEnd = -1
    
    functionName = ""
    functionUses = set()
    functionStart = -1
    functionEnd = -1
    
    derivedTypeName = ""
    derivedTypeStart = -1
    derivedTypeEnd = -1
    
    lineNumber = 0
    for n in range(0,len(lines)):
        lineNumber += 1
        lines[n] = lines[n].split("!")[0] # Split off any comments
        
        if re.match(r"^\s*use", lines[n], re.IGNORECASE): # Use statement
            if (functionStart != -1 and functionEnd == -1): # Used in a function
                functionUses.add(parseUse(n, lines))
            elif (subroutineStart != -1 and subroutineEnd == -1): # Used in a subroutine
                subroutineUses.add(parseUse(n, lines))
            else: # Used in a module
                moduleUses.add(parseUse(n, lines))
                
        if re.match(r"^\s*((integer)|(character)|(logical)|(real)|(type(?=(\(|\s*\())))(?!\sfunction)", lines[n], re.IGNORECASE): # Variable
            if not lines[n-1].split("!")[0].strip().endswith("&"): # If the line before doesn't end with a &
                for variable in getContinuousLine(n, lines).split("::")[-1].split(","): # merge the &'s and split on ","
                    moduleVariables.add(Element(variable.split("=")[0].strip(), None, lineNumber, lineNumber)) # Take left of "=" and add
        
        if re.match(r"^\s*type(?!(\(|\s*\())", lines[n], re.IGNORECASE): # Derived Type Start
            derivedTypeStart = lineNumber
            derivedTypeName = lines[n].split(" ")[-1].strip()
            
        if re.match(r"^\s*end\stype\s",lines[n], re.IGNORECASE): # Derived Type End
            derivedTypeEnd = lineNumber
            moduleDerivedTypes.add(Element(derivedTypeName, None, derivedTypeStart, derivedTypeEnd))
            derivedTypeName = ""
            derivedTypeStart = -1
            derivedTypeEnd = -1
        
        if re.match(r"^\s*((?!end)\b.*)function\s", lines[n], re.IGNORECASE): # Function Start
            functionStart = lineNumber
            functionName = lines[n].partition("function")[-1].split("(")[0].strip()
            
        if re.match(r"^\s*end\sfunction\s", lines[n], re.IGNORECASE): # Function End
            functionEnd = lineNumber
            moduleFunctions.add(Element(functionName, functionUses, functionStart, functionEnd))
            functionName = ""
            functionUses = set()
            functionStart = -1
            functionEnd = -1
        
        if re.match(r"^\s*subroutine\s", lines[n], re.IGNORECASE): # Subroutine Start
            subroutineStart = lineNumber
            subroutineName = lines[n].partition("subroutine")[-1].split("(")[0].strip()
            
        if re.match(r"^\s*end\ssubroutine\s", lines[n], re.IGNORECASE): # Subroutine End
            subroutineEnd = lineNumber
            moduleSubroutines.add(Element(subroutineName, subroutineUses, subroutineStart, subroutineEnd))
            subroutineName = ""
            subroutineUses = set()
            subroutineStart = -1
            subroutineEnd = -1
        
        if re.match(r"^\s*module\s", lines[n], re.IGNORECASE): # Module Start
            moduleStart = lineNumber
            moduleName = lines[n].partition("module")[-1].strip()
            
        if re.match(r"^\s*end\smodule\s", lines[n], re.IGNORECASE): # Module End
            moduleEnd = lineNumber
            retModules.add(Module(moduleName, moduleSubroutines, moduleFunctions, \
                moduleVariables, moduleDerivedTypes, moduleUses, moduleStart, moduleEnd))
            moduleName = ""
            moduleStart = -1
            moduleEnd = -1
            moduleSubroutines = set()
            moduleFunctions = set()
            moduleVariables = set()
            moduleUses = set()
            moduleDerivedTypes = set()
            
    return retModules

# Returns the use case on the indicated line
def parseUse(n, lines):
    moduleName = ""
    uses = set()
    lines[n] = getContinuousLine(n, lines).replace(" ", "")
    moduleName = lines[n].split(",only:")[0][3:]
    for use in lines[n].split(",only:")[-1].split(","):
        uses.add(use)
    return Use(moduleName, uses)

# Reformats lines split using ampersands to one line
def getContinuousLine(n, lines):
    if re.match(r"^\s*!", lines[n], re.IGNORECASE) or lines[n].startswith("#"): #ignore comments and #
        lines[n]="&"
    lines[n] = lines[n].split("!")[0].strip() # Split off any comments
    if lines[n].endswith("&"):
        lines[n] = lines[n].partition("&")[0]
        lines[n] += getContinuousLine(n+1, lines)
    return lines[n]

if __name__ == "__main__":
    module = parseModulesInFile(sys.argv[1]).pop()
    print module.name#, module.startLine, module.endLine
    print "subroutines"
    for subroutine in module.subroutines:
        print "    ",subroutine.name#, subroutine.startLine, subroutine.endLine
    print "functions"
    for function in module.functions:
        print "    ",function.name#, function.startLine, function.endLine
    print "variables"
    for variable in module.variables:
        print "    ", variable.name, variable.startLine
    print "derived types"
    for derivedType in module.derivedTypes:
        print "    ", derivedType.name#, derivedType.startLine, derivedType.endLine
