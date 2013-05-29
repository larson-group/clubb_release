# $Id$
#-------------------------------------------------------------------------------
# Splits the given line into a number of columns, separated by any number of
# spaces.
#
# Parameters:
#  - line: the line to split
# Return value:
#  - an array of strings corresponding to the columns in the line
#-------------------------------------------------------------------------------
def parseLine(line):
    columns = []
    # strip the trailing newline from the line, if it exists
    if len(line) > 0 and line[len(line)-1] == '\n':
        line = line[0:len(line)-1]
    previousSpaceIndex = -1
    endofline = False
    while not endofline:
        # Find a space.
        spaceIndex = line.find(' ', previousSpaceIndex+1)
        if spaceIndex == -1:
            # Parse the last column and we are done
            if len(line) > (previousSpaceIndex+1):
                columns.append(line[previousSpaceIndex+1:len(line)])
            endofline = True
        elif spaceIndex-previousSpaceIndex == 1:
            # This space was right after the previous space, so ignore it
            pass
        else:
            # Capture this
            columns.append(line[previousSpaceIndex+1:spaceIndex])
        previousSpaceIndex = spaceIndex
    return columns

#-------------------------------------------------------------------------------
# Formats the given columns into a string so that they take up exactly the
# specified number of places, adding extra spaces to make up the difference.
#
# Parameters:
#  - tokens: an array of entries for the columns
#  - sizes:  an array (same size) that gives the sizes for the columns
#-------------------------------------------------------------------------------
def formatOutput(tokens, sizes):
    outstr = ''
    for i in range(0, len(tokens)):
        outstr += tokens[i]
        spaces = sizes[i] - len(tokens[i])
        if spaces > 0:
            outstr += (' ' * spaces)
    outstr += '\n'
    return outstr
