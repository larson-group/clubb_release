#! /usr/bin/python3
# Author: Cavyn VonDeylen
# Date: August 2010
# Larson-Group UWM

from __future__ import print_function
from __future__ import division
from __future__ import unicode_literals
import struct # Handles binary data

#--------------------------------------------------------------------------------------------------
def readGradsData(fileName, numLevels, begTime, endTime, varNum, numVars):
    """
    Reads a GrADS *.dat file and obtains a single-time or time-averaged profile.
    Input: filename: A GrADS *.dat file
           numLevels: Number of z levels in profile
           begTime: Iteration to start averaging at
           endTime: Iteration to end averaging at
           varNum: Which variable to read (see .ctl file)
           numVars: Total number of variables in grads file (see .ctl file)
    """
    
    timeInterval = (endTime-begTime) + 1
        
    # Open in read-binary mode
    dataFile = open(fileName, "rb")

    # Declare array with one slot per z level
    avgField = [0] * numLevels

    # Add data from each time iteration to avgField
    time = begTime
    while True: # Strange loop construct because python doesn't have do-while loops
        byte_position = 4*( (varNum-1)*numLevels+numVars*numLevels*(time-1) )
        dataFile.seek(byte_position)
        
        # Read data in for each z level
        zLevel = 0
        while zLevel < numLevels:
            # Read 4 bytes
            binaryData = dataFile.read(4)
            # Translate binary data to a float.
            avgField[zLevel] = avgField[zLevel] + list(struct.unpack("f", binaryData))[0]
            zLevel += 1
            
        time += 1
        if time >= endTime:
            break

    # Divide by total number of iterations to come up
    # with average value across all iterations for each z level
    zLevel = 0
    while zLevel < numLevels:
        avgField[zLevel] = avgField[zLevel]//timeInterval
        zLevel += 1

    dataFile.close()
    return avgField
    
#--------------------------------------------------------------------------------------------------
# Allows this module to be run as a script
if __name__ == "__main__":

    import sys
    
    # If wrong arguments were given, print a helpful message
    if len(sys.argv) != 7:
        print('Arguments must be: filename z_levels beg_time end_time var_number number_vars')
        sys.exit(0)
    
    print(readNetcdfData( sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), \
        int(sys.argv[5]), int(sys.argv[6]) ))
