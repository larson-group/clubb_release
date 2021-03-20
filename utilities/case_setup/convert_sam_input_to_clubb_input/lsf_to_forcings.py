#!/usr/bin/env python
# $Id$
#-------------------------------------------------------------------------------
# The purpose of this script is to convert a SAM "lsf" file to a "*_forcings.in"
# file for use by CLUBB.
#
# The seven columns of the lsf file are expected to be:
#  z[m], p[mb], tpls[K/s], qls[kg/kg/s], uls, vls, wls[m/s]
#
# This script produces the following output columns (in the *_forcings.in file):
#  Press[Pa], thlm_f[K\s], rtm_f[kg\kg\s], um_ref[m\s], vm_ref[m\s], ...
#  um_f[m\s^2], vm_f[m\s^2], omega[mb\hr], ug[m\s], vg[m\s]
#
# The script produces output in the forcings file for every time listed in the
# lsf file.
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Main code
#-------------------------------------------------------------------------------
import sys
import common_utils
# Check arguments!
if len(sys.argv) != 4:
    print "Usage:"
    print "./sfc_to_sfc.py inputFile outputFile day0"
    print "Where inputFile is the lsf file, outputFile is the *_forcings.in \
file, and day0 is the value for the day field that corresponds to t=0 in CLUBB."
    sys.exit(1)

inputFile = sys.argv[1]
outputFile = sys.argv[2]
day0 = float(sys.argv[3])

# This variable will be used to store the entire file.
inLines = None

with open(inputFile, 'r') as f:
    # Read the entire file into memory!!
    inLines = f.readlines()

# This variable stores all of the lines that will be written to the output file.
outLines = []

# This variable specifies the sizes of the field columns in the output file.
outSizes = [17,16,17,14,14,14,14,15,14,0]

# Write out our units
outLines.append(common_utils.formatOutput(["Press[Pa]", "thlm_f[K\\s]", \
                "rtm_f[kg\\kg\\s]", "um_ref[m\\s]", "vm_ref[m\\s]", \
                "um_f[m\\s^2]", "vm_f[m\\s^2]", "omega[mb\\hr]", "ug[m\\s]", \
                "vg[m\\s]"], outSizes))

# This specifies the sizes of the time declaration rows in the output file.
outTimeSizes = [17,0]

# Loop through all of the input lines. Skip the first line.
for inLine in inLines[1:]:
    inTokens = common_utils.parseLine(inLine)
    if len(inTokens) == 3:
        # This is a time declaration
        day = float(inTokens[0])
        clubb_time = str(int((day - day0) * 86400.0)) + '.'
        # This number is the number of rows for this time value in the input
        # table.
        numRows = inTokens[1]
        outLine = common_utils.formatOutput([clubb_time, numRows], outTimeSizes)
        outLines.append(outLine)
    elif len(inTokens) == 7:
        press_Pa = str(float(inTokens[1])*100.0)
        thlm_f = inTokens[2]
        rtm_f = inTokens[3]
        um_ref = inTokens[4]
        vm_ref = inTokens[5]
        um_f = "-999.9"
        vm_f = "-999.9"
        omega = "-999.9"
        ug = "-999.9"
        vg = "-999.9"
        outLine = common_utils.formatOutput([press_Pa, thlm_f, rtm_f, um_ref, \
                  vm_ref, um_f, vm_f, omega, ug, vg], outSizes)
        outLines.append(outLine)
    else:
        raise Exception("An input line did not have the correct number of \
columns.")

# Write the output lines!!
with open(outputFile, 'w') as f:
    f.writelines(outLines)
