#!/usr/bin/env python
# $Id$
#-------------------------------------------------------------------------------
# This script can be used to convert a SAM "sfc" file to a "*_sfc.in" file
# for use by CLUBB.
#
# The sfc file must be formatted as follows:
# - The input file must contain five columns. Each line of the file is a row.
# - In each row, columns are separated by spaces.
# - The first row contains the names of the columns. They should be:
#   day,  sst(K),  H(W/m2),  LE(W,m2),  TAU(m2/s2)
#
# The output file has columns in the same format with the following names:
#   Time[s], latent_ht[W\m^2], sens_ht[W\m^2], T_sfc[K]
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Main code
#-------------------------------------------------------------------------------
import sys
import common_utils
# Check arguments
if len(sys.argv) != 4:
    print "Usage:"
    print "./sfc_to_sfc.py inputFile outputFile day0"
    print "Where inputFile is the sfc file, outputFile is the *_sfc.in file\
, and day0 is the value for the day field that corresponds to t=0 in CLUBB."
    sys.exit(1)

inputFile = sys.argv[1]
outputFile = sys.argv[2]
day0 = float(sys.argv[3])

# This variable will be used to store the entire file.
inLines = None

with open(inputFile, 'r') as f:
    # Read entire file into memory!
    inLines = f.readlines()

if inLines == None:
    # Some kind of I/O error occured. Bail out.
    print "Something went wrong!"
    sys.exit(2)

outLines = ["Time[s]    latent_ht[W\\m^2]   sens_ht[W\\m^2]    T_sfc[K]\n"]

T_sfc_const = "299.27"

# Loop through all the input lines. Skip the first line.
for i in range(1, len(inLines)):
    inElements = common_utils.parseLine(inLines[i])
    if len(inElements) != 5:
        raise Exception("A row of the file did not have five columns!")
    day = float(inElements[0])
    day = str(int((day - day0) * 86400.0)) + '.'
    outElements = [day, inElements[3], inElements[2], T_sfc_const]
    sizes = [11,19,18,0]
    outLines.append(common_utils.formatOutput(outElements, sizes))

# Write the output lines
with open(outputFile, 'w') as f:
    f.writelines(outLines)
