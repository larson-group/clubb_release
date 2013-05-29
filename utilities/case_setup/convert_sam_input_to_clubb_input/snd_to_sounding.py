#!/usr/bin/env python
# $Id$
#-------------------------------------------------------------------------------
# This script shall be used to convert SAM's snd file to a *_sounding.in file
# for input to CLUBB. This script only takes data from one point in time.
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Main code
#-------------------------------------------------------------------------------
import common_utils
import sys
# Check arguments!
if len(sys.argv) != 3:
    print "Usage:"
    print "./snd_to_sounding.py inputFile outputFile"
    print "Where inputFile is the snd file, and outputFile is the *_sounding.in\
 file to write to"
    sys.exit(1)

inputFile = sys.argv[1]
outputFile = sys.argv[2]

# This variable will be used to store the entire input file.
inLines = None

# Will throw an exception if inputFile could not be opened/read
with open(inputFile, 'r') as f:
    # Read entire file into memory
    inLines = f.readlines()

outLines = []
sizes = [16, 15, 15, 13, 13, 14, 13, 0]
outLines.append(common_utils.formatOutput(["Press[Pa]", "thm[K]", "rt[kg\kg]", \
      "u[m\s]", "v[m\s]", "omega[Pa\s]", "ug[m\s]", "vg[m\s]" ], sizes))

# The second input line should specify the time, and the time row has three
# columns.
if len(common_utils.parseLine(inLines[1])) != 3:
    print "Error: second line of file does not give time information"
    sys.exit(2)

# Now loop through the rest of the lines
for i in range(2, len(inLines)):
    inTokens = common_utils.parseLine(inLines[i])
    if len(inTokens) == 3:
        # We are done here, as this is another row of time information
        break
    elif len(inTokens) != 6:
        print "Error: a file row did not have 6 columns."
        sys.exit(2)
    else:
        press_Pa = str(float(inTokens[1])*100.0)
        thm_K = inTokens[2]
        rt_kg_kg = str(float(inTokens[3]) / 1000.0)
        u_m_s = inTokens[4]
        v_m_s = inTokens[5]
        omega_Pa_s = "-999.9"
        ug_m_s = u_m_s
        vg_m_s = v_m_s
        outTokens = [press_Pa, thm_K, rt_kg_kg, u_m_s, v_m_s, omega_Pa_s, \
                     ug_m_s, vg_m_s]
        outLines.append(common_utils.formatOutput(outTokens, sizes))

# Write output lines
with open(outputFile, 'w') as f:
    f.writelines(outLines)
