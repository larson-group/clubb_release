#!/usr/bin/python3
###############################################################################################
# Script print_tunable_parameters_table.py                                                    #
# Author: Steffen Domke (Nov 2023)                                                            #
# Description:                                                                                #
# Since the default version of tunable_parameter.in files contain comments and                #
# the ones produced by the tuner do not, it is difficult to diff those files easily.          #
# This script tries to fix that by running some regex patterns over the files that are        #
# passed in. The output of this script is either a list of all parameters listed in those     #
# files and their values, or (not implemented yet) a list of those parameters                 #
# that are different between the files.                                                       #
#                                                                                             #
# The script runs through the lines of the files and matches them to a pattern indicating a   #
# valid parameter definition. Those parameter names and their values are saved in a dict      #
# which is eventually used to create the output.                                              #
# The output is formatted as a markdown table that can be directly copied and pasted into     #
# a github comment and will be formatted as a table there.                                    #
# I have not found  way to nicely name the files in the table so for now they are named as    #
# "file <number>" with a reference to the filename above the table.                           #
#                                                                                             #
# TODO:                                                                                       #
#       - Implement diff feature                                                              #
###############################################################################################

import re
import argparse
# import sys
from os.path import exists

# Define output formats for both int and float values
# TODO: Is this distinction necessary? We could just use float formatting for ints as well
# 1. Define line template for table output depending on the type of the parameter value
table_template = {}
table_template[float] = '{comment}{pval:10.5f}'
table_template[int] = '{comment}{pval:4d}      '
# 2. Define line template for bare file output depending on the type of the parameter value
file_template = {}
file_template[float] = '{comment}{pname:>30s} = {pval:10.5f}\n'
file_template[int] = '{comment}{pname:>30s} = {pval:4d}\n'



parser = argparse.ArgumentParser(description="Parse tunable_parameters for parameter-value pairs")

parser.add_argument('infiles', nargs='+', help='tunable_parameters.in file you want to parse')
parser.add_argument('-o', '--outfile', action='store', help='output file to write the parsed information into')
parser.add_argument('-f', '--format', action='store_true', help='For each input file, create a new tunable_parameters.in file that strips away superfluous content and forces a unified formatting.'+
                    'The new file will be in the same folder as the corresponding input file and be marked by the postfix "-std".')
parser.add_argument('-b', '--bare', action='store_true', help='Reduced output only listing the differences between the given files.')
parser.add_argument('-c', '--commented', action='store', choices=['col', 'sym', 'skip'], default='skip', help='Include commented out paramters in output. If the choice is "skip" (default) the commented out parameters will not be listed. Otherwise, they will be included and either noted in an additional column "commented" or as a symbol "!" in the value column.')

args = parser.parse_args()

# print(args)

# Eliminate duplicate file names
# if len(args.infiles) > len(set(args.infiles)):
#     print("Warning! Duplicates in file names. Each file will only be parsed once.")
#     infiles = [args.infiles[0]]
#     for f in args.infiles[1:]:
#         if f not in infiles:
#             infiles.append(f)
#         else:
#             print("'{}' appeared in the file list more than once. Duplicates were eliminated.".format(f))

lines = []

# Read content of each input file
for inf in args.infiles:
    # print("Read lines from file ", inf)
    with open(inf) as file:
        lines.append(file.readlines())

# For each file, parse params
# Create key-val pair where val is empty list
# Then fill up in each iteration
# Careful about missing params!
params = {}
maxl = 0
# Iterate through input files
for i,line in enumerate(lines):
    # print("File ", i)
    # Iterate through lines in <i>th file
    for l in line:
        # print("Line: ", l)
        l = l.strip()
        # Ignore comment lines
        if l.startswith('!'):
            # print("Commented")
            commented = True
            l = l[1:]
        else:
            # print("Not commented")
            commented = False
        # Find parameter name and corresponding value
        res = re.search('\s*(\w+)\s*=\s*(\d*(\.\d*)?)', l)
        if res: # If match found
            # print("Parameter entry found: ", res)
            # Assign key and val
            key, val = res.groups()[0:2]
            # print("(key, val) = ", (key, val))
            if '.' in val:
                val = float(val)
            else:
                val = int(val)
            # print('(key, val) = ', (key, val))
            # Adjust max length of parameter name (for later formatting)
            if maxl < len(key):
                maxl = len(key)
            # If new key, add entry to dict
            if key not in params:
                # print("New entry for param ", key)
                params[key] = []
            # Fill up list with dummies for missing entries if necessary
            for k in range(i-len(params[key])):
                # print("Adding dummy")
                params[key].append(None)    # TODO: Do we want None or something more specific to avoid having to check for None everywhere?

            # Check if params[key] has an ith entry (only if key already appeared in current file)
            if len(params[key])<=i:
                # This is the first time this parameter has appeared in the current file
                # -> Just append
                # print("No {}th entry yet -> append".format(i))
                params[key].append((val,commented))
            else: # There was a previous entry for <key> in the current file
                # print("Previous entry for param in this file -> What to do?")
                if not commented: # New entry is not commented
                    # print("New entry is uncommented")
                    if not params[key][-1][1]: # Previous entry is not commented
                        print("Error! Parameter {} appears uncommented twice in file {}. Only the first entry will be used.".format(key, args.infiles[i]))
                        continue
                    else: # Old entry is commented -> Use new one
                        # Replace commented entry with uncommented entry
                        # print("Old commented, new umcommented -> Use new entry")
                        params[key][-1] = (val, False)
                else:
                    print("Warning! Parameter {} was already listed in the file. The later entry will be ignored.".format(key))

# Fill up lists if necessary
for p in params:
    for k in range(len(args.infiles)-len(params[p])):
        params[p].append(None)

# print(params)

# Generate output
# Create header for table first
header = '|params|'
separator = '\n|---|'
for i, inf in enumerate(args.infiles):
    print('file {}: {}'.format(i+1,inf))
    header += 'file {}|'.format(i+1)
    separator += '---|'

output = header + separator

# Write lines to output for each parameter
for p in params:
    # Skip parameter if its value is equal between all files
    # The check is needed because absent parameters are represented by a single None
    vals = set([item[0] if item else item for item in params[p]])
    if args.bare and len(vals)==1:
        continue
    # Convert params[p] from list of tuples (<name>, <commented>) to string representation of value with '!' in front if parameter was commented in file
    # pname = params[p][0]
    # comment = '!' if params[p][1] else ' '
    stringlist = [table_template[type(param[0])].format(pval=param[0], comment='!' if param[1] else ' ') if param else ' '*11 for param in params[p]]
    # Generate and append string for line showing info parameter <p>
    output += '\n|{key:>{wid}s}|'.format(key=p,wid=maxl) + '|'.join(stringlist) + '|'

if args.outfile:
    if exists(args.outfile):
        print('File {} already exists and will not be overwritten.'.format(args.outfile))
    else:
        with open(args.outfile,'w') as f:
            f.write(output)
        print('Output written to file {}'.format(args.outfile))

if args.format:
    for i,outf in enumerate(args.infiles):
        with open(outf+"-std","w") as file:
            file.write('! Original tunable_parameters.in file at '+outf+'\n')
            file.write('&clubb_params_nl\n')
            for p in params:
                if params[p][i]:
                    # Is the parameter supposed to be commented out?
                    # If yes, prepend a '!'
                    comment = '!' if params[p][i][1] else ' '
                    # Write line to file for parameter <p>
                    file.write(file_template[type(params[p][i][0])].format(pname=p, pval=params[p][i][0], comment=comment))
            file.write('/')

print(output)