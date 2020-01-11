# $Id$

# check_for_errors.py
# Author: Kenneth Connor
# Date: July 2011
# Larson-Group UWM

# This script checks CLUBB code for several common mistakes. It checks for
# uninitialized output variables, magic numbers, and magic flags. This script
# references methods in check_uninitialized_output_variables.py and 
# check_magic_numbers.py.
#
# Note that whenever this code refers to 'subroutine', this refers to a 
# subroutine OR a function.


import sys # Handles command line arguments
sys.path.append( './check_scripts/' ) # use modules from this subdirectory
import check_uninitialized_output_variables
import check_magic_numbers
import check_exponents
import glob

#----------------------------------------------------------------------------
def split_into_subroutines_and_functions(lines):

# Splits a list of lines from a file into a list of lists, where each list
# contains the lines of one subroutine or function.
#
# INPUT
#
# lines: a list of strings containing all of the lines from a .F90 file.
#
# OUTPUT
#
# A list of subroutines. Each subroutine is a list of strings containing the
# lines of a subroutine.
#----------------------------------------------------------------------------

  subroutines = []
  current_subroutine = []

  i = 0
  while( i < len(lines)):
    line = lines[i].rstrip()
    line_number = i + 1

    known_item = False

    # Append this line if it is the start of a subroutine or if the start
    # of the subroutine has already been found
    if( (((line.find("subroutine") != -1) or (line.find("function") != -1)) and
         (line.find("end") == -1 )) or (len(current_subroutine) > 0) ):
      # Remove all comments
      if( line.find("!") != -1 ):
        # Mark any lines that are marked as known magic numbers
        if(line.lower().find("known magic number") != -1 or
           line.lower().find("known magic flag") != -1 ):
          known_item = True
        line = line[:line.find("!")].strip()
      # ignore ifdefs and endifs
      if( line.find("#ifdef") != -1 ):
        line = ""
      if( line.find("#endif") != -1 ):
        line = ""
      
      # ignore blank lines
      if( line != "" ):

        # while the last character in the line is '&', the next line is
        # technically part of this one, so include it
        while( line[len(line) - 1] == '&' ):
          line = line.rstrip('&') # remove the '&'

          nextline = lines[i+1].strip()
          # continue over empty lines
          if( nextline == "" ):
            nextline = "&"
          # Remove all comments
          if( nextline.find("!") != -1 ):
            # Mark any lines that are marked as known magic numbers
            if(nextline.lower().find("known magic number") != -1 or
                nextline.lower().find("known magic flag") != -1 ):
              known_item = True

            index = nextline.find("!")
            nextline = nextline[:nextline.find("!")].strip()
            if( index == 0 ):
              nextline += "&"
          if( nextline.find("#ifdef") != -1 ):
            nextline = "&"
          if( nextline.find("#endif") != -1 ):
            nextline = "&"

          line += nextline.strip() # add the next line and remove whitespace
          
          i += 1
     
        # re-insert a comment indicating a known magic item if needed
        if( known_item ):
          line += " ! known magic item"

        # add the line number to the begining of the line and append it
        current_subroutine.append(str(line_number) + ": " + line)

    # If this line is the end of the subroutine, add the subroutine to
    # the list of subroutines and clear current_subroutine
    if( line.find("end subroutine") != -1 or line.find("end function") != -1 ):
      subroutines.append(current_subroutine)
      current_subroutine = []

    i += 1

  return subroutines
# END split_into_subroutines_and_functions



#-------------------BEGIN MAIN CODE-------------------------

# local variables
total_not_set = 0
warnings = 0
magic_numbers = 0
bad_exponents = 0

# if True, warnings will be printed when an intent(out) variable is set by
# a subroutine or function call
show_warnings = False

files_to_check = []

# If there are no arguments other than the show warnings option, 
# check all of the CLUBB code
if( len(sys.argv) == 2 and 
     (sys.argv[1] == "-w" or sys.argv[1] == "--show-warnings") ):
  files_to_check.append("-w")
  files_to_check += glob.glob("../src/*.F90")
  files_to_check += glob.glob("../src/Benchmark_cases/*.F90")
  files_to_check += glob.glob("../src/CLUBB_core/*.F90")
  files_to_check += glob.glob("../src/KK_microphys/*.F90")
  files_to_check += glob.glob("../src/Latin_hypercube/*.F90")
  files_to_check.remove("../src/gfdl_activation.F90")
elif( len(sys.argv) == 1 ):
  files_to_check += glob.glob("../src/*.F90")
  files_to_check += glob.glob("../src/Benchmark_cases/*.F90")
  files_to_check += glob.glob("../src/CLUBB_core/*.F90")
  files_to_check += glob.glob("../src/KK_microphys/*.F90")
  files_to_check += glob.glob("../src/Latin_hypercube/*.F90")
  files_to_check.remove("../src/gfdl_activation.F90")
else:
  files_to_check = sys.argv[1:]

# sys.argv[0] is this file, so skip it
for arg in files_to_check:
  # check for the option to show warnings
  if( arg == "-w" or arg == "--show-warnings" ):
    show_warnings = True
  # ignore any files that are not .F90 files
  elif( arg.find(".F90") != -1 ):
    # open the file for reading
    f = open(arg, "r")
    try:
      # read all lines
      lines = f.readlines()
    finally:
      f.close()

    file_name_printed = False

    # split the lines into subroutines
    subroutines = split_into_subroutines_and_functions(lines) 

    for subroutine in subroutines:
      # check for uninitialized output variables in this subroutine
      outputs = check_uninitialized_output_variables.check_output_variables(subroutine,
        show_warnings)
      # print the current file name if it has not yet been printed
      if( len(outputs) > 0 and not file_name_printed ):
        print "\n\nFile: " + f.name
        file_name_printed = True
      # print out each line of output and increment counters
      for output in outputs:
        if( output.find("WARNING") != -1 ):
          warnings += 1
        elif output.find("Not Set") != -1:
          total_not_set += 1

        print output

      # check for magic numbers and magic flags
      magic_output = check_magic_numbers.check_magic_numbers(subroutine)
      if( len(magic_output) > 0 ):
        # print out the current file name if it has not yet been printed
        if( not file_name_printed ):
          print "\n\nFile: " + f.name
          file_name_printed = True
        # print each line of output and increment the counter
        for line in magic_output:
          print line
          magic_numbers += 1

      # check for bad exponents
      exponent_output = check_exponents.check_exponents(subroutine)
      if( len(exponent_output) > 0 ):
        # print out the current file name if it has not yet been printed
        if( not file_name_printed ):
          print "\n\nFile: " + f.name
          file_name_printed = True
        # print each line of output and increment the counter
        for line in exponent_output:
          print line
          bad_exponents += 1


# print out totals
print "\nTotal Not Set: " + str(total_not_set)
if( show_warnings ):
  print "Total Warnings: " + str(warnings)
print "Total Magic Numbers: " + str(magic_numbers)
print "Total Bad Exponents: " + str(bad_exponents)

