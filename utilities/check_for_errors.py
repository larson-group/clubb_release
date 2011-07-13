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
import check_uninitialized_output_variables
import check_magic_numbers

#----------------------------------------------------------------------------
def split_into_subroutines_and_functions(lines):

# Splits a list of lines from a file into a list of lists, where each list
# contains the lines of one subroutine or function
#----------------------------------------------------------------------------
  subroutines = []
  current_subroutine = []

  i = 0
  while( i < len(lines)):
    line = lines[i].rstrip()
    line_number = i + 1

    # Append this line if it is the start of a subroutine or if the start
    # of the subroutine has already been found
    if( (((line.find("subroutine") != -1) or (line.find("function") != -1)) and
         (line.find("end") == -1 )) or (len(current_subroutine) > 0) ):
      # Remove all comments
      if( line.find("!") != -1 ):
        line = line[:line.find("!")].strip()
      if( line.find("#ifdef") != -1 ):
        line = ""
      if( line.find("#endif") != -1 ):
        line = ""
      

      if( line != "" ):

        # while the last character in the lind is '&', the next line is
        # technically part of this one, so include it
        while( line[len(line) - 1] == '&' ):
          line = line.rstrip('&') # remove the '&'

          nextline = lines[i+1].strip()
          # Remove all comments
          if( nextline.find("!") != -1 ):
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



total_not_set = 0
warnings = 0
magic_numbers = 0

# sys.argv[0] is this file, so skip it
for arg in sys.argv[1:]:
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
    outputs = check_uninitialized_output_variables.check_output_variables(subroutine)
    # print the current file name if it has not yet been printed
    if len(outputs) > 0 and not file_name_printed:
      print "\n\nFile: " + f.name
      file_name_printed = True
    # print out each line of output and increment counters
    for output in outputs:
      if output.find("WARNING") != -1:
        warnings += 1
      elif output.find("Not Set") != -1:
        total_not_set += 1

      print output

    # checkf or magic numbers and magic flags
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


# print out totals
print "\nTotal Not Set: " + str(total_not_set)
print "Total Warnings: " + str(warnings)
print "Total Magic Numbers: " + str(magic_numbers)

