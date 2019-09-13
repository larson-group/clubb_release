# $Id$

# check_division.py
# Author: Kenneth Connor
# Date: July 2011
# Larson-Group UWM

# This script finds and displays lines of code with division. One common
# error/coding standard violation is dividing by a variable without
# including a max() statement to ensure the variable is never very small or
# zero. This script prints out all cases of division other than those that
# include a max() or just a constant number in the denominator. The user can
# then look through these results and ensure that they haven't caused any
# potential errors.


import sys # for command line arguments
import glob # to get all files in a directory

#---------------------------------------------------------------------------
def remove_comments_and_continuations(lines):

# Removes all comments from a file and removes continuation characters,
# putting concatenating the next line onto the one with the continuation.
#
# INPUT
#
# lines: a list of strings containing all of the lines from a .F90 file.
#
# OUTPUT
#
# A list of strings containing all of the lines in a file with comments and
# continuations removed and line numbers added.
#---------------------------------------------------------------------------

  return_lines = []

  i = 0
  while( i < len(lines)):
    line = lines[i].rstrip()
    line_number = i + 1

    # Remove all comments
    if( line.find("!") != -1 ):
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
      return_lines.append(str(line_number) + ": " + line)

    
    i += 1

  return return_lines
# END split_into_subroutines_and_functions

#------------------------------------------------------------------------
def remove_quotes(line):

# Removes anything withing double or single quotes from a line. Division
# never occurs between quotes, so they need to be removed.
#
# INPUT
#
# line: the line to remove quotes from
#
# OUTPUT
#
# A string containing line with all text between quotes, including the
# quotes, removed.
#-----------------------------------------------------------------------

  return_line = line

  # remove text between double quotes
  quote_index = return_line.find('"')
  # loop until there are no more double quotes
  while quote_index != -1:
    # ingore quotes that are escaped
    if( line[quote_index - 1] != "\\" ):
      next_index = return_line.find('"', quote_index+1)
      # if another quote is not found, just go to the end of the line
      if( next_index == -1 ):
        next_index = len(line)
      # cut out the unwanted text
      return_line = (return_line[:quote_index] + return_line[next_index+1:])
    # look up the next double quote
    quote_index = return_line.find('"')

  # remove text between single quotes
  quote_index = return_line.find("'")
  # loop until there are no more single quotes
  while quote_index != -1:
    # ingore quotes that are escaped
    if( line[quote_index - 1] != "\\" ):
      next_index = return_line.find("'", quote_index+1)
      # if another quote is not found, just go to the end of the line
      if( next_index == -1 ):
        next_index = len(line)
      # cut out the unwanted text
      return_line = (return_line[:quote_index] + return_line[next_index+1:])
    # look up the next single quote
    quote_index = return_line.find("'")

  return return_line
# END remove_quotes


#----------------------------------------------------------------------
def line_has_division(line):

# Checks if a line contains division. Ignores any division that has a
# max() or just a constant in the denominator.
#
# INPUT
#
# line: the line to check for division
#
# OUTPUT
#
# True if line contains division. False otherwise.
#----------------------------------------------------------------------

  # remove any text in between quotes
  clean_line = remove_quotes(line)

  has_division = False
  equals_index = clean_line.find("=")
  divide_index = clean_line.find("/", equals_index)

  # Only continue if there is an '=' in this line and continue as long
  # as there is another divide symbol
  while( equals_index != -1 and divide_index != -1 ):
    # ignore any double slashes (//)
    if( clean_line[divide_index:divide_index+2] != '//' ):
      has_division = True

      # The denominator ends at one of the end_chars characters. Find the
      # index of this end.
      end_chars = [',', ')']
      after_index = len(clean_line)
      for end_char in end_chars:
        index = clean_line.find(end_char, divide_index)
        if( index > -1 and index < after_index ):
          after_index = index 

      # cut out the denominator
      denom = clean_line[divide_index+1:after_index].strip()
    
      # if the denomenator is a number, ignore this division
      try:
        float(denom)
        has_division = False
      except:
        pass

      # if the denominator contains a max(), ignore this division
      if( denom.find("max") == 0 ):
        has_division = False

      # break out if division was found
      if( has_division == True ):
        break
    
    # find the next divide character
    divide_index = clean_line.find("/", divide_index + 2)

  return has_division


#------------------------BEGIN MAIN CODE---------------------------------

total_divisions = 0

files_to_check = []

# If there are no arguments, check all of the CLUBB code
if( len(sys.argv) == 1 ):
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
  # ignore any files that are not .F90 files
  if( arg.find(".F90") != -1 ):
    # open the file for reading
    f = open(arg, "r")
    try:
      # read all lines
      lines = f.readlines()
    finally:
      f.close()

  # clean all of the lines by removing comments and continuations
  clean_lines = remove_comments_and_continuations(lines)

  # loop through all lines and check for division
  for line in clean_lines:
    if( line_has_division(line) ):
      total_divisions += 1
      # print out each line that has division
      print arg + ", " + line

# Print the total number of divisions found. 
print "Total Divisions Found: " + str(total_divisions)
  

