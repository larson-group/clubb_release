# $Id$

# check_magic_numbers.py
# Author: Kenneth Connor
# Date: July 2011
# Larson-Group UWM

# This file contains methods to detect magic numbers and magic flags in a 
# piece of code. A magic number is any number within a function or subroutine 
# call or used in an equation that is not an integer between -6 and +6. A magic
# flag is a .true. or .false. flag found in a subroutine or function call.
# The check_magic_numbers(subroutine) method should be called to use this script
#
# Note that whenever the code refers to 'subroutine', this refers to a 
# subroutine OR a function.
import string

#----------------------------------------------------------------------------
def is_magic_number_used(line):

# Checks to see if a magic number or magic flag is used in a particular line.
#
# INPUT
#
# line: the line to check
#
# OUTPUT
#
# True if there is a magic number or magic flag used in this line. False otherwise.
#----------------------------------------------------------------------------

  found = False

  equals_index = line.find("=")
  call_index = line.find("call")

  # trimmed contains the trimmed line
  trimmed = ""

  # trim off the appropriate portion of the line
  if( (line.find("if") == -1 and line.find("then") == -1) and 
      line.find("::") == -1 ):
    if( equals_index != -1 ):
      trimmed = line[equals_index+1:].strip()
    elif( call_index != -1 ):
      trimmed = line[call_index+4:].strip()
  
  # if trimmed is still blank, then there was no call or assignment on this line
  if( trimmed != "" and line.lower().find("parameter") == -1 ):
    # these are all of the possible characters that could separate two variables
    separators = ['(', ')', '*', '/', '+', '-', ',', '.and.', '.or.']
    
    vars_found = 0

    # just loop until hitting a break
    while( True ):
      # these are the indexes before and after a variable in the line
      before_index = -1
      after_index = -1
      
      # find the separator that shows up first in the line
      # if this is the first variable to be found, just leave it as -1
      if( vars_found > 0 ):
        for separator in separators:
          index = trimmed.find(separator)
          if( index > -1 and (index < before_index or before_index == -1) ):
            before_index = index 

      # if there was no separator found and there was already a variable found,
      # there are no more variables, so just break out of the while loop
      if( before_index == -1 and vars_found > 0 ):
        break

      # increment before_index so the separator isn't included
      before_index += 1

      # find the first separator after before_index
      for separator in separators:
        index = trimmed.find(separator, before_index)
        
        # ignore symbols in scientific notation and ignore - or + signs if
        # this is the first variable on the line
        while( index > -1 and ((separator == '-' or separator == '+') and 
              (trimmed[index-1].lower() == 'e' or vars_found == 0 )) ):
          index = trimmed.find(separator, index+1)

        if( index > -1 and (index < after_index or after_index == -1) ):
          after_index = index

      # if there is no separator after before_index, just use the end of the line
      if( after_index == -1 ):
        # if there was no separator after the first variable, this is just a
        # normal assignment, so leave it alone
        if( vars_found == 0 ):
          break
        after_index = len(trimmed)

      variable = trimmed[before_index:after_index].strip()
      vars_found += 1

      # check for double precision
      double_index = variable.lower().find("d")
      if( double_index != -1 and (variable[double_index-1] == '.' or
          string.digits.find(variable[double_index-1]) != -1) ):
        variable = variable[:double_index]

      # If the variable is blank, ignore it
      if( variable != "" ):
        # if there is a .true. or .false., this is a magic flag
        if( variable.lower() == ".true." or variable.lower() == ".false." ):
          found = True
          break
        # since a Fortran variable must begin with an alpha character, if the
        # first letter of variable is a digit, this is a magic number.
        # If there is a colon, this is not a variable
        elif( (string.digits.find(variable[0]) != -1 or variable[0] == '-') and
              variable.find(":") == -1 ):
          # If there is an underscore at the end of the number, the underscore
          # is indicating a special precision. If this is the case, trim off the
          # precision identifier
          if( variable.find("_") != -1 ):
            variable = variable[:variable.find("_")]
          try:
            num = float(variable)

            # -999 is our special number indicating an undefined value
            # Since this code parses values by minus signs, negative numbers
            # are read as positive numbers, so check if any instance of 999
            # is actually -999.
            if( num == 999 and trimmed[before_index-1] == '-' ):
              num = -999

            # integer values of -6 to 6 are not magic numbers
            if( (int(num) != num or num < -6 or num > 6) and num != 0.5 and
                 num != -999.0 ):
              found = True
              break
          except ValueError:
            # Uncomment the following line for debug information
            #print variable + ": " + line
            pass
          
      # remove this variable from trimmed
      trimmed = trimmed[after_index:]

  return found
# END is_magic_number_used


#----------------------------------------------------------------------
def get_line_number(line):

#
# Get the line number for a given line.
#
# INPUT
#
# line: a string containing a line of code
#
# OUTPUT
#
# The line number as a string
#----------------------------------------------------------------------
 
  return line[:line.find(":")]

#---------------------------------------------------------------------
def check_magic_numbers(subroutine):
#
# Checks a subroutine for magic numbers and magic flags.
#
# INPUT
#
# subroutine: a list of strings containing the lines of a subroutine
#
# OUTPUT
#
# A list of strings containing all of the magic nubmers and magic flags found
# along with the line number that each appears on.
#------------------------------------------------------------------------
  output = []

  for line in subroutine:
    if( line.find("known magic item") == -1 ):
      if( is_magic_number_used(line) ):
        output.append("    Magic Number or Magic Flag found at line " + 
                       get_line_number(line))

  return output

      
      
