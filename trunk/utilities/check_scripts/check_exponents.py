# $Id$

# check_exponents.py
# Author: Kenneth Connor
# Date: August 2011
# Larson-Group UWM

# This script checks all exponents for any floating point exponents that can
# be integer exponents instead, as well as an exponent of 0.5, which could be
# computed with the sqrt function. Fortran computes integer exponents through
# simple multiplication, whereas a floating point exponent is computed using 
# logarithms, which is slower and less accurate. By outputting the lcoations
# of these 'bad' exponents, this script helps the user eliminate them. 

#----------------------------------------------------------------------------
def line_has_bad_exponent(line):
#
# INPUT
#
# line: The line to check for bad exponents.
#
# OUTPUT
#
# True if this line has at least one bad exponent. False otherwise.
#----------------------------------------------------------------------------

  bad_exponent = False

  exp_index = line.find("**")

  # Loop until all exponents are checked or a bad exponent is found
  while( exp_index != -1 and (not bad_exponent)):
    paren_index = line.find("(", exp_index + 2)
    after_index = -1
    exponent = ""
    
    # If an open parenthese is found within 4 characters of the exponent,
    # the exponent is the text between two parenthesis
    if( paren_index != -1 and paren_index < exp_index + 4):
      num_parens = 1
      i = paren_index + 1
      # loop through to find the close parenthese that matches this open parenthese
      while( num_parens > 0 and i < len(line) ):
        if( line[i] == '(' ):
          num_parens += 1
        elif( line[i] == ')' ):
          num_parens -= 1
        i += 1
      # subtract one from the index to get to the right spot
      i -= 1

      exponent = line[paren_index + 1:i]

    # If there was no parenthese found, the exponent is just the text until a
    # separator is found
    else:
      # a list of characters that could mark the end of an exponent
      separators = [')', '*', '/', '+', '-', ',']

      # loop through all separators
      for separator in separators:
        index = line.find(separator, exp_index + 2)
      
        # ignore symbols in scientific notation and ignore - or + signs if
        # this is the first variable on the line
        while( index > -1 and ((separator == '-' or separator == '+') and 
              (line[index-1].lower() == 'e' or index == 0 )) ):
          index = line.find(separator, index+1)

        if( index > -1 and (index < after_index or after_index == -1) ):
          after_index = index

      # If no end parenthese was found, just set after_index to the end of the line
      if( after_index == -1 ):
        after_index = len(line)

      exponent = line[exp_index + 2:after_index].strip()

    # Pull off any parenthesis at the beginning of the exponent
    if( len(exponent) > 0 and exponent[0] == '(' ):
      exponent = exponent[1:]

    # Remove any double precision markers
    if( exponent.find("d0") != -1 ):
      exponent = exponent[:exponent.find("d0")]
      # Make sure there's still a decimal so the script knows this is a double
      if( exponent.find(".") == -1 ):
        exponent += ".0"

    # Remove anything after an underscore
    if( exponent.find("_") != -1 ):
      exponent = exponent[:exponent.find("_")]

    try:
      # If this conversion fails, the exponent is not a single number
      float_value = float(exponent)
      int_value = int(float_value)
      # If float_value is the same as int_value and the exponent is written
      # as an integer, or if the exponent is 0.5, this is a bad exponent
      if( (float_value == int_value and str(int_value) != exponent) or
           abs(float_value) == 0.5 ):
        bad_exponent = True
    except:
      pass

    # Find the next exponent and continue the loop
    exp_index = line.find("**", exp_index + 2)

  return bad_exponent
# END line_has_bad_exponent

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

# END get_line_number

#---------------------------------------------------------------------
def check_exponents(subroutine):
#
# Checks a subroutine for bad exponents.
#
# INPUT
#
# subroutine: a list of strings containing the lines of a subroutine
#
# OUTPUT
#
# A list of strings containing all of bad exponents found along with the
# line number eacg was found on.
#------------------------------------------------------------------------

  output = []

  for line in subroutine:
    if( line_has_bad_exponent(line) ):
      output.append("    Bad exponent found at line " + 
                     get_line_number(line))

  return output

# END check_exponents

