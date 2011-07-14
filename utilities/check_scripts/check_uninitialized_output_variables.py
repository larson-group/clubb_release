# $Id$

# check_uninitialized_output_variables.py
# Author: Kenneth Connor
# Date: June 2011
# Larson-Group UWM

# This file contains methods to check for uninitialized output variables.
# An output variable is uninitialized if it is never initialized within the
# body of a subroutine, or if it is initialized in one part of an if or 
# select statement, but not all of the parts.
#
# The method check_output_variables(subroutine) should be called to use this
# script. 
#
# Note that whenever the code refers to 'subroutine', this refers to a 
# subroutine OR a function.


import string # Used for string comparison


#---------------------------------------------------------------------------
def find_output_vars(subroutine):

# Finds all of the intent(out) variables in a subroutine.
#---------------------------------------------------------------------------

  outputs = []

  i = 0 # an index variable

  # loop through all lines in the subroutine
  for line in subroutine:
    if( line.find("intent(out)") != -1 ):
      # strip the whitespace off of the end
      stripped_line = line.rstrip()
      
#      j = i # temporary index
      # while the last character in the lind is '&', the next line is
      # technically part of this one, so include it
#     while( stripped_line[len(stripped_line) - 1] == '&' ):
#       stripped_line = stripped_line.rstrip('&') # remove the '&'
#       stripped_line += subroutine[j + 1].rstrip() # add the next line and remove whitespace
#       j += 1

      index = stripped_line.find("::") + 2
      outvars = stripped_line[index:] # leave only the variable names

      # Remove anything between parentheses. This isn't part of the variable
      # name.
      while( outvars.find("(") != -1 ):
        index1 = outvars.find("(")
        index2 = outvars.find(")") + 1
        outvars = outvars[:index1] + outvars[index2:]

      # replace all commas with a space and then split by whitespace to
      # separate all variables into a list
      varlist = outvars.replace(',', ' ').split() 
      
      for output in varlist:
        outputs.append(output)

    # increment the index
    i += 1
      
  return outputs
# END find_output_vars

#--------------------------------------------------------------------------
def is_var_set_in_line(var, line):

# Checks if a variable is set in this line. This means that there is an
# equals sign to the right of the variable.
#
# INPUT
#
# var: the varialbe to check
# line: the line to check
#
# OUTPUT
#
# True if var is set in line. False otherwise.
#---------------------------------------------------------------------------

  var_set = False

  # find some indexes
  var_index = line.find(var)
  equals_index = line.find("=", var_index)

  # if var and an equal sign are both found and var is the first thing on the line
  if( var_index != -1 and equals_index != -1 and 
     (line[var_index - 1] == ' ' or line[var_index - 1] == '\t') ):
    paren_index = line.find("(", var_index, equals_index)
    modulo_index = line.find("%", var_index, equals_index)
    # look for parenthesis in case it's a dimension variable
    if( paren_index != -1  and
        line.find(")", paren_index, equals_index) != -1 ):
      var_set = True
    elif( paren_index == -1 and (line[var_index + len(var)] == ' ' or
                                 line[var_index + len(var)] == '=' ) ):
      var_set = True
    # look for a modulo in case it's a type variable
    elif( modulo_index == var_index + len(var) ):
      var_set = True


  return var_set
#END is_var_set_in_line

#--------------------------------------------------------------------------
def is_var_set_by_subroutine(var, line):

# Checks if a variable is set by a subroutine or function call in a line.
#
# INPUT
#
# var: the name of the variable to check
# line: the line to check in
#
# OUTPUT
#
# True if the variable is set by a subroutine of function. False otherwise.
#--------------------------------------------------------------------------

  set = False

  var_index = line.find(var)
  open_paren_index = line.find("(")
  # check each pair of parenthesis in the line
  while( line.find("(", open_paren_index) != -1 ):
    num_parens = 1
    i = open_paren_index + 1
    # loop through to find the close parenthese that matches this open parenthese
    while( num_parens > 0 and i < len(line) ):
      if( line[i] == '(' ):
        num_parens += 1
      elif( line[i] == ')' ):
        num_parens -= 1
      i += 1

    close_paren_index = i

    # if var is in this line and is in between two parenthesis
    if( var_index != -1 and open_paren_index != -1 and close_paren_index != -1 and
        var_index > open_paren_index and var_index < close_paren_index ):
      # trim everything before the parenthesis
      trim = line[:open_paren_index].rstrip()
      # this is a subroutine or function if it is called or if there is
      # a subroutine name before the parenthesis
      if( line.find("call", open_paren_index) != -1 or 
          string.ascii_letters.find(trim[len(trim)-1]) != -1 or
          string.digits.find(trim[len(trim)-1]) ):
        # check to make sure this is actually the variable we're looking for
        # and not part of another variable's name
        if( (line[var_index-1] == ' ' or line[var_index-1] == ',') and
            (line[var_index + len(var)] == ' ' or line[var_index + len(var)] == ',' or
             line[var_index + len(var)] == '(' ) ):
          set = True

          break # stop looking

    # if var wasn't set in this paren block, look again in the next one
    open_paren_index = line.find( "(", open_paren_index + 1 )

  return set
# END is_var_set_by_subroutine
    
#---------------------------------------------------------------------------
def is_var_set(var, subroutine):

# Checks if a variable has been set in a particular subroutine.
#
# INPUT
# 
# var: the name of the variable to check
# subrotuine: a list of the lines that make up the subroutine
#
# OUTPUT
#
# 0 - var is not set
# 1 - var is set
# 2 - var is set in subroutine
#---------------------------------------------------------------------------

  var_set = 0
  
  skip_if = False
  skip_case = False

  i = 1 # index variable
  while( i < len(subroutine) ):
    line = subroutine[i]

    if( line.lower().find("if") != -1 and line.lower().find("then") != -1 and
        line.lower().find("else") == -1 and line.lower().find("end") == -1 ):
      end_if_index = i+1
      num_ifs = 1
      # find the end index for the if block
      while( num_ifs > 0 ):
        temp_line = subroutine[end_if_index].lower()
        if( temp_line.find("if") != -1 and temp_line.find("then") != -1 and
            temp_line.find("else") == -1 and temp_line.find("end") == -1 ):
          num_ifs += 1
        if( temp_line.find("end if") != -1 or temp_line.find("endif") != -1 ):
          num_ifs -= 1

        end_if_index += 1
      
      check = check_var_in_if( var, subroutine[i:end_if_index] )
      if( check != 0 ):
        var_set = check
      #else:
      #  var_set = 0 # not set

      i = end_if_index - 1

    elif( line.lower().find("selectcase") != -1 or line.lower().find("select case") != -1):
      end_select_index = i+1
      num_cases = 1
      # find the end index for the case block
      while( num_cases > 0 ):
        temp_line = subroutine[end_select_index].lower()
        if( temp_line.find("selectcase") != -1 or temp_line.find("select case") != -1 ):
          num_cases += 1
        if( temp_line.find("end select") != -1 or temp_line.find("endselect") != -1 ):
          num_cases -= 1

        end_select_index += 1

      check = check_var_in_case( var, subroutine[i:end_select_index] )
      if( check != 0 ):
        var_set = check
      #else:
      #  var_set = 0 # not set

      i = end_select_index - 1


    # if the variable is set in this line
    elif( is_var_set_in_line(var, line) ):
      var_set = 1

    # If the variable is called in a subroutine, set the subroutine_warning flag
    elif( is_var_set_by_subroutine(var, line) ):
        var_set = 2 # set in subroutine

    # If var has been deemed as set, stop looking
    if( var_set == 1 ):
      break

    i += 1

  return var_set
# END is_var_set
      

#---------------------------------------------------------------------------
def check_var_in_if(var, if_statement):

# Check if var is set in all sections of an if-block
#
# INPUT
#
# var: the name of the variable to check
# if_statement: a list containing the lines of the if statement
#
# OUTPUT
#
# 0 - var is not set
# 1 - var is set
# 2 - var is set in subroutine
#--------------------------------------------------------------------------

  found = 0
  warning = False
  else_found = False

  i = 1
  while( i < len(if_statement) ):
    line = if_statement[i]

    if( line.lower().find("else") != -1 and line.lower().find("if") == -1 ):
      else_found = True

    if( line.lower().find("if") != -1 and line.lower().find("then") != -1 and
        line.lower().find("else") == -1 and line.lower().find("end") == -1 ): 
      end_if_index = i+1
      num_ifs = 1
      # find the end index for the if block
      while( num_ifs > 0 ):
        temp_line = if_statement[end_if_index].lower()
        if( temp_line.find("if") != -1 and temp_line.find("then") != -1 and
            temp_line.find("else") == -1 and temp_line.find("end") == -1 ):
          num_ifs += 1
        if( temp_line.find("end if") != -1 or temp_line.find("endif") != -1 ):
          num_ifs -= 1

        end_if_index += 1

      if( found == 0 ):
        check = check_var_in_if( var, if_statement[i:end_if_index] )
        if( check == 1 ):
          found = 1 # set
        elif( check == 2 ):
          found = 2 # set in subroutine
     
      i = end_if_index - 1

    elif( line.lower().find("selectcase") != -1 or line.lower().find("select case") != -1 ):
      end_select_index = i+1
      num_cases = 1
      # find the end index for the case block
      while( num_cases > 0 ):
        temp_line = if_statement[end_select_index].lower()
        if( temp_line.find("selectcase") != -1 or temp_line.find("select case") != -1 ):
          num_cases += 1
        if( temp_line.find("end select") != -1 or temp_line.find("endselect") != -1 ):
          num_cases -= 1
        end_select_index += 1

      if( found == 0 ):
        check = check_var_in_case( var, if_statement[i:end_select_index] )
        if( check == 1 ):
          found = 1 # set
        elif( check == 2 ):
          found = 2 # set in subroutine

      i = end_select_index - 1


    # if var is set in this section, set found to 1
    elif( is_var_set_in_line(var, line) ):
      found = 1 # set
    # check if var is set in a subroutine
    elif( is_var_set_by_subroutine(var, line) ):
        found = 2 # set in subroutine
        warning = True
    # if the program stops execution in an if block, ignore this block
    elif( line.lower().find("stop") != -1 ):
      templine = line.lower()[line.find(":")+1:].strip()
      if( templine[0:4] == "stop" ):
        found = 1

    # if we have reached the beginning of a new section
    elif( (line.lower().find("elseif") != -1) or
          (line.lower().find("else if") != -1) or
          (line.lower().find("else") != -1) ):
      # if var was not set in the previous section, break and return false
      if( found == 0 ):
        break
      # if var was set, reset found to False and continue with the next section
      else:
        found = 0
    i += 1

  if( not else_found ):
    found = 0

  if( warning and found == 1 ):
    found = 2


  return found
# END check_var_in_if

#-------------------------------------------------------------------------
def check_var_in_case(var, case_statement):

# Check if var is set in all cases of a select-case block
#
# INPUT
#
# var: the name of the variable to check
# case_statement: a list containing the lines of the case statement
#
# OUTPUT
#
# 0 - var is not set
# 1 - var is set
# 2 - var is set in subroutine
#-------------------------------------------------------------------------

  found = 0
  warning = False
  # loop through all but the first line
  i = 2
  while( i < len(case_statement) ):
    line = case_statement[i]

    if( line.lower().find("if") != -1 and line.lower().find("then") != -1 and
        line.lower().find("else") == -1 and line.lower().find("end") == -1 ):
      end_if_index = i+1
      num_ifs = 1
      # find the end index for the if block
      while( num_ifs > 0 ):
        temp_line = case_statement[end_if_index].lower()
        if( temp_line.find("if") != -1 and temp_line.find("then") != -1 and
            temp_line.find("else") == -1 and temp_line.find("end") == -1 ):
          num_ifs += 1
        if( temp_line.find("end if") != -1 or temp_line.find("endif") != -1 ):
          num_ifs -= 1

        end_if_index += 1

      if( found == 0 ):
        check = check_var_in_if( var, case_statement[i:end_if_index] )
        if( check == 1 ):
          found = 1 # set
        elif( check == 2 ):
          found = 2 # set in subroutine
        #else:
        #  found = 0 # not set
     
      i = end_if_index - 1

    elif( line.lower().find("selectcase") != -1 or line.lower().find("select case") != -1 ):
      end_select_index = i+1
      num_cases = 1
      # find the end index for the case block
      while( num_cases > 0 ):
        temp_line = case_statement[end_select_index].lower()
        if( temp_line.find("selectcase") != -1 or temp_line.find("select case") != -1 ):
          num_cases += 1
        if( temp_line.find("end select") != -1 or temp_line.find("endselect") != -1 ):
          num_cases -= 1
        end_select_index += 1

      if( found == 0 ):
        check = check_var_in_case( var, case_statement[i:end_select_index] )
        if( check == 1 ):
          found = 1 # set
        elif( check == 2 ):
          found = 2 # set in subroutine
        #else:
        #  found = 0 # not set
      i = end_select_index - 1


    # if var was set in this case, set found to 1
    elif( is_var_set_in_line(var, line) ):
      found = 1 # set
    # check if var is set in a subroutine
    elif( is_var_set_by_subroutine(var, line) ):
        found = 2 # set in a subroutine
        warning = True
    # if the program stops execution in a case, ignore this case
    elif( line.lower().find("stop") != -1 ):
      templine = line.lower()[line.find(":")+1:].strip()
      if( templine[0:4] == "stop" ):
        found = 1

    # if this line is the start of a new case
    elif( line.lower().find("case") != -1 ):
      # if var was not set in the previous case, return false
      if( found == 0 ):
        break
      # if var was set in the previous case, set found to 0 and keep looking
      else:
        found = 0
    i += 1
    
  if( warning and found == 1 ):
    found = 2
  return found
# END check_var_in_case
      

#---------------------------------------------------------------------
def get_subroutine_name(subroutine):

# Finds the name of a subroutine.
#
# INPUT
#
# subroutine: a list containing the lines of a subroutine
#
# OUTPUT
#
# the name of the subroutine as a string
#---------------------------------------------------------------------

  first_line = subroutine[0]
  name = first_line[ first_line.find("subroutine") + 11 : first_line.find("(") ]
  return name.strip()
# END get_subroutine_name


#----------------------------------------------------------------------
def get_subroutine_line(subroutine):

# Finds the line number where the subroutine begins
#
# INPUT
#
# subroutine: a list containing the lines of a subroutine
#
# OUTPUT
#
# the line number where the subroutine begins
#---------------------------------------------------------------------

  first_line = subroutine[0]
  line = int( first_line[ 0 : first_line.find(":") ] )
  return line
# END get_subroutine_line


#----------------------------------------------------------------------
def check_output_variables(subroutine):

  subroutine_name_printed = False
  output_text = []

  outputs = find_output_vars(subroutine)

  for output in outputs:
    var_set = is_var_set(output, subroutine)
    if( var_set == 0 or var_set == 2):

      if( not subroutine_name_printed ):
        line_number = get_subroutine_line(subroutine)
        output_text.append( "\n  In subroutine " + get_subroutine_name(subroutine) + 
            " (line " + str(line_number) + ")" )
        subroutine_name_printed = True

      if( var_set == 0 ):
        output_text.append( "    " + output + ": Not Set" )
      elif( var_set == 2 ):
        # If the variable was called in a subroutine and not set anywhere else
        # , print a warning
        output_text.append( "    WARNING: " + output + 
                            " is initialized in a subroutine or function." )
        var_set = True # set var_set to get rid of extra outputa

    #else:
      #print output + ": Set"

  return output_text
