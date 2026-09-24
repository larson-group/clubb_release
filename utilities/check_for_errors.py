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
#
# USAGE: python3 check_for_errors.py [-w or --show-warnings] [<filename>.F90 ...]
# With no files given, every .F90 file under ../src (relative to this script)
# is checked.


import os   # for building paths relative to this script
import re   # for recognizing subroutine/function/interface statements
import sys  # Handles command line arguments
import glob # to get all files in a directory

script_dir = os.path.dirname(os.path.abspath(__file__))
src_dir = os.path.join(script_dir, "..", "src")

# use modules from the check_scripts subdirectory
sys.path.append(os.path.join(script_dir, "check_scripts"))
import check_uninitialized_output_variables
import check_magic_numbers
import check_exponents

# The start of a subroutine or function, including any prefixes such as
# pure, elemental, or a return type, e.g. "real( kind = core_rknd ) function"
# or "COMPLEX*16 FUNCTION"
routine_start = re.compile(
  r"^(?:(?:pure|elemental|recursive|impure|module|integer|logical|real|"
  r"complex|character|double\s+precision|type|class)\b"
  r"(?:\s*\([^)]*\)|\s*\*\s*\d+)?\s*)*(?:subroutine|function)\s+\w+",
  re.IGNORECASE)
routine_end = re.compile(r"^end\s*(?:subroutine|function)\b", re.IGNORECASE)
contains_statement = re.compile(r"^contains$", re.IGNORECASE)
interface_start = re.compile(r"^(?:abstract\s+)?interface\b", re.IGNORECASE)
interface_end = re.compile(r"^end\s*interface\b", re.IGNORECASE)

#----------------------------------------------------------------------------
def clean_line(raw_line):

# Removes whitespace, comments, and preprocessor directives from a line.
#
# INPUT
#
# raw_line: a line read from a .F90 file
#
# OUTPUT
#
# A tuple of the cleaned line and True if the comment marked the line as a
# known magic number or known magic flag (False otherwise).
#----------------------------------------------------------------------------

  line = raw_line.strip()

  # ignore preprocessor directives (#ifdef, #endif, #include, etc.)
  if( line.startswith("#") ):
    return "", False

  # find the start of the comment, ignoring any '!' inside a quoted string
  quote_char = None
  for i, char in enumerate(line):
    if( quote_char ):
      if( char == quote_char ):
        quote_char = None
    elif( char == '"' or char == "'" ):
      quote_char = char
    elif( char == '!' ):
      comment = line[i:].lower()
      known_item = ( comment.find("known magic number") != -1 or
                     comment.find("known magic flag") != -1 )
      return line[:i].strip(), known_item

  return line, False
# END clean_line

#----------------------------------------------------------------------------
def split_into_subroutines_and_functions(lines):

# Splits a list of lines from a file into a list of lists, where each list
# contains the lines of one subroutine or function. Comments are removed and
# continued lines are joined. Internal procedures (those after a "contains"
# inside a subroutine) get their own list. Interface blocks are skipped,
# since the procedures declared in them have no bodies to check.
#
# INPUT
#
# lines: a list of strings containing all of the lines from a .F90 file.
#
# OUTPUT
#
# A list of subroutines, in the order they start in the file. Each subroutine
# is a list of strings containing the lines of a subroutine, each prefixed
# with its line number, e.g. "12: x = 1".
#----------------------------------------------------------------------------

  subroutines = []
  open_subroutines = [] # the subroutines currently open, innermost last
  after_contains = False
  in_interface = False

  i = 0
  while( i < len(lines) ):
    line_number = i + 1
    line, known_item = clean_line(lines[i])

    # while the last character in the line is '&', the next line is
    # technically part of this one, so include it
    while( line.endswith("&") and i + 1 < len(lines) ):
      line = line[:-1].rstrip() # remove the '&'
      i += 1

      nextline, next_known_item = clean_line(lines[i])
      known_item = known_item or next_known_item

      # blank, comment-only, and preprocessor lines inside a continued
      # statement don't end it, so keep going
      if( nextline == "" ):
        line += " &"
        continue

      # a continuation line may optionally begin with '&'
      if( nextline.startswith("&") ):
        nextline = nextline[1:].lstrip()

      line += " " + nextline

    i += 1

    # ignore blank lines
    if( line == "" ):
      continue

    # skip over interface blocks
    if( in_interface ):
      if( interface_end.match(line) ):
        in_interface = False
      continue
    if( interface_start.match(line) ):
      in_interface = True
      continue

    # A new subroutine starts either at the top level or, as an internal
    # procedure, after a "contains". Any other subroutine statement inside a
    # subroutine (e.g. an alternate header in an #ifdef) is just a line of it.
    if( routine_start.match(line) and
        (len(open_subroutines) == 0 or after_contains) ):
      open_subroutines.append([])
      subroutines.append(open_subroutines[-1])
      after_contains = False
    elif( len(open_subroutines) == 0 ):
      # not inside a subroutine, so ignore this line
      continue
    elif( contains_statement.match(line) ):
      after_contains = True

    # re-insert a comment indicating a known magic item if needed
    if( known_item ):
      line += " ! known magic item"

    # add the line number to the beginning of the line and append it
    open_subroutines[-1].append(str(line_number) + ": " + line)

    # If this line is the end of the subroutine, close it. An enclosing
    # subroutine can only have more internal procedures after this one.
    if( routine_end.match(line) ):
      open_subroutines.pop()
      after_contains = len(open_subroutines) > 0

  return subroutines
# END split_into_subroutines_and_functions



#-------------------BEGIN MAIN CODE-------------------------

def main():

  # local variables
  total_not_set = 0
  warnings = 0
  magic_numbers = 0
  bad_exponents = 0

  # if True, warnings will be printed when an intent(out) variable is set by
  # a subroutine or function call
  show_warnings = False

  files_to_check = []

  for arg in sys.argv[1:]:
    # check for the option to show warnings
    if( arg == "-w" or arg == "--show-warnings" ):
      show_warnings = True
    else:
      files_to_check.append(arg)

  # If no files were given, check all of the CLUBB code
  if( len(files_to_check) == 0 ):
    files_to_check = sorted(glob.glob(os.path.join(src_dir, "**", "*.F90"),
                                      recursive=True))
    files_to_check = [os.path.relpath(f) for f in files_to_check]

  for filename in files_to_check:
    # ignore any files that are not .F90 files
    if( not filename.endswith(".F90") ):
      continue

    try:
      with open(filename, "r", errors="replace") as f:
        lines = f.readlines()
    except OSError as error:
      print("Could not read " + filename + ": " + error.strerror,
            file=sys.stderr)
      continue

    file_name_printed = False

    # split the lines into subroutines
    subroutines = split_into_subroutines_and_functions(lines)

    for subroutine in subroutines:
      # check for uninitialized output variables in this subroutine
      outputs = check_uninitialized_output_variables.check_output_variables(
        subroutine, show_warnings)
      # print the current file name if it has not yet been printed
      if( len(outputs) > 0 and not file_name_printed ):
        print("\n\nFile: " + filename)
        file_name_printed = True
      # print out each line of output and increment counters
      for output in outputs:
        if( output.find("WARNING") != -1 ):
          warnings += 1
        elif output.find("Not Set") != -1:
          total_not_set += 1

        print(output)

      # check for magic numbers and magic flags
      magic_output = check_magic_numbers.check_magic_numbers(subroutine)
      if( len(magic_output) > 0 ):
        # print out the current file name if it has not yet been printed
        if( not file_name_printed ):
          print("\n\nFile: " + filename)
          file_name_printed = True
        # print each line of output and increment the counter
        for line in magic_output:
          print(line)
          magic_numbers += 1

      # check for bad exponents
      exponent_output = check_exponents.check_exponents(subroutine)
      if( len(exponent_output) > 0 ):
        # print out the current file name if it has not yet been printed
        if( not file_name_printed ):
          print("\n\nFile: " + filename)
          file_name_printed = True
        # print each line of output and increment the counter
        for line in exponent_output:
          print(line)
          bad_exponents += 1

  # print out totals
  print("\nTotal Not Set: " + str(total_not_set))
  if( show_warnings ):
    print("Total Warnings: " + str(warnings))
  print("Total Magic Numbers: " + str(magic_numbers))
  print("Total Bad Exponents: " + str(bad_exponents))

if __name__ == "__main__":
  main()
