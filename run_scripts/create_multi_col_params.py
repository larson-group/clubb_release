#!/bin/python
import sys
import re
import numpy as np
import argparse

# Parse a parameter file
# - this will match anything with the form var = n1, n2, n3, so ensure
#   that the file contains only parameters with 1 value definition
def parse_clubb_params(file_name):

    params_dict = {}
    
    # Regular expression to match lines like 'var = n1, n2, n3, ...'
    pattern = re.compile(r'(\w+)\s*=\s*(.*)')

    with open(file_name, 'r') as file:
        for line in file:

            # Remove anything after the '!' symbol (trailing comments)
            line = line.split('!')[0].strip()
            
            # If the line is empty after removing the comment, skip it
            if not line: continue

            # Match the pattern and extract variable name and values
            match = pattern.match(line)
            if match:
                var_name = match.group(1)
                number   = match.group(2)
                params_dict[var_name] = float(number)

    return params_dict

# Copy the parameters stored in parsed_params ngrdcol times
# For exampe:
#   C1 = 1.0
# becomes
#   C1 = 1.0, 1.0, 1.0, ...   <- ngrdcol entries
def duplicate_params( ngrdcol, parsed_params ):

    clubb_params = {}

    for var_name, value in parsed_params.items():
        clubb_params[var_name] = np.tile(value, ngrdcol)

    return clubb_params


# Duplicate the parameters, then change them slightly based
# on a predefined list and tweaking behavior
def duplicate_and_tweak( ngrdcol, parsed_params ):

    tweak_list = ['C7', 'C11']

    clubb_params = duplicate_params( ngrdcol, parsed_params )

    for param in tweak_list:
        # Tweak parameters in the list to have ngrdcol sample points
        # in the range initial_param_val*[1/2, 2]
        initial_param_val   = parsed_params[param]
        clubb_params[param] = initial_param_val * np.linspace( 0.5, 2, ngrdcol)

    return clubb_params



def write_multi_col_file( ngrdcol, clubb_params, output_file_name, l_multi_col_output ):
    
    with open(output_file_name, 'w') as file:
    
        file.write(f"&multicol_def\n")
        file.write(f"ngrdcol = "+str(ngrdcol)+"\n")
        file.write(f"l_output_multi_col = {'.true.' if l_multi_col_output else '.false.'}\n/\n")

        file.write(f"&clubb_params_nl\n")

        for var_name, values in clubb_params.items():
            
            # Convert the array back to a space-separated string
            values_str = ', '.join(map(str, values))
            
            # Write the variable name and duplicated values to the new file
            file.write(f"{var_name} = {values_str}\n")

        file.write(f"/\n\n")



if __name__ == "__main__":

    # Set up argument parser
    parser = argparse.ArgumentParser(description="Process CLUBB parameters.")

    parser.add_argument( "-n", type=int, required=True, help="Number of grid columns (ngrdcol)")

    parser.add_argument( "-param_file", type=str, help="Path to the CLUBB parameters file",
                         default = None )

    parser.add_argument( "-l_multi_col_output", type=str, help="Enable or disable l_output_multi_col (True/False)",
                         default = "true" )

    parser.add_argument( "-out_file", type=str, help="Enable or disable l_output_multi_col (True/False)",
                         default = "clubb_params_multi_col.in" )

    parser.add_argument( "-mode", type=str, help="Enable or disable l_output_multi_col (True/False)",
                         default = "duplicate" )

    # Parse the arguments
    args = parser.parse_args()

    # Assign parsed arguments to variables
    ngrdcol                 = args.n
    clubb_params_file       = args.param_file
    l_multi_col_output       = args.l_multi_col_output.lower() == "true"
    output_file_name        = args.out_file
    param_creation_mode     = args.mode
    
    if param_creation_mode in [ "duplicate", "dup_tweak" ] and clubb_params_file is None:
        print(f"Defining '-param_file FILE' is required  When using '-mode duplicate'")
        exit(1)
    
    if clubb_params_file is not None:
        # Parse the clubb params file
        parsed_params = parse_clubb_params(clubb_params_file)
    
    if param_creation_mode == "duplicate":
        print(f"Duplicating params from '{clubb_params_file}' '{ngrdcol}' times.")
        new_clubb_params = duplicate_params( ngrdcol, parsed_params )
    elif param_creation_mode == "dup_tweak":
        print(f"Duplicating then tweaking params from '{clubb_params_file}' '{ngrdcol}' times.")
        new_clubb_params = duplicate_and_tweak( ngrdcol, parsed_params )
    else:
        print(f"Mode '{param_creation_mode}' not recognized")
        exit(1)

    
    # Write the new_clubb_params to output_file_name
    write_multi_col_file( ngrdcol, new_clubb_params, output_file_name, l_multi_col_output )
    
    print(f"File '{output_file_name}' has been created")
    print(f"l_output_multi_col is set to '{l_multi_col_output}'")
