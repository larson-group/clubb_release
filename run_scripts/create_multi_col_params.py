#!/bin/python
import sys
import re
import os
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
def duplicate_and_tweak( ngrdcol, parsed_params, tweak_list, mirror=False ):

    clubb_params = duplicate_params( ngrdcol, parsed_params )

    for param in tweak_list:
        # Tweak parameters in the list to have ngrdcol sample points
        # in the range initial_param_val*[1/2, 2]
        initial_param_val   = parsed_params[param]
        clubb_params[param] = initial_param_val * np.linspace( 0.5, 2, ngrdcol)

        # If mirror is set we mirror the first half of the parameter set to the last half
        # This is mainly for testing purposes
        if mirror:
            first_half = clubb_params[param][:ngrdcol//2]
            clubb_params[param][(ngrdcol+1)//2:] = first_half[::-1]

    return clubb_params


# Create a hypergrid from a list a variables to tweak
# For exmaple:
#   -tweak_list C7,C8  and  -ngrdcol 10
# will result in a grid array of points with each parameter 
# axis being 10 possible values
def hypergrid(samples_per_param, parsed_params, tweak_list):
    ngrdcol = samples_per_param ** len(tweak_list)

    # Create points for each parameter in tweak_list
    param_grids = []
    for param in tweak_list:
        param_points = parsed_params[param] * np.linspace(0.5, 2, samples_per_param)
        param_grids.append(param_points)

    # Generate a meshgrid for all parameters
    mesh = np.meshgrid(*param_grids, indexing='ij')

    # Flatten each parameter's grid and store in clubb_params
    clubb_params = duplicate_params(ngrdcol, parsed_params)
    for i, param in enumerate(tweak_list):
        clubb_params[param][:] = mesh[i].ravel()

    return clubb_params, ngrdcol


if __name__ == "__main__":

    # Get the directory of the current script
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # Set up argument parser
    parser = argparse.ArgumentParser(description="Process CLUBB parameters.")

    parser.add_argument( "-n", type=int, required=True, help="Number of grid columns (ngrdcol)")

    tunable_parameters_default_path = os.path.join(script_dir, "../input/tunable_parameters/tunable_parameters.in")
    parser.add_argument( "-param_file", type=str, help="Path to the CLUBB parameters file",
                         default = str(tunable_parameters_default_path) )

    parser.add_argument( "-l_multi_col_output", type=str, help="Set precision or disable l_output_multi_col (no/single/double)",
                         default = "double" )

    parser.add_argument( "-out_file", type=str, help="Enable or disable l_output_multi_col (True/False)",
                         default = "clubb_params_multi_col.in" )

    parser.add_argument( "-mode", type=str, help="Enable or disable l_output_multi_col (True/False)",
                         default = "dup_tweak" )

    parser.add_argument( "-calls_per_out", type=int, help="Number of timesteps between multi_col output call",
                         default = 1 )

    parser.add_argument( "-mirror", type=str, help="mirror param lists",
                         default = "false" )

    parser.add_argument(
        "-tweak_list",
        type=lambda s: s.split(','),  # Split the input string on commas
        help="Comma-separated list of tweaks (e.g., C1,C2,C3)"
    )

    # Parse the arguments
    args = parser.parse_args()

    # Assign parsed arguments to variables
    ngrdcol                 = args.n
    clubb_params_file       = args.param_file
    l_multi_col_output      = args.l_multi_col_output
    output_file_name        = args.out_file
    param_creation_mode     = args.mode
    calls_per_out           = args.calls_per_out
    mirror                  = args.mirror == "true"

    if args.tweak_list is None:
        tweak_list = ['C7','C11']
    else:
        tweak_list = args.tweak_list
    
    # Parse the clubb params file
    parsed_params = parse_clubb_params(clubb_params_file)
    
    # Create new parameter set dependin on selected creation mode
    if param_creation_mode == "duplicate":

        print(f"Duplicating params")
        print(f" - Initial values file: '{clubb_params_file}'")
        print(f" - Duplicating '{ngrdcol}' times.")

        clubb_params = duplicate_params( ngrdcol, parsed_params )

    elif param_creation_mode == "dup_tweak":

        print(f"Duplicating then tweaking params")
        print(f" - Initial values file: '{clubb_params_file}'")
        print(f" - Duplicating '{ngrdcol}' times and tweaking: {tweak_list}")

        clubb_params = duplicate_and_tweak( ngrdcol, parsed_params, tweak_list, mirror )

    elif param_creation_mode == "hypergrid":

        print(f"Creating a hypergrid of parameters")
        print(f" - Initial values file: '{clubb_params_file}'")
        print(f" - Gridizing {tweak_list} into {ngrdcol}^{len(tweak_list)}={ngrdcol**len(tweak_list)} total points.")

        clubb_params, ngrdcol = hypergrid( ngrdcol, parsed_params, tweak_list )

    else:
        print(f"Mode '{param_creation_mode}' not recognized")
        exit(1)
    
    #--------------------  Write the clubb_params to output_file_name --------------------  
    print(f"Writing to '{output_file_name}':")
    
    with open(output_file_name, 'w') as file:
    
        file.write(f"&multicol_def\n")

        file.write(f"ngrdcol = {ngrdcol}\n")
        print(f" - ngrdcol = {ngrdcol}")

        if l_multi_col_output == "double":
            file.write(f"l_output_multi_col = .true.\n")
            file.write(f"l_output_double_prec = .true.\n")
            print(f" - l_output_multi_col = .true.")
            print(f" - l_output_double_prec = .true.")
        elif l_multi_col_output == "single":
            file.write(f"l_output_multi_col = .true.\n")
            file.write(f"l_output_double_prec = .false.\n")
            print(f" - l_output_multi_col = .true.")
            print(f" - l_output_double_prec = .false.")
        else:
            file.write(f"l_output_multi_col = .false.\n")
            print(f" - l_output_multi_col = .false.")

        file.write(f"calls_per_out = {calls_per_out}\n")
        print(f" - calls_per_out = {calls_per_out}")
        file.write(f"/\n")

        file.write(f"&clubb_params_nl\n")

        for var_name, values in clubb_params.items():
            
            # Convert the array back to a space-separated string
            values_str = ', '.join(map(str, values))
            
            # Write the variable name and duplicated values to the new file
            file.write(f"{var_name} = {values_str}\n")

        file.write(f"/\n\n")
    
