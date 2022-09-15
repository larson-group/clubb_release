#!/usr/bin/python3

import argparse
import netCDF4
import numpy as np
import tabulate

field_threshold = 1.0e-7

parser = argparse.ArgumentParser(description='Run a test')

parser.add_argument("files", nargs=2,
                    help="need two files to diff")
                    
args = parser.parse_args()


# Create dataset from nc file
dset0 = netCDF4.Dataset(args.files[0])
dset1 = netCDF4.Dataset(args.files[1])

table = [['Var', 'Max Abs Diff', 'Max % Diff']]

for var in dset0.variables:

    if dset0[var].ndim == 3:

        # Clip fields to ignore tiny values
        field_1 = np.clip( dset0[var][:,:,:], a_min = field_threshold, a_max = 9999999.0  )
        field_2 = np.clip( dset1[var][:,:,:], a_min = field_threshold, a_max = 9999999.0 )

        max_absolute_diff = np.max(abs( field_1-field_2 ))

        # Calculate the percent difference, 100 * (a-b) / ((a+b)/2)
        max_percent_diff = np.max(200.0 * ( field_1-field_2 ) / ( field_1+field_2 )    )
                        
        table.append( [var, max_absolute_diff, max_percent_diff] )

print(tabulate.tabulate(table, headers='firstrow'))

