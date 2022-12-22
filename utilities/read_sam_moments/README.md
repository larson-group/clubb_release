# NetCDF Variable Finder
##### Author: Tyler Cernik
##### Date: 2022-12-21
##### Version: 1.0.0

This is a simple python script that will allow us to convert variables from a NetCDF file into a CSV file. This is useful for when we want to convert a NetCDF file into a CSV file for use in other programs or analysis.

### Usage

To use this file, the first thing that you will need to do is change the file path in the variable netcdf_file in line 9.

Next you will want to copy the format of the create_df_1020_scalar_dissipation function. You would likely want to either modify it to support the variables that you want out of your file or create a new function that will do so. This function will take in a NetCDF file object, and a pandas dataframe, and fill the dataframe with the variables that you want. After that it will copy it to a csv file.

There is also a graphing function that graphs the average of the last 3 hours of variables. This may need to be modified depending on how your z and timesteps are formatted.