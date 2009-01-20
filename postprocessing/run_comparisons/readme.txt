Find_grads_differences.m is a MATLAB script which will compare three 
different sets of GrADS output, compare each variable within the files, 
and tell you which fields differ.

I wrote this because I was make infinitesimal changes in the code 
between each individual run and the binary data files produced were 
not identical.  All of the "usual" fields I checked (cf, rcm, thlm, 
wp2zt) were identical over the course of the run, so instead checking 
each of 400 fields by hand, I decided to write a script to do it for me.

The script itself contains documentation at its beginning.

-----

The netcdf_reader makes a multi-panel plot of different intercomparison 
cases from a passel of NetCDF files we received from Chris Golaz.  It 
is useful as an example for how to read NetCDF files using MexNC/MexCDF 
and for how to do multi-panel plots.

(For more examples of how to read data from NetCDF files in MATLAB, look 
at netcdf_rico_reader.m in the utilities/output_scripts/rico directory.  
It is not fancy; it was only written to verify [and plot] the NetCDF 
files we wrote.  For examples of how to write to NetCDF files in MATLAB, 
look at any of the netcdf_rico_writer_#.m files in that same directory.)
