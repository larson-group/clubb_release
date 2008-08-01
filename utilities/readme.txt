The Cost_function directory contains MATLAB scripts used to calculate 
the cost function in HOC.

The output_scripts directory contains MATLAB scripts used to create the 
case-specific output files for the GABLS2, mpace_b, and rico 
intercomparisons.  These files were written by Michael Falk.

The profile_comparisons directory contains MATLAB scripts used to 
compare time-averaged profiles between LES and HOC for 18 different 
variable fields.  These files were written by Brian Griffin using some 
functions written by Michael Falk.

The run_comparisons directory contains a MATLAB script which will 
compare two runs of a given case and give you a list of which fields 
differ (which is useful assuming most of the fields are the same).  
There is also a script which reads NetCDF files into MATLAB, which may 
be useful as a model for NetCDF (MexCDF/MexNC)-based programs you may 
write.

The CLUBBStandardsCheck.pl perl script can be used to check Fortran 
source files to determine if they follow certain good software engineering 
practices which are meant to be enforced for CLUBB source files. This 
script was written by Joshua Fasching.

The differences_binary_casefiles.sh bash script can be used to quickly check
two folders containing GrADS files for binary differences.
