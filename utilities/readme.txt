The CLUBBStandardsCheck.pl perl script can be used to check Fortran 
source files to determine if they follow certain good software engineering 
practices which are meant to be enforced for CLUBB source files. This 
script was written by Joshua Fasching.

The differences_binary_casefiles.sh bash script can be used to quickly check
two folders containing GrADS files for binary differences.

make_nc_grads_readable.m takes most NetCDF files that GrADS can not open
for whatever reason (usually "no discernable x coordinate") and converts
them to a NetCDF file GrADS can open. The resulting file will be close to
COARDS convention. The to_convert and new_file variables will need to be
adjusted for your own use. It takes around an hour to convert most NetCDF
files.
