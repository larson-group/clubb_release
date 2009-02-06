# ------------------------------------------------------------------------------
#
# Machine specific configuration to compile CLUBB on GFDL workstations:
#
# - Intel fortran compiler version 9.1.045
# - netcdf library 3.6.3 compiled using compiled using Intel C++ and Fortran
#   compilers version 9.1.045. The library is located in
#   /home/cjg/local/desktop/netcdf
# - standard lapack and blas libraries located in /usr/lib
#
# ------------------------------------------------------------------------------

# Load desired version of Intel compiler
source /usr/local/intel/fc/9.1.045/bin/ifortvars.sh

# Fortran 90 compiler and linker
FC=ifort
LD=ifort

# Define path to directories
dir="/home/cjg/clubb/clubb_dev/scripts"      # dir where this script resides
bindir="/home/cjg/clubb/clubb_dev/code/bin"  # dir for Makefile and executable
objdir="/home/cjg/clubb/clubb_dev/code/obj"  # dir for *.o and *.mod files
libdir="/home/cjg/clubb/clubb_dev/code/lib"  # dir for *.a library files
# Define include directories. 
# Need location of include and *.mod files for the netcdf library
NETCDF="/home/cjg/local/desktop/netcdf"
CPPFLAGS="-DNETCDF -I$NETCDF/include -Dnooverlap -Dradoffline"

# Compilation flags: optimization, architecture, etc.
FFLAGS="-O2"
#FFLAGS="-g"

# Linker flags: specifiy location and libraries for CLUBB, netcdf and LAPACK
LDFLAGS="-L../lib -lclubb_bugsrad -lclubb_param -lclubb_coamps \
 -L/home/cjg/local/desktop/netcdf/lib -lnetcdf \
 -L/usr/lib -llapack -lblas"


# Directories that contain modified source code. If a file is present
# under both src/ and mods/ the one in mods/ will take precedence when
# building tha makefile. Useful for development purposes.
clubb_param_mods="/home/cjg/clubb/clubb_dev/code/src_mods/clubb_param"
clubb_standalone_mods="/home/cjg/clubb/clubb_dev/code/src_mods/clubb_standalone"

# Location of 'mkmf' utility
mkmf=$dir/mkmf

# gmake command to use and options: '-j 2' enables parallel compilation
gmake="gmake -j 2"

