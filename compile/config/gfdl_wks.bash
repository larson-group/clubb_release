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
source /usr/local/intel/cc/9.1.045/bin/iccvars.sh
source /usr/local/intel/fc/9.1.045/bin/ifortvars.sh

# Fortran 95 compiler and linker
FC=ifort
LD=ifort

# Define path to directories
dir=`pwd` # dir where this script resides
bindir="$dir/../bin"  # dir for Makefile and executable
objdir="$dir/../obj"  # dir for *.o and *.mod files
libdir="$dir/../lib"  # dir for *.a library files
srcdir="$dir/../src"  # dir where the source files reside

DEBUG="-g -traceback -check bounds -check uninit"
#DEBUG="-g -traceback -check bounds -check uninit -no-vec"

# == Warnings ==
WARNINGS="-warn -warn notruncated_source"

# == Machine specific options ==
ARCH=""
#ARCH="-xW" # This should work on carson/steele (AMD Opteron)
#ARCH="-xT"# Core2 Duo (overlie)

# == Optimization ==
OPTIMIZE="-O2"

# == NetCDF Location ==
NETCDF="/home/cjg/local/desktop/netcdf"

# == LAPACK libraries ==
LAPACK="-llapack -lblas"
# == Linking Flags ==
# Use -s to strip (no debugging); 
# Use -L<library path> -l<lib> to link in an external library
LDFLAGS="-L$NETCDF/lib -lnetcdf $LAPACK"

FFLAGS="$ARCH $OPTIMIZE"
#FFLAGS="$ARCH $DEBUG"

# Preprocessing Directives:
#   -DNETCDF enables netCDF output
#   -Dradoffline and -Dnooverlap (see bugsrad documentation)
# You will need to `make clean' if you change these
# Use -I<include path> to set a module or header file directory
# Define include directories. 
# Need location of include and *.mod files for the netcdf library

CPPFLAGS="-DNETCDF -I$NETCDF/include -Dnooverlap -Dradoffline"

# == Static library processing ==
AR=ar
ARFLAGS=cru
RANLIB=ranlib

# == Shared library processing ==
SHARED=$FC
SHAREDFLAGS="-fPIC -shared"

# Location of 'mkmf' utility
mkmf=$dir/mkmf

# gmake command to use and options: '-j 2' enables parallel compilation
gmake="gmake -j 2"
