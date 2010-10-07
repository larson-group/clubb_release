# $Id$
# Configuration file for a Linux machine using GNU compiler collection Fortran
# Note that the version of gfortran that comes with RHEL5 cannot compile clubb.
# However, following options did work on Ubuntu 8.04 LTS (with the packaged
# versions of netcdf and netcdf-dev)


# Fortran 95 compiler and linker
FC=/usr/local/gcc/4.4.3/bin/gfortran
LD=/usr/local/gcc/4.4.3/bin/gfortran

# Define path to directories
dir=`pwd` # dir where this script resides
bindir="$dir/../bin"  # dir for Makefile and executable
objdir="$dir/../obj"  # dir for *.o and *.mod files
libdir="$dir/../lib"  # dir for *.a library files
srcdir="$dir/../src"  # dir where the source files reside


# It is sometimes helpful to turn on floating-point trapping for the 
#  standalone program, but this will not work when using the tuner.
# These are the options for debugging symbols, bounds checking & IEEE-754 
# floating point arithmetic
DEBUG="-g -fbounds-check -mieee-fp"

# == Machine specific flags ==
# Note: some of these are 64 bit architectures, so make sure NetCDF is
# compiled accordingly.
ARCH="-march=native -msse3 -mfpmath=sse -fconvert=big-endian"

# == Optimization ==
OPTIMIZE="-O2"

# == NetCDF Location ==
NETCDF="/sharedapps/LS/vlarson_group/local/netcdf-gcc"

# == LAPACK libraries ==
#LAPACK="-llapack -lblas" # The netlib reference LAPACK/BLAS
LAPACK="-L/sharedapps/LS/vlarson_group/local/atlas-gcc/lib -llapack -lf77blas -lcblas -latlas"


# == Linking Flags ==
# Use -s to strip (no debugging); 
# Use -L<library path> -l<lib> to link in an external library
# Use -Wl,-rpath <library path> to set a search path for shared libs
LDFLAGS="-L$NETCDF/lib -lnetcdf -lnetcdff $LAPACK"

# == Compiler flags ==
# You will need to `make clean' if you change these
FFLAGS="$ARCH $DEBUG"

# Preprocessing Directives:
#   -DNETCDF enables netCDF output
#   -Dradoffline and -Dnooverlap (see BUGSrad documentation)
# You will need to `make clean' if you change these
# Use -I<include path> to set a module or header file directory
CPPFLAGS="-DNETCDF -I$NETCDF/include -D__GFORTRAN__ -Dnooverlap -Dradoffline"

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
gmake="make"
