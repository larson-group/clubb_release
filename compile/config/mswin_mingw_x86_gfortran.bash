# $Id: linux_x86_64_gfortran.bash 4998 2011-01-27 20:44:26Z dschanen@uwm.edu $
# Configuration file for using CLUBB with MinGW and Microsoft Windows.
# This setup uses MinGW-get-inst from 16 March 2011.
# Currently the code compiles but our scripts don't work properly.

# Fortran 95 compiler and linker
FC=gfortran
LD=gfortran

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

# == Warnings ==
WARNINGS="-Wall -pedantic"

# == Machine specific flags ==
# Cygwin currently is 32 bit only
ARCH="-march=native -msse3 -mfpmath=sse -fopenmp"

# == Optimization ==
OPTIMIZE="-O2"

# == NetCDF Location ==
#NETCDF="/usr" 

# == LAPACK libraries ==
# Since there's no LAPACK in mingw, you have to compile LAPACK from source.
LAPACK="-llapack -lblas" # The netlib reference LAPACK/BLAS

# == Linking Flags ==
# Use -s to strip (no debugging); 
# Use -L<library path> -l<lib> to link in an external library
# Use -Wl,-rpath <library path> to set a search path for shared libs
#LDFLAGS="$ARCH -L$NETCDF/lib -lnetcdf -lnetcdff $LAPACK" 
LDFLAGS="$ARCH $LAPACK" 

# == Compiler flags ==
# You will need to `make clean' if you change these
FFLAGS="$ARCH $DEBUG $OPTIMIZE"

# Preprocessing Directives:
#   -DNETCDF enables netCDF output
# You will need to `make clean' if you change these
# Use -I<include path> to set a module or header file directory
#NETCDF_INCLUDE="$NETCDF/include/" 
#CPPDEFS="-DNETCDF -D__GFORTRAN__ -DCLUBB_REAL_TYPE=8"
CPPDEFS="-D__GFORTRAN__ -DCLUBB_REAL_TYPE=8"
#CPPFLAGS="-I$NETCDF_INCLUDE"

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
