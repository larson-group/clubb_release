# $Id$
# Configuration file for using CLUBB with Cygwin and Microsoft Windows.

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
# standalone program, but this will not work when using the tuner.
# These are the options for debugging symbols, bounds checking & IEEE-754 
# floating point arithmetic
DEBUG="-g -fbounds-check -mieee-fp -finit-real=nan -finit-integer=-99999 -finit-logical=false"
#DEBUG="-g -fbounds-check -mieee-fp -ffpe-trap=invalid,zero,overflow -finit-real=nan -finit-integer=-99999 -finit-logical=false" # Floating point trapping enabled

# == Warnings ==
WARNINGS="-Wall -Wconversion -Wunderflow -Wcharacter-truncation -pedantic"

# == Machine specific flags ==
# The stack size cannot be increased in Cygwin.  Running with -fopenmp results
# in stack overflow at runtime.  Any of the following three options will allow
# CLUBB to run successfully. 
#ARCH="-march=native -msse3 -mfpmath=sse -fno-automatic"
#ARCH="-march=native -msse3 -mfpmath=sse"
ARCH="-march=native -msse3 -mfpmath=sse -fmax-stack-var-size=2000000"

# == Used to promote all real variable types to double precision ==
DOUBLE_PRECISION="-fdefault-real-8"

# == Optimization ==
OPTIMIZE="-O2"

# == NetCDF Location ==
NETCDF="/usr" 

# == LAPACK libraries ==
LAPACK="-llapack -lblas" # The netlib reference LAPACK/BLAS

# == Linking Flags ==
# Use -s to strip (no debugging); 
# Use -L<library path> -l<lib> to link in an external library
# Use -Wl,-rpath <library path> to set a search path for shared libs
LDFLAGS="$ARCH -L$NETCDF/lib -lnetcdf -lnetcdff $LAPACK" 

# == Compiler flags ==
# You will need to `make clean' if you change these
FFLAGS="$ARCH $DEBUG $OPTIMIZE"

# Preprocessing Directives:
#   -DNETCDF enables netCDF output
# You will need to `make clean' if you change these
# Use -I<include path> to set a module or header file directory
NETCDF_INCLUDE="$NETCDF/include/" 
CPPDEFS="-DNETCDF -D__GFORTRAN__ -DCLUBB_REAL_TYPE=8"
#CPPDEFS="-DNETCDF -D__GFORTRAN__ -DCLUBB_REAL_TYPE=4"
CPPFLAGS="-I$NETCDF_INCLUDE"

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
