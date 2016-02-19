# $Id$
# Configuration file for a Linux machine using GNU compiler collection Fortran
# Note that the version of gfortran that comes with RHEL5 (4.1.1) cannot compile clubb.
# Instead you must install gcc44-gfortran and use the command gfortran44.

# With some small modifications this configuration should work with newer 
# versions of netCDF and GNU Fortran, such as those that come with Ubuntu or
# Fedora Core.

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
#  (It works for the FIRE case as of Jan 2014)
# These are the options for debugging symbols, bounds checking & IEEE-754 
# floating point arithmetic
# Note:  floating point trapping is not currently enabled for overflow because
# of the code in l_var_covar_src.  This will be changed once all special cases
# for that code are complete.
#DEBUG="-g -fbounds-check -mieee-fp -finit-real=nan -finit-integer=-99999 -finit-logical=false -fall-intrinsics -fbacktrace" # Floating point trapping disabled
#DEBUG="-g -fbounds-check -mieee-fp -ffpe-trap=invalid,zero,overflow -finit-real=nan -finit-integer=-99999 -finit-logical=false -fall-intrinsics -fbacktrace" # Floating point trapping enabled for invalid operations, divide-by-zero, and overflow
DEBUG="-g -fbounds-check -mieee-fp -ffpe-trap=invalid,zero -finit-real=nan -finit-integer=-99999 -finit-logical=false -fall-intrinsics -fbacktrace" # Floating point trapping enabled for invalid operations and divide-by-zero

# == Warnings ==
WARNINGS="-Wall -Wextra -Wconversion -Wunderflow -Wcharacter-truncation -std=f2003 -pedantic"

# == Machine specific flags ==
ARCH="-march=native -msse3 -mfpmath=sse -fopenmp"

# == Used to promote all real's to double precision ==
DOUBLE_PRECISION="-fdefault-real-8"

# == Optimization ==
OPTIMIZE="-O3"

# == NetCDF Location ==
#NETCDF="/usr" # Ubuntu / Fedora
NETCDF="/usr/local/netcdf-gfortran" # RHEL5

# == LAPACK libraries ==
#LAPACK="-llapack -lblas" # The netlib reference LAPACK/BLAS
#LAPACK="-L/usr/lib64 -llapack -L/usr/local/atlas/lib -lf77blas -lcblas -latlas" # ATLAS BLAS (faster)
#LAPACK="-L/usr/lib/atlas-sse3 -llapack -lf77blas -lcblas -latlas" # Fedora 11 setup
ACML="/opt/acml4.4.0/gfortran64/lib"
LAPACK="-L$ACML -lacml -lacml_mv -Wl,-rpath $ACML" # AMD Core Math Library

# == Linking Flags ==
# Use -s to strip (no debugging); 
# Use -L<library path> -l<lib> to link in an external library
# Use -Wl,-rpath <library path> to set a search path for shared libs
#LDFLAGS="-L$NETCDF/lib -lnetcdf -lnetcdff $LAPACK" # Ubuntu
LDFLAGS="$ARCH -L$NETCDF/lib -lnetcdf $LAPACK" # RHEL5

# == Compiler flags ==
# You will need to `make clean' if you change these
FFLAGS="$ARCH $DEBUG $OPTIMIZE"

# Preprocessing Directives:
#   -DNETCDF enables netCDF output
# You will need to `make clean' if you change these
# Use -I<include path> to set a module or header file directory
NETCDF_INCLUDE="$NETCDF/include/" # Ubuntu 10 LTS location
#NETCDF_INCLUDE="$NETCDF/lib/gfortran/modules/" # Fedora Core 11 location
CPPDEFS="-DNETCDF -D__GFORTRAN__ -DCLUBB_REAL_TYPE=8"
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
gmake="make -j5"
