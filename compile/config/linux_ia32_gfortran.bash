# $Id$
# Configuration file for a Linux machine using GNU compiler collection Fortran
# Note that the version of gfortran that comes with RHEL5 cannot compile clubb


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

# == Machine specific flags ==
# Note: some of these are 64 bit architectures, so make sure NetCDF is
# compiled accordingly.
#ARCH="-march=nocona -msse3 -mfpmath=sse"    # New P4s
ARCH="-march=core2 -msse3 -mfpmath=sse" # Core2 Duo's

# == Optimization ==
OPTIMIZE="-O2"

# == NetCDF Location ==
NETCDF="/usr/local/netcdf-gfortran"

# == LAPACK libraries ==
LAPACK="-L/usr/local/lib -llapack -lblas" # The netlib reference LAPACK/BLAS
#LAPACK="-L/usr/lib64 -llapack -L/usr/local/atlas/lib -lf77blas -lcblas -latlas"# ATLAS BLAS (faster)

# == Linking Flags ==
# Use -s to strip (no debugging); 
# Use -L<library path> -l<lib> to link in an external library
# Use -Wl,-rpath <library path> to set a search path for shared libs
LDFLAGS="-L$libdir -Wl,-rpath,$libdir -lclubb_param -lclubb_bugsrad -lclubb_coamps -L$NETCDF/lib -lnetcdf $LAPACK"

# == Compiler flags ==
# You will need to `make clean' if you change these
FFLAGS="$ARCH $DEBUG"

# Preprocessing Directives:
#   -DNETCDF enables netCDF output
#   -Dradoffline and -Dnooverlap (see bugsrad documentation)
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
gmake="make -j 2"

