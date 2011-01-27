# $Id$
# Configuration file for a Linux machine using the Pathscale Fortran compiler.
# Note: this configuration is no longer maintained since UWM doesn't have a
# Pathscale compiler license.  It may or may not work correctly.

# Fortran 95 compiler and linker
FC=pathf95
LD=pathf95

# Define path to directories
dir=`pwd` # dir where this script resides
bindir="$dir/../bin"  # dir for Makefile and executable
objdir="$dir/../obj"  # dir for *.o and *.mod files
libdir="$dir/../lib"  # dir for *.a library files
srcdir="$dir/../src"  # dir where the source files reside


# It is sometimes helpful to turn on floating-point trapping for the 
#  standalone program, but this will not work when using the tuner.
DEBUG="-g -ffortran-bounds-check"

# == Warnings ==
# This is the set of preferred warnings to use when compiling CLUBB.
WARNINGS="-Wall -ansi"

# == Machine specific flags ==
# Note: some of these are 64 bit architectures, so make sure NetCDF is
# compiled accordingly.
ARCH="-march=auto" # Use the native processor

# == Optimization ==
# These are all pretty conservative options. Check your compiler manual
# for information on using more aggressive techniques (inlining, etc.)
OPTIMIZE="-O3"

# == NetCDF Location ==
NETCDF="/usr/local/netcdf-pathscale"

# == LAPACK libraries ==
#LAPACK="-L/usr/local/lib -llapack -L/usr/local/atlas/lib -lf77blas -lcblas -latlas" # ATLAS BLAS (faster)
LAPACK="-L/opt/acml4.2.0/pathscale64/lib -Wl,-rpath,/opt/acml4.2.0/pathscale64/lib -lacml -lacml_mv" # AMD Core Math Library

# == Linking Flags ==
# Use -s to strip (no debugging); 
# Use -L<library path> -l<lib> to link in an external library
# Use -Wl,-rpath <library path> to set a search path for shared libs
LDFLAGS="-L$NETCDF/lib -lnetcdf $LAPACK"

# == Compiler flags ==
# You will need to `make clean' if you change these
FFLAGS="$ARCH $OPTIMIZE $DEBUG"

# Preprocessing Directives:
#   -DNETCDF enables netCDF output
#   -Dradoffline and -Dnooverlap (see bugsrad documentation)
# Define include directories. 
# Need location of include and *.mod files for the netcdf library
CPPDEFS="-DNETCDF -Dnooverlap -Dradoffline"
CPPFLAGS="-I$NETCDF/include"

# == Static library processing ==
AR=ar
ARFLAGS=cru
RANLIB=ranlib

# Location of 'mkmf' utility
mkmf=$dir/mkmf

# gmake command to use and options: '-j 2' enables parallel compilation
gmake="make"

