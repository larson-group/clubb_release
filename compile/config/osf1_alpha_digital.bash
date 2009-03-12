# $Id$
# Configuration File for HP/Compaq/DEC Fortran on the UWM alpha

# Fortran 95 compiler and linker
FC=f95
LD=f95

# Define path to directories
dir=`pwd` # dir where this script resides
bindir="$dir/../bin"  # dir for Makefile and executable
objdir="$dir/../obj"  # dir for *.o and *.mod files
libdir="$dir/../lib"  # dir for *.a library files
srcdir="$dir/../src"  # dir where the source files reside


# It is sometimes helpful to turn on floating-point trapping for the 
#  standalone program, but this will not work when using the tuner.
# The usual options:
DEBUG="-g -check bounds -no_fp_reorder" # Digital f95

# ==  Machine specific flags ==
ARCH="-arch host"

# == Optimization ==
OPTIMIZE="-O2"

# == NetCDF Flags ==
NETCDF="$HOME/netcdf-3.6.1"

# == LAPACK libraries ==
LAPACK="-L/usr/lib -lcxml -L/usr/ccs/lib/cmplrs/cc -lexc"

# == Linking Flags ==
# Use -s to strip (no debugging); 
# Use -L<library path> -l<lib> to link in an external library
LDFLAGS="-L$libdir -Wl,-rpath,$libdir -lclubb_param -lclubb_bugsrad -L$NETCDF/lib -lnetcdf $LAPACK"

# == Compiler flags ==
# You will need to `make clean' if you change these
FFLAGS="$ARCH $OPTIMIZE"

# Preprocessing Directives:
#   -DNETCDF enables netCDF output
#   -Dradoffline and -Dnooverlap (see bugsrad documentation)
# Define include directories. 
# Need location of include and *.mod files for the netcdf library

CPPFLAGS="-DNETCDF -I$NETCDF/include -Dnooverlap -Dradoffline"

# == Static library processing ==
AR=ar
ARFLAGS=cru
RANLIB=ranlib

# == Shared library processing ==
SHARED=ld
SHAREDFLAGS="-shared -KPIC"

# Location of 'mkmf' utility
mkmf=$dir/mkmf

# gmake command to use and options: '-j 2' enables parallel compilation
gmake="gmake"

