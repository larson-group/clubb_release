# $Id$
# Makefile definitions customized for Linux IA64 using the Intel Fortran 
# compiler

# Fortran 95 compiler and linker
FC=ifort
LD=ifort

# Define path to directories
dir=`pwd` # dir where this script resides
bindir="$dir/../bin"  # dir for Makefile and executable
objdir="$dir/../obj"  # dir for *.o and *.mod files
libdir="$dir/../lib"  # dir for *.a library files
srcdir="$dir/../src"  # dir where the source files reside

# == Debugging ==
DEBUG="-g -check bounds"

# == Warnings ++
WARNINGS="-warn -warn notruncated_source"

# == Machine specific options ==
ARCH="-tp p1" # Itanium
#ARCH="-tp p2"# Itanium 2

# == Optimization ==
OPTIMIZE="-O2"

# == NetCDF Location ==
NETCDF="$HOME/netcdf-3.6.3"

# == LAPACK libraries ==
# Intel Math Kernel Library (v8)
LAPACK="-L/opt/intel/mkl/8.1/lib/64 -Wl,-rpath,/opt/intel/mkl/8.1/lib/64 -lmkl_lapack64 -lmkl_i2p -lmkl -lmkl_vml_i2p -lmkl_vml -lvml -lguide -lpthread"
# Generic library
# LAPACK = -llapack -lblas -lg2c

# == Linking Flags ==
# Use -s to strip (no debugging); 
# Use -L<library path> -l<lib> to link in an external library
LDFLAGS="-L$libdir -lclubb_param -lclubb_bugsrad -L$NETCDF/lib -lnetcdf $LAPACK"

# == Compiler flags ==
# You will need to `make clean' if you change these
FFLAGS="$ARCH $DEBUG"

# Preprocessing Directives:
#   -DNETCDF enables netCDF output
#   -Dradoffline and -Dnooverlap (see bugsrad documentation)
# You will need to `make clean' if you change these
# Use -I<include path> to set a module or header file directory

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
gmake="make -j 2"

