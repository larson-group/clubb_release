# $Id$
# Makefile definitions customized for IBM AIX on the Power architecture using the XL
# Fortran compiler on Bluefire.

# Fortran 95 compiler and linker
F90=xlf90_r
F77=xlf_r
LD=xlf90_r

# Define path to directories
dir=`pwd` # dir where this script resides
bindir="$dir/../bin"  # dir for Makefile and executable
objdir="$dir/../obj"  # dir for *.o and *.mod files
libdir="$dir/../lib"  # dir for *.a library files
srcdir="$dir/../src"  # dir where the source files reside

DEBUG="-g"

# == Warnings ==
WARNINGS=""

# == Machine specific options ==
ARCH="-qnoescape" 

# == Optimization ==
OPTIMIZE="-O3"

# == NetCDF Location ==
NETCDF="/usr/local/netcdf"

# == LAPACK libraries ==
# Use Bluefire's version of LAPACK: 
LAPACK="-L/contrib/lapack/3.2.2/lib"
BLAS="-L/contrib/blas/lib"


# == Linking Flags ==
# Use -s to strip (no debugging); 
# Use -L<library path> -l<lib> to link in an external library
LDFLAGS="-L$NETCDF/lib -lnetcdf $LAPACK -llapack $BLAS -lblas"

#FFLAGS="$ARCH $DEBUG $LAPACK"
FFLAGS="$ARCH $OPTIMIZE"

# Preprocessing Directives:
#   -DNETCDF enables netCDF formet for stats output
# You will need to `make clean' if you change these
# Use -I<include path> to set a module or header file directory
# Define include directories. 
# Need location of include and *.mod files for the netcdf library
CPPFLAGS="-I$NETCDF/include"
CPPDEFS="-DNETCDF -DCLUBB_REAL_TYPE=4"

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
