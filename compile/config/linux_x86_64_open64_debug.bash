# $Id$
# Makefile definitions customized for Linux and the Open64 compiler

# Fortran 95 compiler and linker
FC=openf95
LD=openf95

# Define path to directories
dir=`pwd` # dir where this script resides
bindir="$dir/../bin"  # dir for Makefile and executable
objdir="$dir/../obj"  # dir for *.o and *.mod files
libdir="$dir/../lib"  # dir for *.a library files
srcdir="$dir/../src"  # dir where the source files reside

# == Debugging ==
DEBUG="-g -C -trapuv -DEBUG"

# == Machine specific options ==
ARCH="-msse3 -byteswapio" # Most AMD and Intel CPU's these days
#ARCH="-msse3 -mp" # Add OpenMP support

# == Used to promote all real's to double precision ==
DOUBLE_PRECISION="-r8"

# == Optimization ==
OPTIMIZE="-O0"

# == NetCDF Location ==
#NETCDF="$HOME/netcdf-3.6.3"
NETCDF="/usr/local/netcdf-open64"

# == LAPACK libraries ==
#LAPACK="-llapack -lblas"
ACML_PATH="/opt/acml5.0.0/open64_64/lib"
LAPACK="-L$ACML_PATH -Wl,-rpath $ACML_PATH -lacml"

# == Linking Flags ==
LDFLAGS="$ARCH -L$NETCDF/lib -lnetcdf $LAPACK"

# == Compiler flags ==
FFLAGS="$ARCH $OPTIMIZE $DEBUG"
#FFLAGS="$ARCH $OPTIMIZE"

# Preprocessing Directives:
#   -DNETCDF enables netCDF output
# Define include directories. 
# Need location of include and *.mod files for the netcdf library
CPPDEFS="-DNETCDF -DBYTESWAP_IO -DCLUBB_REAL_TYPE=4"
CPPFLAGS="-I$NETCDF/include"

# == Static library processing ==
AR=ar
ARFLAGS=cru
RANLIB=ranlib

# == Shared library processing ==
SHARED=$FC
SHAREDFLAGS="-KPIC -shared"

# Location of 'mkmf' utility
mkmf=$dir/mkmf

# gmake command to use and options: '-j 2' enables parallel compilation
gmake="make -j5"
