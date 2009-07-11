# $Id$
# Makefile definitions customized for Linux IA32 using the Portland Group
# compiler

# Fortran 95 compiler and linker
FC=pgf95
LD=pgf95

# Define path to directories
dir=`pwd` # dir where this script resides
bindir="$dir/../bin"  # dir for Makefile and executable
objdir="$dir/../obj"  # dir for *.o and *.mod files
libdir="$dir/../lib"  # dir for *.a library files
srcdir="$dir/../src"  # dir where the source files reside

# == Debugging ==
DEBUG="-g -C -Kieee"

# == Machine specific options ==
#ARCH="-tp piii"# PGF90, Pentium III
#ARCH="-tp p7"#	PGF90, Pentium IV
ARCH="-tp amd64" # PGF90, amd64

# == Optimization ==
OPTIMIZE="-O2"

# == NetCDF Location ==
#NETCDF="$HOME/netcdf-3.6.3"
NETCDF="/usr/local/netcdf-pgi"

# == LAPACK libraries ==
# Portland group usually has static versions of these
LAPACK="-llapack -lblas"

# == Linking Flags ==
LDFLAGS="-L$NETCDF/lib -lnetcdf $LAPACK"

FFLAGS="$ARCH $DEBUG"

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
SHARED=$FC
SHAREDFLAGS="-KPIC -shared"



# Location of 'mkmf' utility
mkmf=$dir/mkmf

# gmake command to use and options: '-j 2' enables parallel compilation
gmake="make"

