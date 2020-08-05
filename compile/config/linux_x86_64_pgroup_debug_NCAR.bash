# $Id: linux_x86_64_pgroup.bash 5834 2012-05-29 23:06:55Z charlass@uwm.edu $
# Makefile definitions customized for Linux IA32 using the Portland Group
# compiler

# Fortran 95 compiler and linker
FC=pgfortran
LD=pgfortran

# Define path to directories
dir=`pwd` # dir where this script resides
bindir="$dir/../bin"  # dir for Makefile and executable
objdir="$dir/../obj"  # dir for *.o and *.mod files
libdir="$dir/../lib"  # dir for *.a library files
srcdir="$dir/../src"  # dir where the source files reside

# == Debugging ==
DEBUG="-g -C -Kieee -Ktrap=fp"
#DEBUG="-g -C -Kieee"


# == Machine specific options ==
#ARCH="-tp piii"# PGF90, Pentium III
#ARCH="-tp p7"#	PGF90, Pentium IV
ARCH="-Mcache_align" # PGF90, amd64

# == Used to promote all real's to double precision ==
DOUBLE_PRECISION="-r8"

# == Optimization ==
# No optimization
OPTIMIZE=""
#OPTIMIZE="-O2"

# == NetCDF Location ==
#NETCDF="$HOME/netcdf-3.6.3"
NETCDF="/usr/local/netcdf-pgi"

# == LAPACK libraries ==
# Portland group usually has static versions of these
LAPACK="-DNO_LAPACK_ISNAN -llapack -lblas" #NCAR 
#LAPACK="-L/opt/acml5.0.0/pgi64 -lacml"

# == Linking Flags ==
LDFLAGS="$ARCH -L$NETCDF/lib -lnetcdf -lnetcdff $LAPACK" #NCAR

FFLAGS="$ARCH $OPTIMIZE $DEBUG -Mbackslash"

# Preprocessing Directives:
#   -DNETCDF enables netCDF output
# Define include directories. 
# Need location of include and *.mod files for the netcdf library
CPPDEFS="-DNETCDF -DNO_LAPACK_ISNAN -DCLUBB_REAL_TYPE=4" #NCAR -DNO_LAPACK_ISNAN
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

