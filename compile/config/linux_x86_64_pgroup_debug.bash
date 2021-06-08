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
DEBUG="-g -C -Kieee -Ktrap=fp -Mstandard"
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
#Variable defined in larson-group.sh, see here (https://github.com/larson-group/sys_admin/blob/master/set_larson-group_paths/larson-group.sh)
NETCDF="$PGI_NETCDF_FORTRAN"

# == LAPACK libraries ==
# Portland group usually has static versions of these
#LAPACK="-llapack -lblas"
#LAPACK="-lacml"

# == Linking Flags ==
LDFLAGS="-L$NETCDF/lib -lnetcdff $LAPACK"

FFLAGS="$ARCH $OPTIMIZE $DEBUG -Mbackslash"

# Preprocessing Directives:
#   -DNETCDF enables netCDF output
# Define include directories. 
# Need location of include and *.mod files for the netcdf library
CPPDEFS="-DNETCDF -DCLUBB_REAL_TYPE=4"
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
