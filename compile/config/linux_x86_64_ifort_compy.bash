# $Id$
# Makefile definitions customized for Linux x86_64 using the Intel Fortran 
# compiler 


module purge
module load intel/19.0.3
module load netcdf/4.6.3

# Fortran 95 compiler and linker
FC=ifort
LD=ifort

# Define path to directories
dir=`pwd` # dir where this script resides
bindir="$dir/../bin"  # dir for Makefile and executable
objdir="$dir/../obj"  # dir for *.o and *.mod files
libdir="$dir/../lib"  # dir for *.a library files
srcdir="$dir/../src"  # dir where the source files reside


# == Debug ==
DEBUG=""
#DEBUG="-debug full -traceback -check bounds -check uninit -fpe0 -ftz -ftrapuv"

# == Warnings ==
WARNINGS="-warn -warn notruncated_source"

# == Machine specific options ==
ARCH="-xHost" # This should work on most modern AMD/Intel computers
# == Used to promote all real's to double precision ==
DOUBLE_PRECISION="-real-size 64"

# == Optimization ==
OPTIMIZE="-O3"
#OPTIMIZE="-O3 -ipo" # Interprocedural optimization

# == NetCDF Location ==
NETCDF="/share/apps/netcdf/4.6.3/intel/19.0.3" # RHEL5

# == LAPACK libraries ==
# AMD Core Math Library
#ACML="/opt/acml5.1.0/ifort64/lib"
#LAPACK="-L$ACML -Wl,-rpath,$ACML -lacml"

# == Linking Flags ==
# Use -s to strip (no debugging); 
# Use -L<library path> -l<lib> to link in an external library
LDFLAGS="-L$NETCDF/lib -lnetcdff -lnetcdf $LAPACK"

FFLAGS="$ARCH $OPTIMIZE $DEBUG -fp-model strict"

# Preprocessing Directives:
#   -DNETCDF enables netCDF output
#   -DMKL enables MKL solver (PARDISO/GMRES) support
# You will need to `make clean' if you change these
# Use -I<include path> to set a module or header file directory
# Define include directories. 
# Need location of include and *.mod files for the netcdf library

CPPDEFS="-DNETCDF -DCLUBB_REAL_TYPE=8"
CPPFLAGS="-I$MKLPATH/../../include -I$NETCDF/include"

# == Static library processing ==
AR=ar
#AR=xiar # Interprocedural optimization archive tool
ARFLAGS=cru
RANLIB=ranlib

# == Shared library processing ==
SHARED=$FC
SHAREDFLAGS="-fPIC -shared"

# Location of 'mkmf' utility
mkmf=$dir/mkmf

# gmake command to use and options: '-j 2' enables parallel compilation
gmake="make -j5"
