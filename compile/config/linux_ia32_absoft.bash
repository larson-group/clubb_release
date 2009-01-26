# $Id$
# Configuration file for a Linux machine using Absoft Fortran


# Fortran 95 compiler and linker
FC=af95 
LD=af95

# Define path to directories
dir=`pwd` # dir where this script resides
bindir="$dir/../bin"  # dir for Makefile and executable
objdir="$dir/../obj"  # dir for *.o and *.mod files
libdir="$dir/../lib"  # dir for *.a library files

# It is sometimes helpful to turn on floating-point trapping for the 
#  standalone program, but this will not work when using the tuner.
#DEBUG = -g -trap=ALL -et
DEBUG="-g"

# == Machine specific flags ==
# Note: some of these are 64 bit architectures, so make sure NetCDF is
# compiled accordingly.
ARCH="-m64 -cpu:opteron" # Opterons

# == Optimization ==
OPTIMIZE="-O3"

# == NetCDF Location ==
# Currently not working Absoft Fortran on Carson.  This is because I can't 
# find any flags that work with netCDF 3.6.2 and Absoft 10.  -dschanen 4 June 08
#NETCDF="/usr/local/netcdf-absoft64"

# == LAPACK libraries ==
#LAPACK = -L/usr/local/lib -llapack -lblas #  The netlib reference LAPACK/BLAS
#LAPACK = -L/usr/lib64 -llapack -L/usr/local/atlas/lib -lf77blas -lcblas -latlas # ATLAS BLAS (faster)
LAPACK="-L$ACML_DIR/gfortran64/lib -Wl,-rpath,$ACML_DIR/gfortran64/lib -lacml -lacml_mv -L/usr/lib64 -lgfortran"

# == Linking Flags ==
# Use -s to strip (no debugging); 
# Use -L<library path> -l<lib> to link in an external library
# Use -Wl,-rpath <library path> to set a search path for shared libs
LDFLAGS="-lU77 -L$libdir -Wl,-rpath,$libdir -lclubb_param -lclubb_bugsrad $LAPACK"

# == Compiler flags ==
# You will need to `make clean' if you change these
FFLAGS="$ARCH $DEBUG"

# Preprocessing Directives:
#   -DNETCDF enables netCDF output
#   -Dradoffline and -Dnooverlap (see bugsrad documentation)
# Define include directories. 
# Need location of include and *.mod files for the netcdf library
CPPFLAGS="-Dnooverlap -Dradoffline -DAbsoftUNIXFortran"

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
gmake="gmake -j 2"

