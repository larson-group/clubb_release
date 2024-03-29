# $Id$

# Makefile definitions customized for Linux x86_64 using the Cray compiler
# see https://github.com/larson-group/clubb/issues/1138 
# and https://docs.olcf.ornl.gov/systems/frontier_user_guide.html#id4
# and https://cpe.ext.hpe.com/docs/cce/man7/intro_openacc.7.html#compiling
# for more info. This is intended to be used with the following modules 
# module load craype-accel-amd-gfx90a cce/15.0 rocm cray-mpich cray-hdf5-parallel cray-netcdf-hdf5parallel cray-python

# Fortran 95 compiler and linker
FC=crayftn
LD=crayftn

# Define path to directories
dir=`pwd` # dir where this script resides
bindir="$dir/../bin"  # dir for Makefile and executable
objdir="$dir/../obj"  # dir for *.o and *.mod files
libdir="$dir/../lib"  # dir for *.a library files
srcdir="$dir/../src"  # dir where the source files reside

# == Debug ==
#   "-eD" gives debug info with source, great for gdb, but makes collapse() not work
DEBUG=""

# == Optimization ==
OPTIMIZE="-O2"

# == NetCDF Location ==
#NETCDF="/opt/cray/pe/netcdf/default/CRAYCLANG/14.0/"
NETCDF=$(nf-config --prefix)

# == Linking Flags ==
LDFLAGS="-L$NETCDF/lib -lnetcdff -h noomp -h acc"

FFLAGS="$OPTIMIZE $DEBUG -hnopattern -h noomp -h acc"  

# Preprocessing Directives:
#   -DNETCDF enables netCDF output
CPPDEFS="-DNETCDF -DCLUBB_REAL_TYPE=8 -DCLUBB_GPU"
CPPFLAGS="-I$NETCDF/include"

# == Static library processing ==
AR=ar
ARFLAGS=cru
RANLIB=ranlib

# == Shared library processing ==
SHARED=$FC
SHAREDFLAGS="-fPIC -shared"

# Location of 'mkmf' utility
mkmf=$dir/mkmf

# gmake command that uses "2+num_cores/2" cores for compilation
gmake="make -j$(echo "scale=0 ; 2+$(nproc) / 2.0" | bc)"
#gmake="make -j1"

