# $Id$

# Makefile definitions customized for Linux x86_64 using the Cray compiler
# see https://github.com/larson-group/clubb/issues/1138 
# and https://cpe.ext.hpe.com/docs/cce/man7/intro_openacc.7.html#compiling
# for more info. 
#
# This files was developed mainly for Frontier (https://docs.olcf.ornl.gov/systems/frontier_user_guide.html#id4) 
# and it worked once with the following modules:
#     craype-accel-amd-gfx90a cce/15.0 rocm cray-mpich cray-hdf5-parallel cray-netcdf-hdf5parallel cray-python
# but experience has shown that there is no chance that those will load and work correctly. Not sure why. 
# The problem usually seems to be finding a version of rocm that works with cce/15.0, and it seems like 
# it changes from time to time. 
# 
# This also works on Derecho with cce/15.0.0, and no such problems happen there.
#
# IMPORTANT NOTES:
# Cray treats the acc directive "default(present)" differently than nvfortran, and these
# should be deleted before compiling for GPU. This command will do the trick:
#   find ../src -type f | xargs -P 16 -I {} sed -i "s:default(present)::g" {}
#
# Also cce/15 doesn't work with async directives as implemented with out async.py script. 
# It compiles but encouters a runtime error - a mystery yet to be investigated.


# Fortran 95 compiler and linker
FC=ftn
LD=ftn

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
#NETCDF=$(nf-config --prefix)

# == Linking Flags ==
LDFLAGS="-h noomp -h acc"

FFLAGS="$OPTIMIZE $DEBUG -hnopattern -h noomp -h acc -h list=a -h msgs -rm"  

# Preprocessing Directives:
#   -DNETCDF enables netCDF output
CPPDEFS="-DCLUBB_REAL_TYPE=8 -DCLUBB_GPU"
CPPFLAGS=""

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

