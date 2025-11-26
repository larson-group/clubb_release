set(CMAKE_SYSTEM_NAME       Linux)
set(CMAKE_SYSTEM_PROCESSOR  x86_64)
set(CMAKE_Fortran_COMPILER  gfortran CACHE STRING "CLUBB Fortran Compiler")

add_compile_definitions(__GFORTRAN__)

# Flag groups 
set(WARNINGS        "-Wall -Wextra -Wconversion -Wunderflow -Wcharacter-truncation -pedantic")
set(COMMON_FLAGS    "-mieee-fp -fall-intrinsics -std=gnu -fallow-argument-mismatch")

# CMAKE_Fortran_FLAGS are default, and CMAKE_Fortran_FLAGS_DEBUG are added by setting DCMAKE_BUILD_TYPE=Debug
# in the cmake call (found in compile script, currently clubb/compile.py)
set(CMAKE_Fortran_FLAGS         "${COMMON_FLAGS} ${WARNINGS}")

# Flags that depend on compile type. 
#   - RELEASE is fast 
#   - DEBUG sets debug flags and -O0 to help traceback
set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -march=native -msse3 -mfpmath=sse")
set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -fbounds-check -mieee-fp -ffpe-trap=invalid,zero,overflow -finit-real=snan -finit-integer=-99999 -finit-logical=false -fall-intrinsics -fbacktrace")

# Define GPU options
set(GPU "none" CACHE STRING "GPU backend (none, openacc, openmp)")
set_property(CACHE GPU PROPERTY STRINGS none openacc openmp)

# gfortran can in theory compile GPU code, but as of 2025 the capability is somewhat lacking and
# we have not tried to implement this
if (NOT GPU STREQUAL "none")
    message(FATAL_ERROR "GPU compilation is not enabled for gnu/gfortran")
else()
    set(GPU_DEFINITIONS     "")
    set(GPU_COMPILE_FLAGS   "")
    set(GPU_LINK_FLAGS      "")
endif()