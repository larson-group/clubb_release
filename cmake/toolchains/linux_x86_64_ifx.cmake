set(CMAKE_SYSTEM_NAME       Linux)
set(CMAKE_SYSTEM_PROCESSOR  x86_64)
set(CMAKE_Fortran_COMPILER  ifx CACHE STRING "CLUBB Fortran Compiler")

set(WARNINGS        "-warn -warn notruncated_source")
set(COMMON_FLAGS    "")

# CMAKE_Fortran_FLAGS are default, and CMAKE_Fortran_FLAGS_DEBUG are added by setting DCMAKE_BUILD_TYPE=Debug
# in the cmake call (found in compile script, currently clubb/compile.py)
set(CMAKE_Fortran_FLAGS         "${WARNINGS} ${COMMON_FLAGS}")

# Flags that depend on compile type. 
#   - RELEASE is fast
#   - DEBUG sets debug flags and -O0 to help traceback
set(CMAKE_Fortran_FLAGS_RELEASE "-xHost -O2")
set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -debug full -g -traceback -check bounds -fpe0 -ftz -ftrapuv -init=snan,arrays")

# Define GPU options
set(GPU "none" CACHE STRING "GPU backend (none, openacc, openmp)")
set_property(CACHE GPU PROPERTY STRINGS none openacc openmp)

# Intel's ifx can in thoery compile GPU code, but currently (as of 2025) only works for Intel GPUs
# so this has not been experimented with yet
if (NOT GPU STREQUAL "none")
    message(FATAL_ERROR "GPU compilation is not enabled for intel/ifx")
else()
    set(GPU_DEFINITIONS     "")
    set(GPU_COMPILE_FLAGS   "")
    set(GPU_LINK_FLAGS      "")
endif()



