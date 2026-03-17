set(CMAKE_SYSTEM_NAME       Linux)
set(CMAKE_SYSTEM_PROCESSOR  x86_64)
set(CMAKE_Fortran_COMPILER  crayftn CACHE STRING "CLUBB Fortran Compiler")

set(WARNINGS        "")
set(COMMON_FLAGS    "-hnopattern")

# CMAKE_Fortran_FLAGS are default, and CMAKE_Fortran_FLAGS_DEBUG are added by setting DCMAKE_BUILD_TYPE=Debug
# in the cmake call (found in compile script, currently clubb/compile.py)
set(CMAKE_Fortran_FLAGS         "${COMMON_FLAGS} ${WARNINGS}")

# Flags that depend on compile type.
#   - RELEASE is fast
#   - DEBUG sets debug flags and -O0 to help traceback
#   Note: "-eD" gives source-level debug info (great for gdb) but makes OpenACC collapse() not work
set(CMAKE_Fortran_FLAGS_RELEASE "-O2")
set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -h bounds")

# Define GPU options
set(GPU "none" CACHE STRING "GPU backend (none, openacc, openmp)")
set_property(CACHE GPU PROPERTY STRINGS none openacc openmp)

if (GPU STREQUAL "none")

    message(STATUS "Configuring for CPU build")
    set(GPU_DEFINITIONS     "")
    set(GPU_COMPILE_FLAGS   "-h noacc" "-h noomp")
    set(GPU_LINK_FLAGS      "-h noacc" "-h noomp")

elseif (GPU STREQUAL "openacc")

    message(STATUS "Configuring for OpenACC GPU build")
    set(GPU_DEFINITIONS     CLUBB_GPU)
    set(GPU_COMPILE_FLAGS   "-h acc" "-h noomp")
    set(GPU_LINK_FLAGS      "-h acc" "-h noomp")

elseif (GPU STREQUAL "openmp")

    message(STATUS "Configuring for OpenMP GPU build")
    set(GPU_DEFINITIONS     CLUBB_GPU)
    set(GPU_COMPILE_FLAGS   "-h omp" "-h noacc")
    set(GPU_LINK_FLAGS      "-h omp" "-h noacc")

endif()
