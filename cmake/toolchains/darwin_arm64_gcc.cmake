set(CMAKE_SYSTEM_NAME       Darwin)
set(CMAKE_SYSTEM_PROCESSOR  arm64)
set(CMAKE_SYSTEM_VERSION    "${CMAKE_HOST_SYSTEM_VERSION}")
set(CMAKE_OSX_ARCHITECTURES "arm64")
set(CMAKE_Fortran_COMPILER  gfortran CACHE STRING "CLUBB Fortran Compiler")

add_compile_definitions(__GFORTRAN__)

# Precision is configured centrally in cmake/utils/CPPDefinitions.cmake.
set(WARNINGS        "-Wall -Wextra -Wconversion -Wunderflow -Wcharacter-truncation -pedantic")
set(COMMON_FLAGS    "-fall-intrinsics -std=gnu -fallow-argument-mismatch")

# CMAKE_Fortran_FLAGS are default, and CMAKE_Fortran_FLAGS_DEBUG are added by setting DCMAKE_BUILD_TYPE=Debug
# in the cmake call (found in compile script, currently clubb/compile.py)
set(CMAKE_Fortran_FLAGS         "${COMMON_FLAGS} ${WARNINGS}")

# Flags that depend on compile type.
#   - RELEASE is fast
#   - DEBUG sets debug flags and -O0 to help traceback
set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -march=native")
set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -fbounds-check -ffpe-trap=invalid,zero,overflow -finit-real=snan -finit-integer=-99999 -finit-logical=false -fbacktrace")

# Reserve enough Mach-O header space for install_name_tool to inject rpaths at install time.
set(CMAKE_EXE_LINKER_FLAGS    "${CMAKE_EXE_LINKER_FLAGS} -Wl,-headerpad_max_install_names")
set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,-headerpad_max_install_names")
set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} -Wl,-headerpad_max_install_names")

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
