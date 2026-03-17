set(CMAKE_SYSTEM_NAME       Linux)
set(CMAKE_SYSTEM_PROCESSOR  x86_64)
set(CMAKE_Fortran_COMPILER  nvfortran CACHE STRING "CLUBB Fortran Compiler")

set(WARNINGS        "")
set(COMMON_FLAGS    "-Mcache_align -Mbackslash -Mstandard")

# CMAKE_Fortran_FLAGS are default, and CMAKE_Fortran_FLAGS_DEBUG are added by setting DCMAKE_BUILD_TYPE=Debug
# in the cmake call (found in compile script, currently clubb/compile.py)
set(CMAKE_Fortran_FLAGS         "${COMMON_FLAGS} ${WARNINGS}")

# Flags that depend on compile type. 
#   - RELEASE is fast (uses -O2 because -O3 can break GPU runs)
#   - DEBUG sets debug flags and -O0 to help traceback
set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -Mstack_arrays")
set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -traceback -g -C -Kieee -Ktrap=fp")



# Define GPU options
set(GPU "none" CACHE STRING "GPU backend (none, openacc, openmp)")
set_property(CACHE GPU PROPERTY STRINGS none openacc openmp)

# Figure out which GPU flags to add
if (GPU STREQUAL "none")

    message(STATUS "Configuring for CPU build")
    set(GPU_COMPILE_FLAGS   "-noacc")

else()

    # GPU run needs CLUBB_GPU defined, and CUDA when using nvhpc
    set(GPU_DEFINITIONS CLUBB_GPU CUDA)

    # Default GPU flags
    #   -Minfo=accel enables prints of what is GPUized
    set(GPU_COMPILE_FLAGS   "-Minfo=accel")

    # Different settings are required for OpenACC or OpenMP
    if (GPU STREQUAL "openacc")
        message(STATUS "Configuring for OpenACC GPU build")
        set(GPU_COMPILE_FLAGS   "${GPU_COMPILE_FLAGS}" "-acc")
        set(GPU_LINK_FLAGS      "-acc")
    elseif (GPU STREQUAL "openmp")
        message(STATUS "Configuring for OpenMP GPU build")
        set(GPU_COMPILE_FLAGS   "${GPU_COMPILE_FLAGS}" "-mp=gpu" "-noacc")
        set(GPU_LINK_FLAGS      "-mp=gpu")
    endif()

    # Older versions of nvhpc use -Mcuda to link, but newer versions use -cuda
    # This bit extracts the version number and chooses between -Mcuda or -cuda
    execute_process(
        COMMAND ${CMAKE_Fortran_COMPILER} --version
        OUTPUT_VARIABLE NVHPC_VERSION_RAW
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    string(REGEX MATCH "[0-9]+\\.[0-9]+" NVHPC_VERSION "${NVHPC_VERSION_RAW}")

    if (NVHPC_VERSION VERSION_LESS 24.11)
        set(GPU_LINK_FLAGS  "${GPU_LINK_FLAGS}" "-Mcuda")
    else()
        set(GPU_LINK_FLAGS  "${GPU_LINK_FLAGS}" "-cuda")
    endif()

endif()



