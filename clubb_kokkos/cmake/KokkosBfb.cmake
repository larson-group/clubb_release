# Select the optional Fortran-backed SILHS matrix operations.
# Include after declaring the bridge and Native LAPACK targets. The requested
# cache option stays unchanged when unavailable; KOKKOS_BFB_ENABLED records the
# effective choice used by compiler definitions and build metadata.
set(KOKKOS_BFB_ENABLED OFF)
if(KOKKOS_BFB)
    if(CMAKE_Fortran_COMPILER_LOADED AND TARGET clubb_f2k_lib AND TARGET lapack_lib)
        set(KOKKOS_BFB_ENABLED ON)
    else()
        message(WARNING
            "KOKKOS_BFB=ON requested, but the Fortran/f2k LAPACK backend is unavailable. "
            "Using custom Cholesky; effective KOKKOS_BFB is OFF.")
    endif()
endif()

if(KOKKOS_BFB_ENABLED)
    message(STATUS "Kokkos Cholesky: Fortran/f2k (KOKKOS_BFB enabled)")
else()
    message(STATUS "Kokkos Cholesky: custom (KOKKOS_BFB disabled)")
endif()
