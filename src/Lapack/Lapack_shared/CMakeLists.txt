set(LAPACK_SHARED_SOURCES
                            lsame.f
                            xerbla.f
)

add_library(clubb_lapack_shared OBJECT ${LAPACK_SHARED_SOURCES})
disable_compiler_warnings(clubb_lapack_shared)
add_compile_options("-fpp")
