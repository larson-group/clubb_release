set(LAPACK_SHARED_SOURCES
                            lsame.f
                            xerbla.f
)

add_library(clubb_lapack_shared OBJECT ${LAPACK_SHARED_SOURCES})
add_compile_options("-fpp")