add_subdirectory(Lapack_single)
add_subdirectory(Lapack_double)
add_subdirectory(Lapack_shared)

# Compose comprising lapack objects into central target
add_library(lapack_lib $<TARGET_OBJECTS:clubb_lapack_single>
                         $<TARGET_OBJECTS:clubb_lapack_double>
                         $<TARGET_OBJECTS:clubb_lapack_shared>
)
