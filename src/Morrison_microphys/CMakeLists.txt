set(MORRISON_MICROPHYS_SOURCES 
        module_mp_graupel.F90
)

add_library(morrison_microphys_lib ${MORRISON_MICROPHYS_SOURCES})
target_link_libraries(morrison_microphys_lib PRIVATE clubb_core_lib microphys_utils_lib)
