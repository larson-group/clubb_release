set(P3_MICROPHYS_SOURCES
        module_mp_p3.F90
)

add_library(p3_microphys_lib ${P3_MICROPHYS_SOURCES})
disable_compiler_warnings(p3_microphys_lib)
