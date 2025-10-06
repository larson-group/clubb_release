set(MICROPHYS_UTILS microphys_stats_vars_module.F90)

add_library(microphys_utils_lib ${MICROPHYS_UTILS})
target_link_libraries(microphys_utils_lib PRIVATE clubb_core_lib)

