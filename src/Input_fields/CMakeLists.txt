set(INPUT_FIELDS_SOURCES
        input_fields.F90
        sounding.F90
        input_names.F90
        input_reader.F90
        stat_file_module.F90
        stat_file_utils.F90
        input_netcdf.F90
        hydrostatic_module.F90
        input_interpret.F90
        input_grads.F90
        extrapolation.F90
        corr_varnce_input_reader.F90
)

add_library(input_fields_lib ${INPUT_FIELDS_SOURCES})

target_link_libraries(input_fields_lib PRIVATE clubb_core_lib
                                               radiation_lib)

if(${ENABLE_OMP})
    target_link_libraries(input_fields_lib PRIVATE openmp-fortran)
endif()

if(${USE_NetCDF})
    target_link_libraries(input_fields_lib PRIVATE netcdf-fortran)
endif()
