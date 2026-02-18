set(CLUBB_BENCHMARK_SOURCES
        arm.F90
        arm_0003.F90
        arm_3year.F90
        arm_97.F90
        astex_a209.F90
        atex.F90
        atex_long.F90
        bomex.F90
        clex9_nov02.F90
        clex9_oct14.F90
        cloud_feedback.F90
        cobra.F90
        diag_ustar_module.F90
        dycoms2_rf01.F90
        dycoms2_rf02.F90
        fire.F90
        gabls2.F90
        gabls3.F90
        gabls3_night.F90
        jun25.F90
        lba.F90
        mpace_a.F90
        mpace_b.F90
        nov11.F90
        rico.F90
        neutral_case.F90
        spec_hum_to_mixing_ratio.F90
        sfc_flux.F90
        twp_ice.F90
        wangara.F90
        ekman.F90
        time_dependent_input.F90
        prescribe_forcings.F90
)

add_library(benchmark_cases_lib ${CLUBB_BENCHMARK_SOURCES})

target_link_libraries(benchmark_cases_lib PRIVATE clubb_core_lib input_fields_lib)

if(${ENABLE_OMP})
  target_link_libraries(benchmark_cases_lib PRIVATE openmp-fortran)
endif()
