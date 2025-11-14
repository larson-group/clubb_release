set(G_unit_tests_sources
        KK_integrals_tests.F90
        corr_cholesky_mtx_tests.F90
        hole_filling_tests.F90
        Nc_Ncn_test.F90
        read_corr_mtx_test.F90
        silhs_category_test.F90
        mu_sigma_hm_tests.F90
        pdf_parameter_tests.F90
        spurious_source_test.F90
        tuner_tests.F90
        fill_holes_tests.F90
        w_up_in_cloud_tests.F90
        smooth_heaviside_tests.F90
        smooth_min_max_tests.F90
        rev_direction_grid_test.F90
)

add_library(G_unit_test_lib ${G_unit_tests_sources})
target_link_libraries(G_unit_test_lib PRIVATE clubb_core_lib silhs_lib numerical_recipes_lib )