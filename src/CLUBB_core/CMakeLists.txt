set( CLUBB_CORE_SOURCES
        adg1_adg2_3d_luhar_pdf.F90
        advance_clubb_core_module.F90
        advance_helper_module.F90
        advance_windm_edsclrm_module.F90
        advance_wp2_wp3_module.F90
        advance_xm_wpxp_module.F90
        advance_xp2_xpyp_module.F90
        advance_xp3_module.F90
        array_index.F90
        calc_pressure.F90
        calc_roots.F90
        calendar.F90
        clip_explicit.F90
        clubb_api_module.F90
        clubb_precision.F90
        code_timer_module.F90
        constants_clubb.F90
        corr_varnce_module.F90
        diagnose_correlations_module.F90
        diffusion.F90
        endian.F90
        error_code.F90
        err_info_type_module.F90
        file_functions.F90
        fill_holes.F90
        grid_class.F90
        hydromet_pdf_parameter_module.F90
        index_mapping.F90
        input_names.F90
        input_reader.F90
        interpolation.F90
        lapack_interfaces.F90
        lapack_wrap.F90
        LY93_pdf.F90
        matrix_operations.F90
        matrix_solver_wrapper.F90
        mean_adv.F90
        mixing_length.F90
        model_flags.F90
        mono_flux_limiter.F90
        Nc_Ncn_eqns.F90
        new_hybrid_pdf.F90
        new_hybrid_pdf_main.F90
        new_pdf.F90
        new_pdf_main.F90
        new_tsdadg_pdf.F90
        numerical_check.F90
        output_grads.F90
        output_netcdf.F90
        parameter_indices.F90
        parameters_tunable.F90
        pdf_closure_module.F90
        pdf_parameter_module.F90
        pdf_utilities.F90
        penta_bicgstab_solver.F90
        penta_lu_solver.F90
        pos_definite_module.F90
        precipitation_fraction.F90
        saturation.F90
        setup_clubb_pdf_params.F90
        sfc_varnce_module.F90
        sigma_sqd_w_module.F90
        Skx_module.F90
        sponge_layer_damping.F90
        stat_file_module.F90
        stats_clubb_utilities.F90
        stats_lh_sfc_module.F90
        stats_lh_zt_module.F90
        stats_rad_zm_module.F90
        stats_rad_zt_module.F90
        stats_sfc_module.F90
        stats_type.F90
        stats_type_utilities.F90
        stats_variables.F90
        stats_zm_module.F90
        stats_zt_module.F90
        T_in_K_module.F90
        tridiag_lu_solver.F90
        turbulent_adv_pdf.F90
        remapping_module.F90
        grid_adaptation_module.F90
)

add_library(clubb_core_lib ${CLUBB_CORE_SOURCES})

target_link_libraries(clubb_core_lib PRIVATE lapack_lib)

if(${USE_NetCDF})
    target_link_libraries(clubb_core_lib PRIVATE netcdf-fortran)
endif()

if(${ENABLE_OMP})
    target_link_libraries(clubb_core_lib PRIVATE openmp-fortran)
endif()

if(${USE_GPTL})
    target_link_libraries(clubb_core_lib PRIVATE gptl)
endif()