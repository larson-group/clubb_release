! derived_type_storage.F90 — central Fortran-derived-type storage for the F2PY interface
!
! This module ONLY stores Fortran-derived-type data and provides setters/getters.
! Routine wrappers live in separate *_wrapper.F90 files.
!
module derived_type_storage

  use clubb_precision, only: core_rknd
  use grid_class, only: grid
  use array_index, only: sclr_idx_type
  use model_flags, only: clubb_config_flags_type
  use parameters_tunable, only: nu_vertical_res_dep
  use pdf_parameter_module, only: pdf_parameter, implicit_coefs_terms
  use err_info_type_module, only: err_info_type
  use stats_netcdf, only: stats_type

  implicit none
  private

  ! Public module variables persist across F2PY calls within the same
  ! Python process. This is the bridge between Python objects and
  ! Fortran-derived-types.
  type(grid), save, public :: stored_grid
  type(sclr_idx_type), save, public :: stored_sclr_idx
  type(clubb_config_flags_type), save, public :: stored_config_flags
  type(nu_vertical_res_dep), save, public :: stored_nu_vert_res_dep
  type(pdf_parameter), save, public, target :: stored_pdf_params
  type(pdf_parameter), save, public, target :: stored_pdf_params_zm
  type(implicit_coefs_terms), save, public :: stored_pdf_implicit_coefs_terms
  type(err_info_type), save, public :: stored_err_info
  integer, save, public :: stored_err_info_chunk_idx = 1
  integer, save, public :: stored_err_info_mpi_rank = 0
  real(core_rknd), allocatable, save, public :: stored_err_info_lat(:), stored_err_info_lon(:)
  type(stats_type), save, public :: stored_stats

  public :: set_grid, get_grid
  public :: set_sclr_idx, get_sclr_idx
  public :: set_nu_vert_res_dep, get_nu_vert_res_dep
  public :: get_err_info_dims, get_err_code
  public :: get_pdf_params_dims
  public :: get_pdf_params_fields, set_pdf_params_fields
  public :: get_implicit_coefs_dims
  public :: get_implicit_coefs_fields, set_implicit_coefs_fields
  public :: get_implicit_coefs_fields_2d, set_implicit_coefs_fields_2d

contains

  !---------------------------------------------------------------------------
  ! Grid setter/getter (unchanged from Step 2)
  !---------------------------------------------------------------------------
  subroutine set_grid(nzm, nzt, ngrdcol, &
                      zm, zt, dzm, dzt, invrs_dzm, invrs_dzt, &
                      weights_zt2zm, weights_zm2zt, &
                      k_lb_zm, k_ub_zm, k_lb_zt, k_ub_zt, &
                      grid_dir_indx, grid_dir) &
    bind(C, name="set_grid_")

    integer, intent(in) :: nzm, nzt, ngrdcol
    real(core_rknd), dimension(ngrdcol, nzm), intent(in) :: zm, dzm, invrs_dzm
    real(core_rknd), dimension(ngrdcol, nzt), intent(in) :: zt, dzt, invrs_dzt
    real(core_rknd), dimension(ngrdcol, nzm, 2), intent(in) :: weights_zt2zm
    real(core_rknd), dimension(ngrdcol, nzt, 2), intent(in) :: weights_zm2zt
    integer, intent(in) :: k_lb_zm, k_ub_zm, k_lb_zt, k_ub_zt, grid_dir_indx
    real(core_rknd), intent(in) :: grid_dir

    if (allocated(stored_grid%zm))             deallocate(stored_grid%zm)
    if (allocated(stored_grid%zt))             deallocate(stored_grid%zt)
    if (allocated(stored_grid%dzm))            deallocate(stored_grid%dzm)
    if (allocated(stored_grid%dzt))            deallocate(stored_grid%dzt)
    if (allocated(stored_grid%invrs_dzm))      deallocate(stored_grid%invrs_dzm)
    if (allocated(stored_grid%invrs_dzt))      deallocate(stored_grid%invrs_dzt)
    if (allocated(stored_grid%weights_zt2zm))  deallocate(stored_grid%weights_zt2zm)
    if (allocated(stored_grid%weights_zm2zt))  deallocate(stored_grid%weights_zm2zt)

    stored_grid%nzm           = nzm
    stored_grid%nzt           = nzt
    stored_grid%k_lb_zm       = k_lb_zm
    stored_grid%k_ub_zm       = k_ub_zm
    stored_grid%k_lb_zt       = k_lb_zt
    stored_grid%k_ub_zt       = k_ub_zt
    stored_grid%grid_dir_indx = grid_dir_indx
    stored_grid%grid_dir      = grid_dir

    allocate(stored_grid%zm(ngrdcol, nzm))
    allocate(stored_grid%zt(ngrdcol, nzt))
    allocate(stored_grid%dzm(ngrdcol, nzm))
    allocate(stored_grid%dzt(ngrdcol, nzt))
    allocate(stored_grid%invrs_dzm(ngrdcol, nzm))
    allocate(stored_grid%invrs_dzt(ngrdcol, nzt))
    allocate(stored_grid%weights_zt2zm(ngrdcol, nzm, 2))
    allocate(stored_grid%weights_zm2zt(ngrdcol, nzt, 2))

    stored_grid%zm             = zm
    stored_grid%zt             = zt
    stored_grid%dzm            = dzm
    stored_grid%dzt            = dzt
    stored_grid%invrs_dzm      = invrs_dzm
    stored_grid%invrs_dzt      = invrs_dzt
    stored_grid%weights_zt2zm  = weights_zt2zm
    stored_grid%weights_zm2zt  = weights_zm2zt

  end subroutine set_grid

  subroutine get_grid(ngrdcol, nzm, nzt, &
                      zm, zt, dzm, dzt, invrs_dzm, invrs_dzt, &
                      weights_zt2zm, weights_zm2zt, &
                      k_lb_zm, k_ub_zm, k_lb_zt, k_ub_zt, &
                      grid_dir_indx, grid_dir) &
    bind(C, name="get_grid_")

    integer, intent(in)  :: ngrdcol, nzm, nzt
    real(core_rknd), dimension(ngrdcol, nzm), intent(out) :: zm, dzm, invrs_dzm
    real(core_rknd), dimension(ngrdcol, nzt), intent(out) :: zt, dzt, invrs_dzt
    real(core_rknd), dimension(ngrdcol, nzm, 2), intent(out) :: weights_zt2zm
    real(core_rknd), dimension(ngrdcol, nzt, 2), intent(out) :: weights_zm2zt
    integer, intent(out) :: k_lb_zm, k_ub_zm, k_lb_zt, k_ub_zt, grid_dir_indx
    real(core_rknd), intent(out) :: grid_dir

    zm             = stored_grid%zm
    zt             = stored_grid%zt
    dzm            = stored_grid%dzm
    dzt            = stored_grid%dzt
    invrs_dzm      = stored_grid%invrs_dzm
    invrs_dzt      = stored_grid%invrs_dzt
    weights_zt2zm  = stored_grid%weights_zt2zm
    weights_zm2zt  = stored_grid%weights_zm2zt
    k_lb_zm        = stored_grid%k_lb_zm
    k_ub_zm        = stored_grid%k_ub_zm
    k_lb_zt        = stored_grid%k_lb_zt
    k_ub_zt        = stored_grid%k_ub_zt
    grid_dir_indx  = stored_grid%grid_dir_indx
    grid_dir       = stored_grid%grid_dir

  end subroutine get_grid

  !---------------------------------------------------------------------------
  ! Scalar index setter/getter (6 integers, trivial)
  !---------------------------------------------------------------------------
  subroutine set_sclr_idx(iisclr_rt, iisclr_thl, iisclr_CO2, &
                          iiedsclr_rt, iiedsclr_thl, iiedsclr_CO2) &
    bind(C, name="set_sclr_idx_")

    integer, intent(in) :: iisclr_rt, iisclr_thl, iisclr_CO2
    integer, intent(in) :: iiedsclr_rt, iiedsclr_thl, iiedsclr_CO2

    stored_sclr_idx%iisclr_rt    = iisclr_rt
    stored_sclr_idx%iisclr_thl   = iisclr_thl
    stored_sclr_idx%iisclr_CO2   = iisclr_CO2
    stored_sclr_idx%iiedsclr_rt  = iiedsclr_rt
    stored_sclr_idx%iiedsclr_thl = iiedsclr_thl
    stored_sclr_idx%iiedsclr_CO2 = iiedsclr_CO2

  end subroutine set_sclr_idx

  subroutine get_sclr_idx(iisclr_rt, iisclr_thl, iisclr_CO2, &
                          iiedsclr_rt, iiedsclr_thl, iiedsclr_CO2) &
    bind(C, name="get_sclr_idx_")

    integer, intent(out) :: iisclr_rt, iisclr_thl, iisclr_CO2
    integer, intent(out) :: iiedsclr_rt, iiedsclr_thl, iiedsclr_CO2

    iisclr_rt    = stored_sclr_idx%iisclr_rt
    iisclr_thl   = stored_sclr_idx%iisclr_thl
    iisclr_CO2   = stored_sclr_idx%iisclr_CO2
    iiedsclr_rt  = stored_sclr_idx%iiedsclr_rt
    iiedsclr_thl = stored_sclr_idx%iiedsclr_thl
    iiedsclr_CO2 = stored_sclr_idx%iiedsclr_CO2

  end subroutine get_sclr_idx

  !---------------------------------------------------------------------------
  ! nu_vertical_res_dep setter/getter (7 allocatable 1D real arrays)
  !---------------------------------------------------------------------------
  subroutine set_nu_vert_res_dep(nzm, nu1, nu2, nu6, nu8, nu9, nu10, nu_hm) &
    bind(C, name="set_nu_vert_res_dep_")

    integer, intent(in) :: nzm
    real(core_rknd), dimension(nzm), intent(in) :: nu1, nu2, nu6, nu8, nu9, nu10, nu_hm

    if (allocated(stored_nu_vert_res_dep%nu1))  deallocate(stored_nu_vert_res_dep%nu1)
    if (allocated(stored_nu_vert_res_dep%nu2))  deallocate(stored_nu_vert_res_dep%nu2)
    if (allocated(stored_nu_vert_res_dep%nu6))  deallocate(stored_nu_vert_res_dep%nu6)
    if (allocated(stored_nu_vert_res_dep%nu8))  deallocate(stored_nu_vert_res_dep%nu8)
    if (allocated(stored_nu_vert_res_dep%nu9))  deallocate(stored_nu_vert_res_dep%nu9)
    if (allocated(stored_nu_vert_res_dep%nu10)) deallocate(stored_nu_vert_res_dep%nu10)
    if (allocated(stored_nu_vert_res_dep%nu_hm)) deallocate(stored_nu_vert_res_dep%nu_hm)

    allocate(stored_nu_vert_res_dep%nu1(nzm))
    allocate(stored_nu_vert_res_dep%nu2(nzm))
    allocate(stored_nu_vert_res_dep%nu6(nzm))
    allocate(stored_nu_vert_res_dep%nu8(nzm))
    allocate(stored_nu_vert_res_dep%nu9(nzm))
    allocate(stored_nu_vert_res_dep%nu10(nzm))
    allocate(stored_nu_vert_res_dep%nu_hm(nzm))

    stored_nu_vert_res_dep%nu1   = nu1
    stored_nu_vert_res_dep%nu2   = nu2
    stored_nu_vert_res_dep%nu6   = nu6
    stored_nu_vert_res_dep%nu8   = nu8
    stored_nu_vert_res_dep%nu9   = nu9
    stored_nu_vert_res_dep%nu10  = nu10
    stored_nu_vert_res_dep%nu_hm = nu_hm

  end subroutine set_nu_vert_res_dep

  subroutine get_nu_vert_res_dep(nzm, nu1, nu2, nu6, nu8, nu9, nu10, nu_hm) &
    bind(C, name="get_nu_vert_res_dep_")

    integer, intent(in) :: nzm
    real(core_rknd), dimension(nzm), intent(out) :: nu1, nu2, nu6, nu8, nu9, nu10, nu_hm

    nu1   = stored_nu_vert_res_dep%nu1
    nu2   = stored_nu_vert_res_dep%nu2
    nu6   = stored_nu_vert_res_dep%nu6
    nu8   = stored_nu_vert_res_dep%nu8
    nu9   = stored_nu_vert_res_dep%nu9
    nu10  = stored_nu_vert_res_dep%nu10
    nu_hm = stored_nu_vert_res_dep%nu_hm

  end subroutine get_nu_vert_res_dep

  !---------------------------------------------------------------------------
  ! err_info: expose dimension and err_code to Python
  !---------------------------------------------------------------------------
  subroutine get_err_info_dims(ngrdcol) &
    bind(C, name="get_err_info_dims_")

    integer, intent(out) :: ngrdcol

    if (allocated(stored_err_info%err_code)) then
      ngrdcol = size(stored_err_info%err_code)
    else
      ngrdcol = 0
    end if

  end subroutine get_err_info_dims

  subroutine get_err_code(ngrdcol, err_code) &
    bind(C, name="get_err_code_")

    integer, intent(in)  :: ngrdcol
    integer, intent(out) :: err_code(ngrdcol)

    err_code = stored_err_info%err_code

  end subroutine get_err_code

  !---------------------------------------------------------------------------
  ! pdf_parameter: fieldwise get/set (47 fields, all (ngrdcol, nz))
  !---------------------------------------------------------------------------

  subroutine get_pdf_params_dims(ngrdcol, nz, ngrdcol_zm, nz_zm) &
    bind(C, name="get_pdf_params_dims_")

    integer, intent(out) :: ngrdcol, nz, ngrdcol_zm, nz_zm

    ngrdcol    = stored_pdf_params%ngrdcol
    nz         = stored_pdf_params%nz
    ngrdcol_zm = stored_pdf_params_zm%ngrdcol
    nz_zm      = stored_pdf_params_zm%nz

  end subroutine get_pdf_params_dims

  subroutine get_pdf_params_fields(ngrdcol, nz, which, &
                                       w_1, w_2, varnce_w_1, varnce_w_2, &
                                       rt_1, rt_2, varnce_rt_1, varnce_rt_2, &
                                       thl_1, thl_2, varnce_thl_1, varnce_thl_2, &
                                       corr_w_rt_1, corr_w_rt_2, corr_w_thl_1, corr_w_thl_2, &
                                       corr_rt_thl_1, corr_rt_thl_2, alpha_thl, alpha_rt, &
                                       crt_1, crt_2, cthl_1, cthl_2, chi_1, chi_2, &
                                       stdev_chi_1, stdev_chi_2, stdev_eta_1, stdev_eta_2, &
                                       covar_chi_eta_1, covar_chi_eta_2, &
                                       corr_w_chi_1, corr_w_chi_2, corr_w_eta_1, corr_w_eta_2, &
                                       corr_chi_eta_1, corr_chi_eta_2, rsatl_1, rsatl_2, &
                                       rc_1, rc_2, cloud_frac_1, cloud_frac_2, mixt_frac, &
                                       ice_supersat_frac_1, ice_supersat_frac_2) &
    bind(C, name="get_pdf_params_fields_")

    integer, intent(in) :: ngrdcol, nz, which
    real(core_rknd), dimension(ngrdcol, nz), intent(out) :: &
      w_1, w_2, varnce_w_1, varnce_w_2, &
      rt_1, rt_2, varnce_rt_1, varnce_rt_2, &
      thl_1, thl_2, varnce_thl_1, varnce_thl_2, &
      corr_w_rt_1, corr_w_rt_2, corr_w_thl_1, corr_w_thl_2, &
      corr_rt_thl_1, corr_rt_thl_2, alpha_thl, alpha_rt, &
      crt_1, crt_2, cthl_1, cthl_2, chi_1, chi_2, &
      stdev_chi_1, stdev_chi_2, stdev_eta_1, stdev_eta_2, &
      covar_chi_eta_1, covar_chi_eta_2, &
      corr_w_chi_1, corr_w_chi_2, corr_w_eta_1, corr_w_eta_2, &
      corr_chi_eta_1, corr_chi_eta_2, rsatl_1, rsatl_2, &
      rc_1, rc_2, cloud_frac_1, cloud_frac_2, mixt_frac, &
      ice_supersat_frac_1, ice_supersat_frac_2

    type(pdf_parameter), pointer :: pp

    if (which == 2) then
      pp => stored_pdf_params_zm
    else
      pp => stored_pdf_params
    end if

    w_1 = pp%w_1
    w_2 = pp%w_2
    varnce_w_1 = pp%varnce_w_1
    varnce_w_2 = pp%varnce_w_2
    rt_1 = pp%rt_1
    rt_2 = pp%rt_2
    varnce_rt_1 = pp%varnce_rt_1
    varnce_rt_2 = pp%varnce_rt_2
    thl_1 = pp%thl_1
    thl_2 = pp%thl_2
    varnce_thl_1 = pp%varnce_thl_1
    varnce_thl_2 = pp%varnce_thl_2
    corr_w_rt_1 = pp%corr_w_rt_1
    corr_w_rt_2 = pp%corr_w_rt_2
    corr_w_thl_1 = pp%corr_w_thl_1
    corr_w_thl_2 = pp%corr_w_thl_2
    corr_rt_thl_1 = pp%corr_rt_thl_1
    corr_rt_thl_2 = pp%corr_rt_thl_2
    alpha_thl = pp%alpha_thl
    alpha_rt = pp%alpha_rt
    crt_1 = pp%crt_1
    crt_2 = pp%crt_2
    cthl_1 = pp%cthl_1
    cthl_2 = pp%cthl_2
    chi_1 = pp%chi_1
    chi_2 = pp%chi_2
    stdev_chi_1 = pp%stdev_chi_1
    stdev_chi_2 = pp%stdev_chi_2
    stdev_eta_1 = pp%stdev_eta_1
    stdev_eta_2 = pp%stdev_eta_2
    covar_chi_eta_1 = pp%covar_chi_eta_1
    covar_chi_eta_2 = pp%covar_chi_eta_2
    corr_w_chi_1 = pp%corr_w_chi_1
    corr_w_chi_2 = pp%corr_w_chi_2
    corr_w_eta_1 = pp%corr_w_eta_1
    corr_w_eta_2 = pp%corr_w_eta_2
    corr_chi_eta_1 = pp%corr_chi_eta_1
    corr_chi_eta_2 = pp%corr_chi_eta_2
    rsatl_1 = pp%rsatl_1
    rsatl_2 = pp%rsatl_2
    rc_1 = pp%rc_1
    rc_2 = pp%rc_2
    cloud_frac_1 = pp%cloud_frac_1
    cloud_frac_2 = pp%cloud_frac_2
    mixt_frac = pp%mixt_frac
    ice_supersat_frac_1 = pp%ice_supersat_frac_1
    ice_supersat_frac_2 = pp%ice_supersat_frac_2

  end subroutine get_pdf_params_fields

  subroutine set_pdf_params_fields(ngrdcol, nz, which, &
                                   w_1, w_2, varnce_w_1, varnce_w_2, &
                                   rt_1, rt_2, varnce_rt_1, varnce_rt_2, &
                                   thl_1, thl_2, varnce_thl_1, varnce_thl_2, &
                                   corr_w_rt_1, corr_w_rt_2, corr_w_thl_1, corr_w_thl_2, &
                                   corr_rt_thl_1, corr_rt_thl_2, alpha_thl, alpha_rt, &
                                   crt_1, crt_2, cthl_1, cthl_2, chi_1, chi_2, &
                                   stdev_chi_1, stdev_chi_2, stdev_eta_1, stdev_eta_2, &
                                   covar_chi_eta_1, covar_chi_eta_2, &
                                   corr_w_chi_1, corr_w_chi_2, corr_w_eta_1, corr_w_eta_2, &
                                   corr_chi_eta_1, corr_chi_eta_2, rsatl_1, rsatl_2, &
                                   rc_1, rc_2, cloud_frac_1, cloud_frac_2, mixt_frac, &
                                   ice_supersat_frac_1, ice_supersat_frac_2) &
    bind(C, name="set_pdf_params_fields_")

    integer, intent(in) :: ngrdcol, nz, which
    real(core_rknd), dimension(ngrdcol, nz), intent(in) :: &
      w_1, w_2, varnce_w_1, varnce_w_2, &
      rt_1, rt_2, varnce_rt_1, varnce_rt_2, &
      thl_1, thl_2, varnce_thl_1, varnce_thl_2, &
      corr_w_rt_1, corr_w_rt_2, corr_w_thl_1, corr_w_thl_2, &
      corr_rt_thl_1, corr_rt_thl_2, alpha_thl, alpha_rt, &
      crt_1, crt_2, cthl_1, cthl_2, chi_1, chi_2, &
      stdev_chi_1, stdev_chi_2, stdev_eta_1, stdev_eta_2, &
      covar_chi_eta_1, covar_chi_eta_2, &
      corr_w_chi_1, corr_w_chi_2, corr_w_eta_1, corr_w_eta_2, &
      corr_chi_eta_1, corr_chi_eta_2, rsatl_1, rsatl_2, &
      rc_1, rc_2, cloud_frac_1, cloud_frac_2, mixt_frac, &
      ice_supersat_frac_1, ice_supersat_frac_2

    type(pdf_parameter), pointer :: pp

    if (which == 2) then
      pp => stored_pdf_params_zm
    else
      pp => stored_pdf_params
    end if

    call ensure_pdf_params_allocated(pp, ngrdcol, nz)
    pp%ngrdcol = ngrdcol
    pp%nz = nz

    pp%w_1 = w_1
    pp%w_2 = w_2
    pp%varnce_w_1 = varnce_w_1
    pp%varnce_w_2 = varnce_w_2
    pp%rt_1 = rt_1
    pp%rt_2 = rt_2
    pp%varnce_rt_1 = varnce_rt_1
    pp%varnce_rt_2 = varnce_rt_2
    pp%thl_1 = thl_1
    pp%thl_2 = thl_2
    pp%varnce_thl_1 = varnce_thl_1
    pp%varnce_thl_2 = varnce_thl_2
    pp%corr_w_rt_1 = corr_w_rt_1
    pp%corr_w_rt_2 = corr_w_rt_2
    pp%corr_w_thl_1 = corr_w_thl_1
    pp%corr_w_thl_2 = corr_w_thl_2
    pp%corr_rt_thl_1 = corr_rt_thl_1
    pp%corr_rt_thl_2 = corr_rt_thl_2
    pp%alpha_thl = alpha_thl
    pp%alpha_rt = alpha_rt
    pp%crt_1 = crt_1
    pp%crt_2 = crt_2
    pp%cthl_1 = cthl_1
    pp%cthl_2 = cthl_2
    pp%chi_1 = chi_1
    pp%chi_2 = chi_2
    pp%stdev_chi_1 = stdev_chi_1
    pp%stdev_chi_2 = stdev_chi_2
    pp%stdev_eta_1 = stdev_eta_1
    pp%stdev_eta_2 = stdev_eta_2
    pp%covar_chi_eta_1 = covar_chi_eta_1
    pp%covar_chi_eta_2 = covar_chi_eta_2
    pp%corr_w_chi_1 = corr_w_chi_1
    pp%corr_w_chi_2 = corr_w_chi_2
    pp%corr_w_eta_1 = corr_w_eta_1
    pp%corr_w_eta_2 = corr_w_eta_2
    pp%corr_chi_eta_1 = corr_chi_eta_1
    pp%corr_chi_eta_2 = corr_chi_eta_2
    pp%rsatl_1 = rsatl_1
    pp%rsatl_2 = rsatl_2
    pp%rc_1 = rc_1
    pp%rc_2 = rc_2
    pp%cloud_frac_1 = cloud_frac_1
    pp%cloud_frac_2 = cloud_frac_2
    pp%mixt_frac = mixt_frac
    pp%ice_supersat_frac_1 = ice_supersat_frac_1
    pp%ice_supersat_frac_2 = ice_supersat_frac_2

  end subroutine set_pdf_params_fields

  !---------------------------------------------------------------------------
  ! implicit_coefs_terms: fieldwise get/set
  !   - 19 2D fields (ngrdcol, nz)
  !   - 8 3D fields (ngrdcol, nz, sclr_dim), only if sclr_dim > 0
  !---------------------------------------------------------------------------

  subroutine get_implicit_coefs_dims(ngrdcol, nz, sclr_dim) &
    bind(C, name="get_implicit_coefs_dims_")

    integer, intent(out) :: ngrdcol, nz, sclr_dim

    ngrdcol  = stored_pdf_implicit_coefs_terms%ngrdcol
    nz       = stored_pdf_implicit_coefs_terms%nz
    sclr_dim = stored_pdf_implicit_coefs_terms%sclr_dim

  end subroutine get_implicit_coefs_dims

  subroutine get_implicit_coefs_fields(ngrdcol, nz, sclr_dim, sclr_dim_transport, &
                                       coef_wp4_implicit, &
                                       coef_wp2rtp_implicit, term_wp2rtp_explicit, &
                                       coef_wp2thlp_implicit, term_wp2thlp_explicit, &
                                       coef_wp2up_implicit, term_wp2up_explicit, &
                                       coef_wp2vp_implicit, term_wp2vp_explicit, &
                                       coef_wprtp2_implicit, term_wprtp2_explicit, &
                                       coef_wpthlp2_implicit, term_wpthlp2_explicit, &
                                       coef_wprtpthlp_implicit, term_wprtpthlp_explicit, &
                                       coef_wpup2_implicit, term_wpup2_explicit, &
                                       coef_wpvp2_implicit, term_wpvp2_explicit, &
                                       coef_wp2sclrp_implicit, term_wp2sclrp_explicit, &
                                       coef_wpsclrp2_implicit, term_wpsclrp2_explicit, &
                                       coef_wprtpsclrp_implicit, term_wprtpsclrp_explicit, &
                                       coef_wpthlpsclrp_implicit, term_wpthlpsclrp_explicit) &
    bind(C, name="get_implicit_coefs_fields_")

    integer, intent(in) :: ngrdcol, nz, sclr_dim, sclr_dim_transport
    real(core_rknd), dimension(ngrdcol, nz), intent(out) :: &
      coef_wp4_implicit, &
      coef_wp2rtp_implicit, term_wp2rtp_explicit, &
      coef_wp2thlp_implicit, term_wp2thlp_explicit, &
      coef_wp2up_implicit, term_wp2up_explicit, &
      coef_wp2vp_implicit, term_wp2vp_explicit, &
      coef_wprtp2_implicit, term_wprtp2_explicit, &
      coef_wpthlp2_implicit, term_wpthlp2_explicit, &
      coef_wprtpthlp_implicit, term_wprtpthlp_explicit, &
      coef_wpup2_implicit, term_wpup2_explicit, &
      coef_wpvp2_implicit, term_wpvp2_explicit
    real(core_rknd), dimension(ngrdcol, nz, sclr_dim_transport), intent(out) :: &
      coef_wp2sclrp_implicit, term_wp2sclrp_explicit, &
      coef_wpsclrp2_implicit, term_wpsclrp2_explicit, &
      coef_wprtpsclrp_implicit, term_wprtpsclrp_explicit, &
      coef_wpthlpsclrp_implicit, term_wpthlpsclrp_explicit

    coef_wp4_implicit = stored_pdf_implicit_coefs_terms%coef_wp4_implicit
    coef_wp2rtp_implicit = stored_pdf_implicit_coefs_terms%coef_wp2rtp_implicit
    term_wp2rtp_explicit = stored_pdf_implicit_coefs_terms%term_wp2rtp_explicit
    coef_wp2thlp_implicit = stored_pdf_implicit_coefs_terms%coef_wp2thlp_implicit
    term_wp2thlp_explicit = stored_pdf_implicit_coefs_terms%term_wp2thlp_explicit
    coef_wp2up_implicit = stored_pdf_implicit_coefs_terms%coef_wp2up_implicit
    term_wp2up_explicit = stored_pdf_implicit_coefs_terms%term_wp2up_explicit
    coef_wp2vp_implicit = stored_pdf_implicit_coefs_terms%coef_wp2vp_implicit
    term_wp2vp_explicit = stored_pdf_implicit_coefs_terms%term_wp2vp_explicit
    coef_wprtp2_implicit = stored_pdf_implicit_coefs_terms%coef_wprtp2_implicit
    term_wprtp2_explicit = stored_pdf_implicit_coefs_terms%term_wprtp2_explicit
    coef_wpthlp2_implicit = stored_pdf_implicit_coefs_terms%coef_wpthlp2_implicit
    term_wpthlp2_explicit = stored_pdf_implicit_coefs_terms%term_wpthlp2_explicit
    coef_wprtpthlp_implicit = stored_pdf_implicit_coefs_terms%coef_wprtpthlp_implicit
    term_wprtpthlp_explicit = stored_pdf_implicit_coefs_terms%term_wprtpthlp_explicit
    coef_wpup2_implicit = stored_pdf_implicit_coefs_terms%coef_wpup2_implicit
    term_wpup2_explicit = stored_pdf_implicit_coefs_terms%term_wpup2_explicit
    coef_wpvp2_implicit = stored_pdf_implicit_coefs_terms%coef_wpvp2_implicit
    term_wpvp2_explicit = stored_pdf_implicit_coefs_terms%term_wpvp2_explicit

    coef_wp2sclrp_implicit = 0._core_rknd
    term_wp2sclrp_explicit = 0._core_rknd
    coef_wpsclrp2_implicit = 0._core_rknd
    term_wpsclrp2_explicit = 0._core_rknd
    coef_wprtpsclrp_implicit = 0._core_rknd
    term_wprtpsclrp_explicit = 0._core_rknd
    coef_wpthlpsclrp_implicit = 0._core_rknd
    term_wpthlpsclrp_explicit = 0._core_rknd

    if (sclr_dim > 0) then
      coef_wp2sclrp_implicit(:,:,1:sclr_dim) = stored_pdf_implicit_coefs_terms%coef_wp2sclrp_implicit
      term_wp2sclrp_explicit(:,:,1:sclr_dim) = stored_pdf_implicit_coefs_terms%term_wp2sclrp_explicit
      coef_wpsclrp2_implicit(:,:,1:sclr_dim) = stored_pdf_implicit_coefs_terms%coef_wpsclrp2_implicit
      term_wpsclrp2_explicit(:,:,1:sclr_dim) = stored_pdf_implicit_coefs_terms%term_wpsclrp2_explicit
      coef_wprtpsclrp_implicit(:,:,1:sclr_dim) = stored_pdf_implicit_coefs_terms%coef_wprtpsclrp_implicit
      term_wprtpsclrp_explicit(:,:,1:sclr_dim) = stored_pdf_implicit_coefs_terms%term_wprtpsclrp_explicit
      coef_wpthlpsclrp_implicit(:,:,1:sclr_dim) = stored_pdf_implicit_coefs_terms%coef_wpthlpsclrp_implicit
      term_wpthlpsclrp_explicit(:,:,1:sclr_dim) = stored_pdf_implicit_coefs_terms%term_wpthlpsclrp_explicit
    end if

  end subroutine get_implicit_coefs_fields

  subroutine get_implicit_coefs_fields_2d(ngrdcol, nz, &
                                          coef_wp4_implicit, &
                                          coef_wp2rtp_implicit, term_wp2rtp_explicit, &
                                          coef_wp2thlp_implicit, term_wp2thlp_explicit, &
                                          coef_wp2up_implicit, term_wp2up_explicit, &
                                          coef_wp2vp_implicit, term_wp2vp_explicit, &
                                          coef_wprtp2_implicit, term_wprtp2_explicit, &
                                          coef_wpthlp2_implicit, term_wpthlp2_explicit, &
                                          coef_wprtpthlp_implicit, term_wprtpthlp_explicit, &
                                          coef_wpup2_implicit, term_wpup2_explicit, &
                                          coef_wpvp2_implicit, term_wpvp2_explicit) &
    bind(C, name="get_implicit_coefs_fields_2d_")

    integer, intent(in) :: ngrdcol, nz
    real(core_rknd), dimension(ngrdcol, nz), intent(out) :: &
      coef_wp4_implicit, &
      coef_wp2rtp_implicit, term_wp2rtp_explicit, &
      coef_wp2thlp_implicit, term_wp2thlp_explicit, &
      coef_wp2up_implicit, term_wp2up_explicit, &
      coef_wp2vp_implicit, term_wp2vp_explicit, &
      coef_wprtp2_implicit, term_wprtp2_explicit, &
      coef_wpthlp2_implicit, term_wpthlp2_explicit, &
      coef_wprtpthlp_implicit, term_wprtpthlp_explicit, &
      coef_wpup2_implicit, term_wpup2_explicit, &
      coef_wpvp2_implicit, term_wpvp2_explicit

    coef_wp4_implicit = stored_pdf_implicit_coefs_terms%coef_wp4_implicit
    coef_wp2rtp_implicit = stored_pdf_implicit_coefs_terms%coef_wp2rtp_implicit
    term_wp2rtp_explicit = stored_pdf_implicit_coefs_terms%term_wp2rtp_explicit
    coef_wp2thlp_implicit = stored_pdf_implicit_coefs_terms%coef_wp2thlp_implicit
    term_wp2thlp_explicit = stored_pdf_implicit_coefs_terms%term_wp2thlp_explicit
    coef_wp2up_implicit = stored_pdf_implicit_coefs_terms%coef_wp2up_implicit
    term_wp2up_explicit = stored_pdf_implicit_coefs_terms%term_wp2up_explicit
    coef_wp2vp_implicit = stored_pdf_implicit_coefs_terms%coef_wp2vp_implicit
    term_wp2vp_explicit = stored_pdf_implicit_coefs_terms%term_wp2vp_explicit
    coef_wprtp2_implicit = stored_pdf_implicit_coefs_terms%coef_wprtp2_implicit
    term_wprtp2_explicit = stored_pdf_implicit_coefs_terms%term_wprtp2_explicit
    coef_wpthlp2_implicit = stored_pdf_implicit_coefs_terms%coef_wpthlp2_implicit
    term_wpthlp2_explicit = stored_pdf_implicit_coefs_terms%term_wpthlp2_explicit
    coef_wprtpthlp_implicit = stored_pdf_implicit_coefs_terms%coef_wprtpthlp_implicit
    term_wprtpthlp_explicit = stored_pdf_implicit_coefs_terms%term_wprtpthlp_explicit
    coef_wpup2_implicit = stored_pdf_implicit_coefs_terms%coef_wpup2_implicit
    term_wpup2_explicit = stored_pdf_implicit_coefs_terms%term_wpup2_explicit
    coef_wpvp2_implicit = stored_pdf_implicit_coefs_terms%coef_wpvp2_implicit
    term_wpvp2_explicit = stored_pdf_implicit_coefs_terms%term_wpvp2_explicit
  end subroutine get_implicit_coefs_fields_2d

  subroutine set_implicit_coefs_fields(ngrdcol, nz, sclr_dim, sclr_dim_transport, &
                                       coef_wp4_implicit, &
                                       coef_wp2rtp_implicit, term_wp2rtp_explicit, &
                                       coef_wp2thlp_implicit, term_wp2thlp_explicit, &
                                       coef_wp2up_implicit, term_wp2up_explicit, &
                                       coef_wp2vp_implicit, term_wp2vp_explicit, &
                                       coef_wprtp2_implicit, term_wprtp2_explicit, &
                                       coef_wpthlp2_implicit, term_wpthlp2_explicit, &
                                       coef_wprtpthlp_implicit, term_wprtpthlp_explicit, &
                                       coef_wpup2_implicit, term_wpup2_explicit, &
                                       coef_wpvp2_implicit, term_wpvp2_explicit, &
                                       coef_wp2sclrp_implicit, term_wp2sclrp_explicit, &
                                       coef_wpsclrp2_implicit, term_wpsclrp2_explicit, &
                                       coef_wprtpsclrp_implicit, term_wprtpsclrp_explicit, &
                                       coef_wpthlpsclrp_implicit, term_wpthlpsclrp_explicit) &
    bind(C, name="set_implicit_coefs_fields_")

    integer, intent(in) :: ngrdcol, nz, sclr_dim, sclr_dim_transport
    real(core_rknd), dimension(ngrdcol, nz), intent(in) :: &
      coef_wp4_implicit, &
      coef_wp2rtp_implicit, term_wp2rtp_explicit, &
      coef_wp2thlp_implicit, term_wp2thlp_explicit, &
      coef_wp2up_implicit, term_wp2up_explicit, &
      coef_wp2vp_implicit, term_wp2vp_explicit, &
      coef_wprtp2_implicit, term_wprtp2_explicit, &
      coef_wpthlp2_implicit, term_wpthlp2_explicit, &
      coef_wprtpthlp_implicit, term_wprtpthlp_explicit, &
      coef_wpup2_implicit, term_wpup2_explicit, &
      coef_wpvp2_implicit, term_wpvp2_explicit
    real(core_rknd), dimension(ngrdcol, nz, sclr_dim_transport), intent(in) :: &
      coef_wp2sclrp_implicit, term_wp2sclrp_explicit, &
      coef_wpsclrp2_implicit, term_wpsclrp2_explicit, &
      coef_wprtpsclrp_implicit, term_wprtpsclrp_explicit, &
      coef_wpthlpsclrp_implicit, term_wpthlpsclrp_explicit

    call ensure_implicit_2d_allocated(stored_pdf_implicit_coefs_terms, ngrdcol, nz)
    call ensure_implicit_3d_allocated(stored_pdf_implicit_coefs_terms, ngrdcol, nz, sclr_dim)
    stored_pdf_implicit_coefs_terms%ngrdcol = ngrdcol
    stored_pdf_implicit_coefs_terms%nz = nz
    stored_pdf_implicit_coefs_terms%sclr_dim = sclr_dim

    stored_pdf_implicit_coefs_terms%coef_wp4_implicit = coef_wp4_implicit
    stored_pdf_implicit_coefs_terms%coef_wp2rtp_implicit = coef_wp2rtp_implicit
    stored_pdf_implicit_coefs_terms%term_wp2rtp_explicit = term_wp2rtp_explicit
    stored_pdf_implicit_coefs_terms%coef_wp2thlp_implicit = coef_wp2thlp_implicit
    stored_pdf_implicit_coefs_terms%term_wp2thlp_explicit = term_wp2thlp_explicit
    stored_pdf_implicit_coefs_terms%coef_wp2up_implicit = coef_wp2up_implicit
    stored_pdf_implicit_coefs_terms%term_wp2up_explicit = term_wp2up_explicit
    stored_pdf_implicit_coefs_terms%coef_wp2vp_implicit = coef_wp2vp_implicit
    stored_pdf_implicit_coefs_terms%term_wp2vp_explicit = term_wp2vp_explicit
    stored_pdf_implicit_coefs_terms%coef_wprtp2_implicit = coef_wprtp2_implicit
    stored_pdf_implicit_coefs_terms%term_wprtp2_explicit = term_wprtp2_explicit
    stored_pdf_implicit_coefs_terms%coef_wpthlp2_implicit = coef_wpthlp2_implicit
    stored_pdf_implicit_coefs_terms%term_wpthlp2_explicit = term_wpthlp2_explicit
    stored_pdf_implicit_coefs_terms%coef_wprtpthlp_implicit = coef_wprtpthlp_implicit
    stored_pdf_implicit_coefs_terms%term_wprtpthlp_explicit = term_wprtpthlp_explicit
    stored_pdf_implicit_coefs_terms%coef_wpup2_implicit = coef_wpup2_implicit
    stored_pdf_implicit_coefs_terms%term_wpup2_explicit = term_wpup2_explicit
    stored_pdf_implicit_coefs_terms%coef_wpvp2_implicit = coef_wpvp2_implicit
    stored_pdf_implicit_coefs_terms%term_wpvp2_explicit = term_wpvp2_explicit

    if (sclr_dim > 0) then
      stored_pdf_implicit_coefs_terms%coef_wp2sclrp_implicit = coef_wp2sclrp_implicit
      stored_pdf_implicit_coefs_terms%term_wp2sclrp_explicit = term_wp2sclrp_explicit
      stored_pdf_implicit_coefs_terms%coef_wpsclrp2_implicit = coef_wpsclrp2_implicit
      stored_pdf_implicit_coefs_terms%term_wpsclrp2_explicit = term_wpsclrp2_explicit
      stored_pdf_implicit_coefs_terms%coef_wprtpsclrp_implicit = coef_wprtpsclrp_implicit
      stored_pdf_implicit_coefs_terms%term_wprtpsclrp_explicit = term_wprtpsclrp_explicit
      stored_pdf_implicit_coefs_terms%coef_wpthlpsclrp_implicit = coef_wpthlpsclrp_implicit
      stored_pdf_implicit_coefs_terms%term_wpthlpsclrp_explicit = term_wpthlpsclrp_explicit
    end if

  end subroutine set_implicit_coefs_fields

  subroutine set_implicit_coefs_fields_2d(ngrdcol, nz, &
                                          coef_wp4_implicit, &
                                          coef_wp2rtp_implicit, term_wp2rtp_explicit, &
                                          coef_wp2thlp_implicit, term_wp2thlp_explicit, &
                                          coef_wp2up_implicit, term_wp2up_explicit, &
                                          coef_wp2vp_implicit, term_wp2vp_explicit, &
                                          coef_wprtp2_implicit, term_wprtp2_explicit, &
                                          coef_wpthlp2_implicit, term_wpthlp2_explicit, &
                                          coef_wprtpthlp_implicit, term_wprtpthlp_explicit, &
                                          coef_wpup2_implicit, term_wpup2_explicit, &
                                          coef_wpvp2_implicit, term_wpvp2_explicit) &
    bind(C, name="set_implicit_coefs_fields_2d_")

    integer, intent(in) :: ngrdcol, nz
    real(core_rknd), dimension(ngrdcol, nz), intent(in) :: &
      coef_wp4_implicit, &
      coef_wp2rtp_implicit, term_wp2rtp_explicit, &
      coef_wp2thlp_implicit, term_wp2thlp_explicit, &
      coef_wp2up_implicit, term_wp2up_explicit, &
      coef_wp2vp_implicit, term_wp2vp_explicit, &
      coef_wprtp2_implicit, term_wprtp2_explicit, &
      coef_wpthlp2_implicit, term_wpthlp2_explicit, &
      coef_wprtpthlp_implicit, term_wprtpthlp_explicit, &
      coef_wpup2_implicit, term_wpup2_explicit, &
      coef_wpvp2_implicit, term_wpvp2_explicit

    call ensure_implicit_2d_allocated(stored_pdf_implicit_coefs_terms, ngrdcol, nz)
    call ensure_implicit_3d_allocated(stored_pdf_implicit_coefs_terms, ngrdcol, nz, 0)
    stored_pdf_implicit_coefs_terms%ngrdcol = ngrdcol
    stored_pdf_implicit_coefs_terms%nz = nz
    stored_pdf_implicit_coefs_terms%sclr_dim = 0

    stored_pdf_implicit_coefs_terms%coef_wp4_implicit = coef_wp4_implicit
    stored_pdf_implicit_coefs_terms%coef_wp2rtp_implicit = coef_wp2rtp_implicit
    stored_pdf_implicit_coefs_terms%term_wp2rtp_explicit = term_wp2rtp_explicit
    stored_pdf_implicit_coefs_terms%coef_wp2thlp_implicit = coef_wp2thlp_implicit
    stored_pdf_implicit_coefs_terms%term_wp2thlp_explicit = term_wp2thlp_explicit
    stored_pdf_implicit_coefs_terms%coef_wp2up_implicit = coef_wp2up_implicit
    stored_pdf_implicit_coefs_terms%term_wp2up_explicit = term_wp2up_explicit
    stored_pdf_implicit_coefs_terms%coef_wp2vp_implicit = coef_wp2vp_implicit
    stored_pdf_implicit_coefs_terms%term_wp2vp_explicit = term_wp2vp_explicit
    stored_pdf_implicit_coefs_terms%coef_wprtp2_implicit = coef_wprtp2_implicit
    stored_pdf_implicit_coefs_terms%term_wprtp2_explicit = term_wprtp2_explicit
    stored_pdf_implicit_coefs_terms%coef_wpthlp2_implicit = coef_wpthlp2_implicit
    stored_pdf_implicit_coefs_terms%term_wpthlp2_explicit = term_wpthlp2_explicit
    stored_pdf_implicit_coefs_terms%coef_wprtpthlp_implicit = coef_wprtpthlp_implicit
    stored_pdf_implicit_coefs_terms%term_wprtpthlp_explicit = term_wprtpthlp_explicit
    stored_pdf_implicit_coefs_terms%coef_wpup2_implicit = coef_wpup2_implicit
    stored_pdf_implicit_coefs_terms%term_wpup2_explicit = term_wpup2_explicit
    stored_pdf_implicit_coefs_terms%coef_wpvp2_implicit = coef_wpvp2_implicit
    stored_pdf_implicit_coefs_terms%term_wpvp2_explicit = term_wpvp2_explicit
  end subroutine set_implicit_coefs_fields_2d

  subroutine ensure_pdf_params_allocated(pp, ngrdcol, nz)
    type(pdf_parameter), intent(inout) :: pp
    integer, intent(in) :: ngrdcol, nz
    logical :: needs_alloc

    needs_alloc = .not.allocated(pp%w_1)
    if (.not.needs_alloc) then
      if (size(pp%w_1, 1) /= ngrdcol .or. size(pp%w_1, 2) /= nz) needs_alloc = .true.
    end if

    if (needs_alloc) then
      if (allocated(pp%w_1)) deallocate(pp%w_1)
      if (allocated(pp%w_2)) deallocate(pp%w_2)
      if (allocated(pp%varnce_w_1)) deallocate(pp%varnce_w_1)
      if (allocated(pp%varnce_w_2)) deallocate(pp%varnce_w_2)
      if (allocated(pp%rt_1)) deallocate(pp%rt_1)
      if (allocated(pp%rt_2)) deallocate(pp%rt_2)
      if (allocated(pp%varnce_rt_1)) deallocate(pp%varnce_rt_1)
      if (allocated(pp%varnce_rt_2)) deallocate(pp%varnce_rt_2)
      if (allocated(pp%thl_1)) deallocate(pp%thl_1)
      if (allocated(pp%thl_2)) deallocate(pp%thl_2)
      if (allocated(pp%varnce_thl_1)) deallocate(pp%varnce_thl_1)
      if (allocated(pp%varnce_thl_2)) deallocate(pp%varnce_thl_2)
      if (allocated(pp%corr_w_rt_1)) deallocate(pp%corr_w_rt_1)
      if (allocated(pp%corr_w_rt_2)) deallocate(pp%corr_w_rt_2)
      if (allocated(pp%corr_w_thl_1)) deallocate(pp%corr_w_thl_1)
      if (allocated(pp%corr_w_thl_2)) deallocate(pp%corr_w_thl_2)
      if (allocated(pp%corr_rt_thl_1)) deallocate(pp%corr_rt_thl_1)
      if (allocated(pp%corr_rt_thl_2)) deallocate(pp%corr_rt_thl_2)
      if (allocated(pp%alpha_thl)) deallocate(pp%alpha_thl)
      if (allocated(pp%alpha_rt)) deallocate(pp%alpha_rt)
      if (allocated(pp%crt_1)) deallocate(pp%crt_1)
      if (allocated(pp%crt_2)) deallocate(pp%crt_2)
      if (allocated(pp%cthl_1)) deallocate(pp%cthl_1)
      if (allocated(pp%cthl_2)) deallocate(pp%cthl_2)
      if (allocated(pp%chi_1)) deallocate(pp%chi_1)
      if (allocated(pp%chi_2)) deallocate(pp%chi_2)
      if (allocated(pp%stdev_chi_1)) deallocate(pp%stdev_chi_1)
      if (allocated(pp%stdev_chi_2)) deallocate(pp%stdev_chi_2)
      if (allocated(pp%stdev_eta_1)) deallocate(pp%stdev_eta_1)
      if (allocated(pp%stdev_eta_2)) deallocate(pp%stdev_eta_2)
      if (allocated(pp%covar_chi_eta_1)) deallocate(pp%covar_chi_eta_1)
      if (allocated(pp%covar_chi_eta_2)) deallocate(pp%covar_chi_eta_2)
      if (allocated(pp%corr_w_chi_1)) deallocate(pp%corr_w_chi_1)
      if (allocated(pp%corr_w_chi_2)) deallocate(pp%corr_w_chi_2)
      if (allocated(pp%corr_w_eta_1)) deallocate(pp%corr_w_eta_1)
      if (allocated(pp%corr_w_eta_2)) deallocate(pp%corr_w_eta_2)
      if (allocated(pp%corr_chi_eta_1)) deallocate(pp%corr_chi_eta_1)
      if (allocated(pp%corr_chi_eta_2)) deallocate(pp%corr_chi_eta_2)
      if (allocated(pp%rsatl_1)) deallocate(pp%rsatl_1)
      if (allocated(pp%rsatl_2)) deallocate(pp%rsatl_2)
      if (allocated(pp%rc_1)) deallocate(pp%rc_1)
      if (allocated(pp%rc_2)) deallocate(pp%rc_2)
      if (allocated(pp%cloud_frac_1)) deallocate(pp%cloud_frac_1)
      if (allocated(pp%cloud_frac_2)) deallocate(pp%cloud_frac_2)
      if (allocated(pp%mixt_frac)) deallocate(pp%mixt_frac)
      if (allocated(pp%ice_supersat_frac_1)) deallocate(pp%ice_supersat_frac_1)
      if (allocated(pp%ice_supersat_frac_2)) deallocate(pp%ice_supersat_frac_2)

      allocate(pp%w_1(ngrdcol, nz))
      allocate(pp%w_2(ngrdcol, nz))
      allocate(pp%varnce_w_1(ngrdcol, nz))
      allocate(pp%varnce_w_2(ngrdcol, nz))
      allocate(pp%rt_1(ngrdcol, nz))
      allocate(pp%rt_2(ngrdcol, nz))
      allocate(pp%varnce_rt_1(ngrdcol, nz))
      allocate(pp%varnce_rt_2(ngrdcol, nz))
      allocate(pp%thl_1(ngrdcol, nz))
      allocate(pp%thl_2(ngrdcol, nz))
      allocate(pp%varnce_thl_1(ngrdcol, nz))
      allocate(pp%varnce_thl_2(ngrdcol, nz))
      allocate(pp%corr_w_rt_1(ngrdcol, nz))
      allocate(pp%corr_w_rt_2(ngrdcol, nz))
      allocate(pp%corr_w_thl_1(ngrdcol, nz))
      allocate(pp%corr_w_thl_2(ngrdcol, nz))
      allocate(pp%corr_rt_thl_1(ngrdcol, nz))
      allocate(pp%corr_rt_thl_2(ngrdcol, nz))
      allocate(pp%alpha_thl(ngrdcol, nz))
      allocate(pp%alpha_rt(ngrdcol, nz))
      allocate(pp%crt_1(ngrdcol, nz))
      allocate(pp%crt_2(ngrdcol, nz))
      allocate(pp%cthl_1(ngrdcol, nz))
      allocate(pp%cthl_2(ngrdcol, nz))
      allocate(pp%chi_1(ngrdcol, nz))
      allocate(pp%chi_2(ngrdcol, nz))
      allocate(pp%stdev_chi_1(ngrdcol, nz))
      allocate(pp%stdev_chi_2(ngrdcol, nz))
      allocate(pp%stdev_eta_1(ngrdcol, nz))
      allocate(pp%stdev_eta_2(ngrdcol, nz))
      allocate(pp%covar_chi_eta_1(ngrdcol, nz))
      allocate(pp%covar_chi_eta_2(ngrdcol, nz))
      allocate(pp%corr_w_chi_1(ngrdcol, nz))
      allocate(pp%corr_w_chi_2(ngrdcol, nz))
      allocate(pp%corr_w_eta_1(ngrdcol, nz))
      allocate(pp%corr_w_eta_2(ngrdcol, nz))
      allocate(pp%corr_chi_eta_1(ngrdcol, nz))
      allocate(pp%corr_chi_eta_2(ngrdcol, nz))
      allocate(pp%rsatl_1(ngrdcol, nz))
      allocate(pp%rsatl_2(ngrdcol, nz))
      allocate(pp%rc_1(ngrdcol, nz))
      allocate(pp%rc_2(ngrdcol, nz))
      allocate(pp%cloud_frac_1(ngrdcol, nz))
      allocate(pp%cloud_frac_2(ngrdcol, nz))
      allocate(pp%mixt_frac(ngrdcol, nz))
      allocate(pp%ice_supersat_frac_1(ngrdcol, nz))
      allocate(pp%ice_supersat_frac_2(ngrdcol, nz))
    end if
  end subroutine ensure_pdf_params_allocated

  subroutine ensure_implicit_2d_allocated(ic, ngrdcol, nz)
    type(implicit_coefs_terms), intent(inout) :: ic
    integer, intent(in) :: ngrdcol, nz
    logical :: needs_alloc

    needs_alloc = .not.allocated(ic%coef_wp4_implicit)
    if (.not.needs_alloc) then
      if (size(ic%coef_wp4_implicit, 1) /= ngrdcol .or. &
          size(ic%coef_wp4_implicit, 2) /= nz) needs_alloc = .true.
    end if

    if (needs_alloc) then
      if (allocated(ic%coef_wp4_implicit)) deallocate(ic%coef_wp4_implicit)
      if (allocated(ic%coef_wp2rtp_implicit)) deallocate(ic%coef_wp2rtp_implicit)
      if (allocated(ic%term_wp2rtp_explicit)) deallocate(ic%term_wp2rtp_explicit)
      if (allocated(ic%coef_wp2thlp_implicit)) deallocate(ic%coef_wp2thlp_implicit)
      if (allocated(ic%term_wp2thlp_explicit)) deallocate(ic%term_wp2thlp_explicit)
      if (allocated(ic%coef_wp2up_implicit)) deallocate(ic%coef_wp2up_implicit)
      if (allocated(ic%term_wp2up_explicit)) deallocate(ic%term_wp2up_explicit)
      if (allocated(ic%coef_wp2vp_implicit)) deallocate(ic%coef_wp2vp_implicit)
      if (allocated(ic%term_wp2vp_explicit)) deallocate(ic%term_wp2vp_explicit)
      if (allocated(ic%coef_wprtp2_implicit)) deallocate(ic%coef_wprtp2_implicit)
      if (allocated(ic%term_wprtp2_explicit)) deallocate(ic%term_wprtp2_explicit)
      if (allocated(ic%coef_wpthlp2_implicit)) deallocate(ic%coef_wpthlp2_implicit)
      if (allocated(ic%term_wpthlp2_explicit)) deallocate(ic%term_wpthlp2_explicit)
      if (allocated(ic%coef_wprtpthlp_implicit)) deallocate(ic%coef_wprtpthlp_implicit)
      if (allocated(ic%term_wprtpthlp_explicit)) deallocate(ic%term_wprtpthlp_explicit)
      if (allocated(ic%coef_wpup2_implicit)) deallocate(ic%coef_wpup2_implicit)
      if (allocated(ic%term_wpup2_explicit)) deallocate(ic%term_wpup2_explicit)
      if (allocated(ic%coef_wpvp2_implicit)) deallocate(ic%coef_wpvp2_implicit)
      if (allocated(ic%term_wpvp2_explicit)) deallocate(ic%term_wpvp2_explicit)

      allocate(ic%coef_wp4_implicit(ngrdcol, nz))
      allocate(ic%coef_wp2rtp_implicit(ngrdcol, nz))
      allocate(ic%term_wp2rtp_explicit(ngrdcol, nz))
      allocate(ic%coef_wp2thlp_implicit(ngrdcol, nz))
      allocate(ic%term_wp2thlp_explicit(ngrdcol, nz))
      allocate(ic%coef_wp2up_implicit(ngrdcol, nz))
      allocate(ic%term_wp2up_explicit(ngrdcol, nz))
      allocate(ic%coef_wp2vp_implicit(ngrdcol, nz))
      allocate(ic%term_wp2vp_explicit(ngrdcol, nz))
      allocate(ic%coef_wprtp2_implicit(ngrdcol, nz))
      allocate(ic%term_wprtp2_explicit(ngrdcol, nz))
      allocate(ic%coef_wpthlp2_implicit(ngrdcol, nz))
      allocate(ic%term_wpthlp2_explicit(ngrdcol, nz))
      allocate(ic%coef_wprtpthlp_implicit(ngrdcol, nz))
      allocate(ic%term_wprtpthlp_explicit(ngrdcol, nz))
      allocate(ic%coef_wpup2_implicit(ngrdcol, nz))
      allocate(ic%term_wpup2_explicit(ngrdcol, nz))
      allocate(ic%coef_wpvp2_implicit(ngrdcol, nz))
      allocate(ic%term_wpvp2_explicit(ngrdcol, nz))
    end if
  end subroutine ensure_implicit_2d_allocated

  subroutine ensure_implicit_3d_allocated(ic, ngrdcol, nz, sclr_dim)
    type(implicit_coefs_terms), intent(inout) :: ic
    integer, intent(in) :: ngrdcol, nz, sclr_dim
    logical :: needs_alloc

    if (sclr_dim <= 0) then
      if (allocated(ic%coef_wp2sclrp_implicit)) deallocate(ic%coef_wp2sclrp_implicit)
      if (allocated(ic%term_wp2sclrp_explicit)) deallocate(ic%term_wp2sclrp_explicit)
      if (allocated(ic%coef_wpsclrp2_implicit)) deallocate(ic%coef_wpsclrp2_implicit)
      if (allocated(ic%term_wpsclrp2_explicit)) deallocate(ic%term_wpsclrp2_explicit)
      if (allocated(ic%coef_wprtpsclrp_implicit)) deallocate(ic%coef_wprtpsclrp_implicit)
      if (allocated(ic%term_wprtpsclrp_explicit)) deallocate(ic%term_wprtpsclrp_explicit)
      if (allocated(ic%coef_wpthlpsclrp_implicit)) deallocate(ic%coef_wpthlpsclrp_implicit)
      if (allocated(ic%term_wpthlpsclrp_explicit)) deallocate(ic%term_wpthlpsclrp_explicit)
      return
    end if

    needs_alloc = .not.allocated(ic%coef_wp2sclrp_implicit)
    if (.not.needs_alloc) then
      if (size(ic%coef_wp2sclrp_implicit, 1) /= ngrdcol .or. &
          size(ic%coef_wp2sclrp_implicit, 2) /= nz .or. &
          size(ic%coef_wp2sclrp_implicit, 3) /= sclr_dim) needs_alloc = .true.
    end if

    if (needs_alloc) then
      if (allocated(ic%coef_wp2sclrp_implicit)) deallocate(ic%coef_wp2sclrp_implicit)
      if (allocated(ic%term_wp2sclrp_explicit)) deallocate(ic%term_wp2sclrp_explicit)
      if (allocated(ic%coef_wpsclrp2_implicit)) deallocate(ic%coef_wpsclrp2_implicit)
      if (allocated(ic%term_wpsclrp2_explicit)) deallocate(ic%term_wpsclrp2_explicit)
      if (allocated(ic%coef_wprtpsclrp_implicit)) deallocate(ic%coef_wprtpsclrp_implicit)
      if (allocated(ic%term_wprtpsclrp_explicit)) deallocate(ic%term_wprtpsclrp_explicit)
      if (allocated(ic%coef_wpthlpsclrp_implicit)) deallocate(ic%coef_wpthlpsclrp_implicit)
      if (allocated(ic%term_wpthlpsclrp_explicit)) deallocate(ic%term_wpthlpsclrp_explicit)

      allocate(ic%coef_wp2sclrp_implicit(ngrdcol, nz, sclr_dim))
      allocate(ic%term_wp2sclrp_explicit(ngrdcol, nz, sclr_dim))
      allocate(ic%coef_wpsclrp2_implicit(ngrdcol, nz, sclr_dim))
      allocate(ic%term_wpsclrp2_explicit(ngrdcol, nz, sclr_dim))
      allocate(ic%coef_wprtpsclrp_implicit(ngrdcol, nz, sclr_dim))
      allocate(ic%term_wprtpsclrp_explicit(ngrdcol, nz, sclr_dim))
      allocate(ic%coef_wpthlpsclrp_implicit(ngrdcol, nz, sclr_dim))
      allocate(ic%term_wpthlpsclrp_explicit(ngrdcol, nz, sclr_dim))
    end if
  end subroutine ensure_implicit_3d_allocated

end module derived_type_storage
