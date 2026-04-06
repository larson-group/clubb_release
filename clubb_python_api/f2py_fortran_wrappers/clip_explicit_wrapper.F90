! clip_explicit_wrapper.F90 — wrappers organized by source module

subroutine f2py_clip_rcm(nzt, ngrdcol, rtm, message, rcm)

  use clubb_precision, only: core_rknd
  use clip_explicit, only: clip_rcm

  implicit none

  integer, intent(in) :: nzt, ngrdcol
  real(core_rknd), dimension(ngrdcol, nzt), intent(in) :: rtm
  character(len=*), intent(in) :: message
  real(core_rknd), dimension(ngrdcol, nzt), intent(inout) :: rcm

  call clip_rcm(nzt, ngrdcol, rtm, message, rcm)

end subroutine f2py_clip_rcm

subroutine f2py_clip_covar(nzm, ngrdcol, solve_type, &
    l_first_clip_ts, l_last_clip_ts, dt, xp2, yp2, &
    l_predict_upwp_vpwp, xpyp, xpyp_chnge)

  use clubb_precision, only: core_rknd
  use derived_type_storage, only: stored_stats
  use clip_explicit, only: clip_covar

  implicit none

  integer, intent(in) :: nzm, ngrdcol, solve_type
  logical, intent(in) :: l_first_clip_ts, l_last_clip_ts
  real(core_rknd), intent(in) :: dt
  real(core_rknd), dimension(ngrdcol, nzm), intent(in) :: xp2, yp2
  logical, intent(in) :: l_predict_upwp_vpwp
  real(core_rknd), dimension(ngrdcol, nzm), intent(inout) :: xpyp
  real(core_rknd), dimension(ngrdcol, nzm), intent(out) :: xpyp_chnge

  call clip_covar(nzm, ngrdcol, solve_type, &
    l_first_clip_ts, l_last_clip_ts, &
    dt, xp2, yp2, l_predict_upwp_vpwp, &
    stored_stats, xpyp, xpyp_chnge)

end subroutine f2py_clip_covar

subroutine f2py_clip_variance(nzm, ngrdcol, solve_type, dt, &
    threshold_lo, xp2)

  use clubb_precision, only: core_rknd
  use derived_type_storage, only: stored_grid, stored_stats
  use clip_explicit, only: clip_variance

  implicit none

  integer, intent(in) :: nzm, ngrdcol, solve_type
  real(core_rknd), intent(in) :: dt
  real(core_rknd), dimension(ngrdcol, nzm), intent(in) :: threshold_lo
  real(core_rknd), dimension(ngrdcol, nzm), intent(inout) :: xp2

  call clip_variance(nzm, ngrdcol, stored_grid, solve_type, dt, &
    threshold_lo, stored_stats, xp2)

end subroutine f2py_clip_variance

subroutine f2py_clip_covars_denom(nzm, ngrdcol, sclr_dim, sclr_dim_transport, dt, &
    rtp2, thlp2, up2, vp2, wp2, sclrp2, &
    wprtp_cl_num, wpthlp_cl_num, wpsclrp_cl_num, &
    upwp_cl_num, vpwp_cl_num, &
    l_predict_upwp_vpwp, l_tke_aniso, l_linearize_pbl_winds, &
    wprtp, wpthlp, upwp, vpwp, wpsclrp, upwp_pert, vpwp_pert)

  use clubb_precision, only: core_rknd
  use derived_type_storage, only: stored_stats
  use clip_explicit, only: clip_covars_denom

  implicit none

  integer, intent(in) :: nzm, ngrdcol, sclr_dim, sclr_dim_transport
  real(core_rknd), intent(in) :: dt
  real(core_rknd), dimension(ngrdcol, nzm), intent(in) :: &
    rtp2, thlp2, up2, vp2, wp2
  real(core_rknd), dimension(ngrdcol, nzm, sclr_dim_transport), intent(in) :: sclrp2
  integer, intent(in) :: wprtp_cl_num, wpthlp_cl_num, wpsclrp_cl_num
  integer, intent(in) :: upwp_cl_num, vpwp_cl_num
  logical, intent(in) :: l_predict_upwp_vpwp, l_tke_aniso
  logical, intent(in) :: l_linearize_pbl_winds
  real(core_rknd), dimension(ngrdcol, nzm), intent(inout) :: &
    wprtp, wpthlp, upwp, vpwp
  real(core_rknd), dimension(ngrdcol, nzm, sclr_dim_transport), intent(inout) :: wpsclrp
  real(core_rknd), dimension(ngrdcol, nzm), intent(inout) :: upwp_pert, vpwp_pert

  call clip_covars_denom(nzm, ngrdcol, sclr_dim, dt, &
    rtp2, thlp2, up2, vp2, wp2, sclrp2, &
    wprtp_cl_num, wpthlp_cl_num, wpsclrp_cl_num, &
    upwp_cl_num, vpwp_cl_num, &
    l_predict_upwp_vpwp, l_tke_aniso, &
    l_linearize_pbl_winds, &
    stored_stats, &
    wprtp, wpthlp, upwp, vpwp, wpsclrp, upwp_pert, vpwp_pert)

end subroutine f2py_clip_covars_denom

subroutine f2py_clip_skewness(nzt, ngrdcol, dt, sfc_elevation, &
    skw_max_mag, wp2_zt, l_use_wp3_lim_with_smth_heaviside, wp3)

  use clubb_precision, only: core_rknd
  use derived_type_storage, only: stored_grid, stored_stats
  use clip_explicit, only: clip_skewness

  implicit none

  integer, intent(in) :: nzt, ngrdcol
  real(core_rknd), intent(in) :: dt
  real(core_rknd), dimension(ngrdcol), intent(in) :: sfc_elevation, skw_max_mag
  real(core_rknd), dimension(ngrdcol, nzt), intent(in) :: wp2_zt
  logical, intent(in) :: l_use_wp3_lim_with_smth_heaviside
  real(core_rknd), dimension(ngrdcol, nzt), intent(inout) :: wp3

  call clip_skewness(nzt, ngrdcol, stored_grid, dt, sfc_elevation, &
    skw_max_mag, wp2_zt, l_use_wp3_lim_with_smth_heaviside, stored_stats, wp3)

end subroutine f2py_clip_skewness
