!$Id$
!===============================================================================
module fixed_correlations

  use clubb_precision, only: &
      core_rknd

  implicit none

  private ! Default scope

  ! Prescribed correlations (read in from a correlation array) for grid levels
  ! with at least some cloud.
  real( kind = core_rknd ), public :: &
    corr_wrr_NL_cloud,  & ! Prescribed in-cloud correlation of w and rr     [-]
    corr_wNr_NL_cloud,  & ! Prescribed in-cloud correlation of w and Nr     [-]
    corr_wNcn_NL_cloud, & ! Prescribed in-cloud correlation of w and Ncn    [-]
    corr_sw_NN_cloud,   & ! Prescribed in-cloud correlation of s and w      [-]
    corr_srr_NL_cloud,  & ! Prescribed in-cloud correlation of s and rr     [-]
    corr_sNr_NL_cloud,  & ! Prescribed in-cloud correlation of s and Nr     [-]
    corr_sNcn_NL_cloud, & ! Prescribed in-cloud correlation of s and Ncn    [-]
    corr_trr_NL_cloud,  & ! Prescribed in-cloud correlation of t and rr     [-]
    corr_tNr_NL_cloud,  & ! Prescribed in-cloud correlation of t and Nr     [-]
    corr_tNcn_NL_cloud, & ! Prescribed in-cloud correlation of t and Ncn    [-]
    corr_rrNr_LL_cloud    ! Prescribed in-cloud correlation of rr and Nr    [-]

!$omp threadprivate( corr_wrr_NL_cloud, corr_wNr_NL_cloud, corr_wNcn_NL_cloud, &
!$omp                corr_srr_NL_cloud, corr_sNr_NL_cloud, corr_sNcn_NL_cloud, &
!$omp                corr_trr_NL_cloud, corr_tNr_NL_cloud, corr_tNcn_NL_cloud, &
!$omp                corr_rrNr_LL_cloud, corr_sw_NN_cloud )

  ! Prescribed correlations (read in from a correlation array) for grid levels
  ! without any cloud.
  real( kind = core_rknd ), public :: &  ! RF02 value
    corr_wrr_NL_below,  & ! Prescribed below-cloud correlation of w and rr  [-]
    corr_wNr_NL_below,  & ! Prescribed below-cloud correlation of w and Nr  [-]
    corr_wNcn_NL_below, & ! Prescribed below-cloud correlation of w and Ncn [-]
    corr_sw_NN_below,   & ! Prescribed below-cloud correlation of s and w   [-]
    corr_srr_NL_below,  & ! Prescribed below-cloud correlation of s and rr  [-]
    corr_sNr_NL_below,  & ! Prescribed below-cloud correlation of s and Nr  [-]
    corr_sNcn_NL_below, & ! Prescribed below-cloud correlation of s and Ncn [-]
    corr_trr_NL_below,  & ! Prescribed below-cloud correlation of t and rr  [-]
    corr_tNr_NL_below,  & ! Prescribed below-cloud correlation of t and Nr  [-]
    corr_tNcn_NL_below, & ! Prescribed below-cloud correlation of t and Ncn [-]
    corr_rrNr_LL_below    ! Prescribed below-cloud correlation of rr and Nr [-]

!$omp threadprivate( corr_wrr_NL_below, corr_wNr_NL_below, corr_wNcn_NL_below, &
!$omp                corr_srr_NL_below, corr_sNr_NL_below, corr_sNcn_NL_below, &
!$omp                corr_trr_NL_below, corr_tNr_NL_below, corr_tNcn_NL_below, &
!$omp                corr_rrNr_LL_below, corr_sw_NN_below )

  ! Only needed when l_fix_s_t_correlations is true
  real( kind = core_rknd ), public :: &
    corr_st_NN_cloud, & ! Prescribed in-cloud correlation of s and t        [-]
    corr_st_NN_below    ! Prescribed below-cloud correlation of s and t     [-]

!$omp threadprivate( corr_st_NN_cloud, corr_st_NN_below )

end module fixed_correlations
