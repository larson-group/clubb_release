!-----------------------------------------------------------------------
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
    corr_wrr_NL_cloud,  & ! Prescribed in-cloud correlation of w and rr ip   [-]
    corr_wNr_NL_cloud,  & ! Prescribed in-cloud correlation of w and Nr ip   [-]
    corr_wNcn_NL_cloud, & ! Prescribed in-cloud correlation of w and Ncn     [-]
    corr_sw_NN_cloud,   & ! Prescribed in-cloud correlation of s and w       [-]
    corr_srr_NL_cloud,  & ! Prescribed in-cloud correlation of s and rr ip   [-]
    corr_sNr_NL_cloud,  & ! Prescribed in-cloud correlation of s and Nr ip   [-]
    corr_sNcn_NL_cloud, & ! Prescribed in-cloud correlation of s and Ncn     [-]
    corr_trr_NL_cloud,  & ! Prescribed in-cloud correlation of t and rr ip   [-]
    corr_tNr_NL_cloud,  & ! Prescribed in-cloud correlation of t and Nr ip   [-]
    corr_tNcn_NL_cloud, & ! Prescribed in-cloud correlation of t and Ncn     [-]
    corr_rrNr_LL_cloud    ! Prescribed in-cloud correlation of rr and Nr ip  [-]

!$omp threadprivate( corr_wrr_NL_cloud, corr_wNr_NL_cloud, corr_wNcn_NL_cloud, &
!$omp                corr_srr_NL_cloud, corr_sNr_NL_cloud, corr_sNcn_NL_cloud, &
!$omp                corr_trr_NL_cloud, corr_tNr_NL_cloud, corr_tNcn_NL_cloud, &
!$omp                corr_rrNr_LL_cloud, corr_sw_NN_cloud )

  ! Prescribed correlations (read in from a correlation array) for grid levels
  ! without any cloud.
  real( kind = core_rknd ), public :: &  ! RF02 value
    corr_wrr_NL_below,  & ! Prescribed below-cloud correlation: w and rr ip  [-]
    corr_wNr_NL_below,  & ! Prescribed below-cloud correlation: w and Nr ip  [-]
    corr_wNcn_NL_below, & ! Prescribed below-cloud correlation: w and Ncn    [-]
    corr_sw_NN_below,   & ! Prescribed below-cloud correlation: s and w      [-]
    corr_srr_NL_below,  & ! Prescribed below-cloud correlation: s and rr ip  [-]
    corr_sNr_NL_below,  & ! Prescribed below-cloud correlation: s and Nr ip  [-]
    corr_sNcn_NL_below, & ! Prescribed below-cloud correlation: s and Ncn    [-]
    corr_trr_NL_below,  & ! Prescribed below-cloud correlation: t and rr ip  [-]
    corr_tNr_NL_below,  & ! Prescribed below-cloud correlation: t and Nr ip  [-]
    corr_tNcn_NL_below, & ! Prescribed below-cloud correlation: t and Ncn    [-]
    corr_rrNr_LL_below    ! Prescribed below-cloud correlation: rr and Nr ip [-]

!$omp threadprivate( corr_wrr_NL_below, corr_wNr_NL_below, corr_wNcn_NL_below, &
!$omp                corr_srr_NL_below, corr_sNr_NL_below, corr_sNcn_NL_below, &
!$omp                corr_trr_NL_below, corr_tNr_NL_below, corr_tNcn_NL_below, &
!$omp                corr_rrNr_LL_below, corr_sw_NN_below )

  ! Only needed when l_fix_s_t_correlations is true
  real( kind = core_rknd ), public :: &
    corr_st_NN_cloud, & ! Prescribed in-cloud correlation of s and t         [-]
    corr_st_NN_below    ! Prescribed below-cloud correlation of s and t      [-]

!$omp threadprivate( corr_st_NN_cloud, corr_st_NN_below )

  integer, public :: &
    iipcorr_w,   & ! Row/Column index of w in prescribed correlation array
    iipcorr_s,   & ! Row/Column index of s in prescribed correlation array
    iipcorr_t,   & ! Row/Column index of t in prescribed correlation array
    iipcorr_Ncn, & ! Row/Column index of N_cn in prescribed correlation array
    iipcorr_rr,  & ! Row/Column index of r_r in prescribed correlation array
    iipcorr_Nr,  & ! Row/Column index of N_r in prescribed correlation array
    iipcorr_ri,  & ! Row/Column index of r_i in prescribed correlation array
    iipcorr_Ni,  & ! Row/Column index of N_i in prescribed correlation array
    iipcorr_rs,  & ! Row/Column index of r_s in prescribed correlation array
    iipcorr_Ns,  & ! Row/Column index of N_s in prescribed correlation array
    iipcorr_rg,  & ! Row/Column index of r_g in prescribed correlation array
    iipcorr_Ng     ! Row/Column index of N_g in prescribed correlation array

!$omp threadprivate( iipcorr_w, iipcorr_s, iipcorr_t, iipcorr_Ncn, &
!$omp                iipcorr_rr, iipcorr_Nr, iipcorr_ri, iipcorr_Ni, &
!$omp                iipcorr_rs, iipcorr_Ns, iipcorr_rg, iipcorr_Ng )

end module fixed_correlations
