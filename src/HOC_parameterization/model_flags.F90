!-----------------------------------------------------------------------
! $Id: model_flags.F90,v 1.3 2008-07-30 21:22:27 faschinj Exp $

module model_flags

! Description:
!   Various model options that can be toggled off and on as desired.

! References:
!   None
!-----------------------------------------------------------------------

  implicit none

  public :: setup_model_flags
  private ! Default Scope

  logical, parameter, public ::  & 
    LH_on      = .false., & ! Latin hypercube calculation
    local_kk   = .false., & ! Local drizzle for rain microphysics
    lpos_def   = .true.,  & ! Flux limiting pos. def. scheme on rtm
    lhole_fill = .true.     ! Hole filling pos. def. scheme on wp2,up2,rtp2,etc

  logical, parameter, public ::  & 
    lsingle_C2_Skw = .false.,  & ! Use a single Skw dependent value for C2
    lgamma_Skw     = .true.,   & ! Use a Skw dependent gamma parameter
    lbyteswap_io   = .false.     ! Swap byte order in GrADS output

  logical, public ::  & 
    l_bugsrad,      & ! BUGSrad interactive radiation scheme
    l_kk_rain,       & ! Khairoutdinov and Kogan (2000) drizzle scheme. - Brian
    l_licedfs,       & ! Simplified ice scheme
    l_coamps_micro, & ! COAMPS rain microphysics
    l_cloud_sed,     & ! Cloud water droplet sedimentation. - Brian
    l_uv_nudge,     & ! For wind speed nudging. - Michael Falk
    l_Khm_aniso    ! For anisotropic Khm, as in GABLS2.

!$omp threadprivate(l_bugsrad, l_kk_rain, l_licedfs)
!$omp threadprivate(l_coamps_micro, l_cloud_sed, l_uv_nudge)
!$omp threadprivate(l_Khm_aniso)

  contains
!-----------------------------------------------------------------------
  subroutine setup_model_flags & 
             ( l_bugsrad_in, l_kk_rain_in, l_cloud_sed_in,  & 
               l_licedfs_in, l_coamps_micro_in, & 
               l_uv_nudge_in, l_Khm_aniso_in )

! Description:
!   Setup model flags
!-----------------------------------------------------------------------
    use constants, only:  & 
      fstderr ! Variable(s)

    implicit none

    ! Input Variables
    logical, intent(in) ::  & 
      l_bugsrad_in, l_kk_rain_in, l_cloud_sed_in, & 
      l_licedfs_in, l_coamps_micro_in, l_uv_nudge_in, & 
      l_Khm_aniso_in

!-----------------------------------------------------------------------
    l_bugsrad      = l_bugsrad_in
    l_kk_rain       = l_kk_rain_in
    l_cloud_sed     = l_cloud_sed_in
    l_coamps_micro = l_coamps_micro_in
    l_licedfs       = l_licedfs_in
    l_uv_nudge     = l_uv_nudge_in
    l_Khm_aniso    = l_Khm_aniso_in

        ! Make sure only one microphysical scheme is enabled.
!       if ( .not.( l_kk_rain .and. l_coamps_micro ) .and. &
!            .not.( l_kk_rain .and. l_licedfs ) .and. &
!            .not.( l_licedfs .and. l_coamps_micro ) ) then

    if ( count( (/l_kk_rain, l_coamps_micro, l_licedfs/) ) > 1 ) then

      write(unit=fstderr, fmt='(3(a18,l1,a1))')  & 
        "l_kk_rain = ", l_kk_rain, ",", & 
        "l_coamps_micro = ", l_coamps_micro,",", & 
        "l_licedfs = ", l_licedfs, "."
      stop "Only one microphysics scheme may be enabled per run"

    end if ! More than one microphysical scheme enabled

    return
  end subroutine setup_model_flags

end module model_flags
