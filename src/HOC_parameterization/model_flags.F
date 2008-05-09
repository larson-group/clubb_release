!-----------------------------------------------------------------------
! $Id: model_flags.F,v 1.12 2008-05-09 16:26:32 faschinj Exp $

        module model_flags

!       Description:
!       Various model options that can be toggled off and on as desired.

!       References:
!       None
!-----------------------------------------------------------------------

        implicit none

        public ! All elements of this module are public scope

        logical, parameter :: 
     .  LH_on    = .false., ! Latin hypercube calculation
     .  local_kk = .false.  ! Local drizzle for rain microphysics

        logical :: 
     .  lbugsrad,     ! BUGSrad interactive radiation scheme
     .  kk_rain,      ! Khairoutdinov and Kogan (2000) drizzle scheme. - Brian
     .  licedfs,      ! Simplified ice scheme
     .  lcoamps_micro,! COAMPS rain microphysics
     .  cloud_sed,    ! Cloud water droplet sedimentation. - Brian
     .  luv_nudge,    ! For wind speed nudging. - Michael Falk
     .  lKhm_aniso    ! For anisotropic Khm, as in GABLS2.

!$omp  threadprivate(lbugsrad, kk_rain, licedfs)
!$omp  threadprivate(lcoamps_micro, cloud_sed, luv_nudge)
!$omp  threadprivate(lKhm_aniso)

       end module model_flags

!-----------------------------------------------------------------------
        subroutine setup_model_flags
     .             ( lbugsrad_in, kk_rain_in, cloud_sed_in, 
     .               licedfs_in, lcoamps_micro_in,
     .               luv_nudge_in, lKhm_aniso_in )
!      Description:
!      Setup model flags
!-----------------------------------------------------------------------
        use constants, only: fstderr
        use model_flags

        implicit none

        ! Input Variables
        logical, intent(in) :: 
     .  lbugsrad_in, kk_rain_in, cloud_sed_in,
     .  licedfs_in, lcoamps_micro_in, luv_nudge_in,
     .  lKhm_aniso_in

        lbugsrad      = lbugsrad_in
        kk_rain       = kk_rain_in
        cloud_sed     = cloud_sed_in
        lcoamps_micro = lcoamps_micro_in
        licedfs       = licedfs_in
        luv_nudge     = luv_nudge_in
        lKhm_aniso    = lKhm_aniso_in

        ! Make sure only one microphysical scheme is enabled.
!       if ( .not.( kk_rain .and. lcoamps_micro ) .and. 
!    .       .not.( kk_rain .and. licedfs ) .and. 
!    .       .not.( licedfs .and. lcoamps_micro ) ) then
        if ( count( (/kk_rain, lcoamps_micro, licedfs/) ) > 1 ) then

          write(unit=fstderr, fmt='(3(a18,l1,a1))') 
     .      "kk_rain = ", kk_rain, ",",
     .      "lcoamps_micro = ", lcoamps_micro,",",
     .      "licedfs = ", licedfs, "."
          stop "Only one microphysics scheme may be enabled per run"

        end if ! More than one microphysical scheme enabled

        return
        end subroutine setup_model_flags
