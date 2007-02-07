!-----------------------------------------------------------------------
! $Id: model_flags.F,v 1.7 2007-02-07 22:36:30 mjfalk Exp $

! Various model options that can be toggled off and on as desired
!-----------------------------------------------------------------------
        module model_flags

        logical, parameter :: 
     .  LH_on    = .false., ! Latin hypercube calculation
     .  local_kk = .false.  ! Local drizzle for rain microphysics

        logical :: 
     .  lbugsrad,     ! BUGSrad interactive radiation scheme
     .  kk_rain,      ! Khairoutdinov and Kogan (2000) drizzle scheme. - Brian
     .  licedfs,      ! Simplified ice scheme
     .  lcoamps_micro,! COAMPS rain microphysics
     .  cloud_sed,    ! Cloud water droplet sedimentation. - Brian
     .  lmpace_nudge, ! For wind speed nudging. - Michael Falk
     .  lKhm_aniso    ! For anisotropic Khm, as in GABLS.

       ! Flag for the lowest model level being at the surface (true)
       ! or not at the surface (false).  Example of case where it is not
       ! at the surface:  Nov. 11 altocumulus case.
       logical :: bottom_at_sfc

!$omp threadprivate(lbugsrad, kk_rain, licedfs)
!$omp threadprivate(lcoamps_micro, cloud_sed, lmpace_nudge)
!$omp threadprivate(bottom_at_sfc)

       end module model_flags

!-----------------------------------------------------------------------
        subroutine setup_model_flags
     .             ( lbugsrad_in, kk_rain_in, cloud_sed_in, 
     .               licedfs_in, lcoamps_micro_in, bottom_at_sfc_in,
     .               lmpace_nudge_in, lKhm_aniso_in )
!      Description:
!      Setup model flags
!-----------------------------------------------------------------------
        use constants, only: fstderr
        use model_flags

        implicit none

        ! Input Variables
        logical, intent(in) :: 
     .  lbugsrad_in, kk_rain_in, cloud_sed_in,
     .  bottom_at_sfc_in, licedfs_in, lcoamps_micro_in, lmpace_nudge_in,
     .  lKhm_aniso_in

        lbugsrad      = lbugsrad_in
        kk_rain       = kk_rain_in
        cloud_sed     = cloud_sed_in
        bottom_at_sfc = bottom_at_sfc_in
        lcoamps_micro = lcoamps_micro_in
        licedfs       = licedfs_in
        lmpace_nudge  = lmpace_nudge_in
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
