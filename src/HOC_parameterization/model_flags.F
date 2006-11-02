!-----------------------------------------------------------------------
! $Id: model_flags.F,v 1.4 2006-11-02 01:44:41 dschanen Exp $

! Various model options that can be toggled off and on as desired
!-----------------------------------------------------------------------
       module model_flags

       logical, parameter :: 
     . LH_on    = .false., ! latin hypercube calculation
     . local_kk = .false.  ! local drizzle for rain microphysics

       logical ::
     . lbugsrad,    ! BUGSrad interactive radiation scheme
     . kk_rain,     ! Khairoutdinov and Kogan (2000) drizzle scheme. - Brian
     . cloud_sed,   ! Cloud water droplet sedimentation. - Brian
     . lmpace_nudge ! For wind speed nudging. - Michael Falk

       ! Flag for the lowest model level being at the surface (true)
       ! or not at the surface (false).  Example of case where it is not
       ! at the surface:  Nov. 11 altocumulus case.
       logical :: bottom_at_sfc

       end module model_flags

!-----------------------------------------------------------------------
       subroutine setup_model_flags
     .            ( lbugsrad_in, kk_rain_in, 
     .              cloud_sed_in, bottom_at_sfc_in, lmpace_nudge_in)
!      Description:
!      Setup model flags

!      References:
!      None
!-----------------------------------------------------------------------
       use model_flags
       implicit none

       logical, intent(in) :: 
     . lbugsrad_in, 
     . kk_rain_in,
     . cloud_sed_in,
     . bottom_at_sfc_in,
     . lmpace_nudge_in

       lbugsrad      = lbugsrad_in
       kk_rain       = kk_rain_in
       cloud_sed     = cloud_sed_in
       bottom_at_sfc = bottom_at_sfc_in
       lmpace_nudge  = lmpace_nudge_in

       return
       end subroutine setup_model_flags

