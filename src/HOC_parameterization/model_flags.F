!-----------------------------------------------------------------------
! $Id: model_flags.F,v 1.3 2006-11-01 20:42:56 mjfalk Exp $

! Various model options that can be toggled off and on as desired
!-----------------------------------------------------------------------
       module model_flags

       logical, parameter :: LH_on    = .false. ! latin hypercube calculation
       logical, parameter :: local_kk = .false. ! local drizzle for rain microphysics
       logical, parameter :: lnudging = .false. ! for wind speed nudging.  Michael Falk, 1 Nov 2006

!      BUGSrad interactive radiation scheme
       logical :: lbugsrad

!      Khairoutdinov and Kogan (2000) drizzle scheme. - Brian
       logical :: kk_rain

!      Cloud water droplet sedimentation. - Brian
       logical :: cloud_sed

!      Flag for the lowest model level being at the surface (true)
!      or not at the surface (false).  Example of case where it is not
!      at the surface:  Nov. 11 altocumulus case.
       logical :: bottom_at_sfc

       end module model_flags
!-----------------------------------------------------------------------

       subroutine setup_model_flags( lbugsrad_in, kk_rain_in, 
     .                               cloud_sed_in, bottom_at_sfc_in )
       use model_flags
       implicit none

       logical, intent(in) :: lbugsrad_in, kk_rain_in, cloud_sed_in
       logical, intent(in) :: bottom_at_sfc_in

       lbugsrad      = lbugsrad_in
       kk_rain       = kk_rain_in
       cloud_sed     = cloud_sed_in
       bottom_at_sfc = bottom_at_sfc_in

       return
       end subroutine setup_model_flags

