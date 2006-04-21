!-----------------------------------------------------------------------
! $Id: model_flags.F,v 1.1 2006-04-21 21:50:56 dschanen Exp $

! Various model options that can be toggled off and on as desired
!-----------------------------------------------------------------------
       module model_flags

       logical, parameter :: LH_on    = .false. ! latin hypercube calculation
       logical, parameter :: local_kk = .false. ! local drizzle for rain microphysics

!      BUGSrad interactive radiation scheme
       logical :: lbugsrad

!      Cloud water droplet sedimentation. - Brian
       logical :: kk_rain

!      Cloud water droplet sedimentation. - Brian
       logical :: cloud_sed

       end module model_flags
!-----------------------------------------------------------------------

       subroutine setup_model_flags( lbugsrad_in, kk_rain_in, 
     .                               cloud_sed_in )
       use model_flags
       implicit none

       logical, intent(in) :: lbugsrad_in, kk_rain_in, cloud_sed_in

       lbugsrad  = lbugsrad_in
       kk_rain   = kk_rain_in
       cloud_sed = cloud_sed_in

       return
       end subroutine setup_model_flags

