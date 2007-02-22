!-----------------------------------------------------------------------
        subroutine compute_Khm_anisotropic
     .             ( C4, C13, c_K, time_current,
     .               time_initial, wtol, lmin, 
     .               um, vm, Lscale,
     .               wp2zt, up2_plus_vp2,
     .               Khm_anisotropic )
!      Description:
!      Recompute Kh (Eddy diffusivity) anisotropically.
!      
!      Author:
!      Michael Falk 11 February 2007.
!-----------------------------------------------------------------------

        use grid_class

        implicit none

        ! External
        intrinsic :: max, sqrt

        ! Input variables
        real, intent(in) :: 
     .  C4,  ! Constant C4 from model namelist                          [-]
     .  C13, ! Constant parameter                                       [-]
     .  c_K  ! Eddy diffusivity coefficient(Duynkerke & Driedonks 1987) [-]

        real, intent(in) :: 
     .  time_current,! Time since zero  [s]
     .  time_initial,! Initial time     [s]   eMFc
     .  wtol,        ! w tolerance      [m/s]
     .  lmin         ! lmin             [m]

        real, dimension(gr%nnzp), intent(in) ::
     .  um,     ! u-direction wind speed        [m s^-1]
     .  vm,     ! v-direction wind speed        [m s^-1]
     .  Lscale, ! Turbulent length scale        [m]
     .  wp2zt   ! wp2 interpolated to t grid    [m^2 s^-2] eMFc

        ! Output variables
        real, dimension(gr%nnzp), intent(out) ::
     .  up2_plus_vp2,    ! up2_plus_vp2                         [m^2 s^-2]
     .  Khm_anisotropic  ! new Khm based on new e_anisotropic   [m^2 s^-1]

        ! Local variables
        real, dimension(gr%nnzp) ::
     .  du_dz,  ! Derivative of um with respect to height       [s^-1]
     .  dv_dz   ! Derivative of vm with respect to height       [s^-1]

        ! The new TKE, computed without the assumption that up2=vp2=wp2
        real, dimension(gr%nnzp) ::
     .  e_anisotropic    ! TKE is em in the statistics.         [m^2/s^2]

        real, dimension(gr%nnzp) ::
     .  Lscale_with_min, ! Turbulent length scale               [m]
     .  wp2zt_with_min   ! wp2 interpolated to t grid           [m^2 s^-2]

        integer :: k ! Index

        du_dz = ddzm( zt2zm( um ) )
        dv_dz = ddzm( zt2zm( vm ) )

        do k=1,gr%nnzp
          wp2zt_with_min(k)  = max( wp2zt(k), wtol**2 )
          Lscale_with_min(k) = max( Lscale(k), lmin )
! These are the old values.  7 Feb 2007, Michael Falk
!          wp2zt_with_min(k) = max(wp2zt(k), 0.3)
!          Lscale_with_min(k) = max(Lscale(k), 30.0)
          up2_plus_vp2(k) = 1. / (C13 + (1./3.) * C4)
     .                      * ( 2 * c_K * Lscale_with_min(k)
     .                          * Lscale_with_min(k)
     .                          *   (  (du_dz(k) * du_dz(k))
     .                               + (dv_dz(k) * dv_dz(k)) )
     .                          + (2./3.) * C4 * wp2zt_with_min(k)
     .                        )
          e_anisotropic(k) = (0.5 * up2_plus_vp2(k))
     .                     + (0.5 * wp2zt_with_min(k))

          Khm_anisotropic(k) = c_K * Lscale_with_min(k) *
     .                         sqrt( e_anisotropic(k) )
        end do

        return
        end subroutine compute_Khm_anisotropic
