! $Id: compute_Khm_anisotropic.F,v 1.1 2007-02-09 00:52:20 dschanen Exp $
!-----------------------------------------------------------------------
        subroutine compute_Khm_anisotropic
     .             ( C4, C_k, time, time_initial, um, vm, Lscale, 
     .               wp2zt, up2_plus_vp2, Khm_anisotropic )
!       Description:
!
!       References:
!
!-----------------------------------------------------------------------
        use grid_class

        implicit none

        ! Local constant
        real, parameter ::
     .  C13 = 0.1     ! Constant for up2_plus_vp2

        ! Input variables
        real, intent(in) :: 
     .  C4,          ! Constant C4 from model standalone file   [-]   start Michael Falk, 11 Jan 2007
     .  C_k,         ! Coefficient of something                 [-]
     .  time,        ! time since zero                          [s]
     .  time_initial ! initialization time                      [s]   eMFc

        real, dimension(gr%nnzp), intent(in) ::
     .  um,     ! u-direction wind speed                    [m s^-1]  start Michael Falk, 11 Jan 2007
     .  vm,     ! v-direction wind speed                    [m s^-1]
     .  Lscale, ! Turbulent length scale                    [m]
     .  wp2zt   ! wp2 interpolated to t grid                [m^2 s^-2] eMFc

        ! Output variables
        real, dimension(gr%nnzp), intent(out) ::
     .  up2_plus_vp2,    ! up2_plus_vp2                   [m^2 s^-2]
     .  Khm_anisotropic  ! new Khm based on new e_anisotropic

        ! Local variables
        real, dimension(gr%nnzp) ::
     .  du_dz,           ! Derivative of um with respect to height [s^-1]
     .  dv_dz,           ! Derivative of vm with respect to height [s^-1]
     .  e_anisotropic,   ! new TKE, computed without the assumption that up2=vp2=wp2
     .  Lscale_with_min, ! Turbulent length scale                    [m]
     .  wp2zt_with_min   ! wp2 interpolated to t grid                [m^2 s^-2]

        integer :: k ! Index

        du_dz = ddzm( zt2zm( um ) )
        dv_dz = ddzm( zt2zm( vm ) )

        do k=1,gr%nnzp
          wp2zt_with_min(k) = max(wp2zt(k), 2.e-2*2.e-2) !wtol^2
          Lscale_with_min(k) = max(Lscale(k), 5.0)
! These are the old values.  7 Feb 2007, Michael Falk
!          wp2zt_with_min(k) = max(wp2zt(k), 0.3)
!          Lscale_with_min(k) = max(Lscale(k), 30.0)
          up2_plus_vp2(k) = 1 / (C13 + (1./3.) * C4)
     .                      * ( 2 * C_k * Lscale_with_min(k)
     .                          * Lscale_with_min(k)
     .                          *   ((du_dz(k) * du_dz(k))
     .                              + (dv_dz(k) * dv_dz(k)))
     .                          + (2./3.) * C4 * wp2zt_with_min(k)
     .                        )
          e_anisotropic(k) = (0.5 * up2_plus_vp2(k))
     .                       + (0.5 * wp2zt_with_min(k))

          Khm_anisotropic(k) = C_k * Lscale_with_min(k) *
     .                         sqrt(e_anisotropic(k))
        end do

        return
        end subroutine compute_Khm_anisotropic
