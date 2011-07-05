!----------------------------------------------------------------------
! $Id$
module clex9_nov02

! Description:
!   Contains subroutines for the Clex9 Nov 02 case.

! References:
!   See below.
!----------------------------------------------------------------------

  implicit none

  public :: clex9_nov02_read_t_dependent

  ! Note: bottom point not at surface, so there is no sfc
  ! subroutine

  private ! Default Scope

  contains

!-----------------------------------------------------------------------
  subroutine clex9_nov02_read_t_dependent( time, &
                                           T_sfc, sens_ht, latent_ht )

! Description:
!   This subroutine reads in the values from the _surface.in file for
!   this case.
!
! References:
!   None
!--------------------------------------------------------------------------

    use stats_precision, only: &
      time_precision ! Variable(s)

    use time_dependent_input, only: time_sfc_given, T_sfc_given, & ! Variable(s)
                                    sens_ht_given, latent_ht_given, &
                                    time_select                   ! Procedure(s)

    use interpolation, only: factor_interp ! Procedure(s)

    implicit none

    ! Input variables
    real(kind=time_precision), intent(in) :: &
      time             ! Current time          [s]

    ! Output variables
    real, intent(out) :: &
      T_sfc, &   ! surface temperature [K]
      sens_ht,    &   ! sensible heat flux [W/m^2]
      latent_ht         ! latent heat flux [W/m^2]

    ! Local variables
    real :: &
      time_frac ! time fraction used for interpolation

    integer :: &
      before_time, after_time  ! time indexes used for interpolation


    ! ---- Begin Code ----

    ! interpolate T_sfc from time_dependent_input

    call time_select( time, size(time_sfc_given), time_sfc_given, &
                      before_time, after_time, time_frac )

    T_sfc = factor_interp( time_frac, T_sfc_given(after_time), &
                                      T_sfc_given(before_time) )

    sens_ht = factor_interp( time_frac, sens_ht_given(after_time), &
                                   sens_ht_given(before_time) )

    latent_ht = factor_interp( time_frac, latent_ht_given(after_time), &
                                   latent_ht_given(before_time) )

    return

  end subroutine clex9_nov02_read_t_dependent 


end module clex9_nov02
