!-----------------------------------------------------------------------
! $Id: explicit_clip.F90,v 1.8 2008-08-07 16:10:13 griffinb Exp $
!===============================================================================
module explicit_clip

implicit none

private

public :: covariance_clip, & 
          variance_clip, & 
          skewness_clip

contains

!===============================================================================
subroutine covariance_clip( solve_type, l_first_clip_ts,  & 
                            l_last_clip_ts, dt, xp2, yp2,  & 
                            xpyp )

! Description:
! Clipping the value of covariance x'y' based on the correlation between x 
! and y.
!
! The correlation between variables x and y is:
!
! corr_(x,y) = x'y' / [ sqrt(x'^2) * sqrt(y'^2) ];
!
! where x'^2 is the variance of x, y'^2 is the variance of y, and x'y' is the 
! covariance of x and y.
!
! The correlation of two variables must always have a value between -1 and 1, 
! such that:
!
! -1 <= corr_(x,y) <= 1.
!
! Therefore, there is an upper limit on x'y', such that:
!
! x'y' <=  [ sqrt(x'^2) * sqrt(y'^2) ];
!
! and a lower limit on x'y', such that:
!
! x'y' >= -[ sqrt(x'^2) * sqrt(y'^2) ].
!
! The values of x'y', x'^2, and y'^2 are all found on momentum levels.
!
! The value of x'y' may need to be clipped whenever x'y', x'^2, or y'^2 is 
! updated.
!
! The following covariances are found in the code:  
!
! w'r_t', w'th_l', w'sclr', (computed in mixing.F);
! r_t'th_l', sclr'r_t', sclr'th_l', (computed in diag_var.F);
! u'w', v'w', w'edsclr' (computed in compute_um_edsclrm_mod.F).

! References:
!-----------------------------------------------------------------------

use grid_class, only: & 
    gr ! Variable(s)

use stats_precision, only: & 
    time_precision ! Variable(s)
 
use stats_type, only: & 
    stat_begin_update,  & ! Procedure(s)
    stat_modify, & 
    stat_end_update

use stats_variables, only: & 
    zm,  & ! Variable(s)
    iwprtp_cl, & 
    iwpthlp_cl, & 
    irtpthlp_cl, & 
    l_stats_samp

implicit none

! Input Variables
character(len=*), intent(in) :: & 
  solve_type       ! Variable being solved; used for STATS.

logical, intent(in) :: & 
  l_first_clip_ts, & ! First instance of clipping in a timestep.
  l_last_clip_ts     ! Last instance of clipping in a timestep.

real(kind=time_precision), intent(in) ::  & 
  dt     ! Model timestep; used here for STATS           [s]

real, dimension(gr%nnzp), intent(in) :: & 
  xp2, & ! Variance of x, x'^2 (momentum levels)         [{x units}^2]
  yp2    ! Variance of y, y'^2 (momentum levels)         [{y units}^2]

! Output Variable
real, dimension(gr%nnzp), intent(inout) :: & 
  xpyp   ! Covariance of x and y, x'y' (momentum levels) [{x units}*{y units}]

 
! Local Variable
integer :: & 
  ixpyp_cl


select case ( trim( solve_type ) )
case ( "wprtp" )   ! wprtp clipping budget term
   ixpyp_cl = iwprtp_cl
case ( "wpthlp" )   ! wpthlp clipping budget term
   ixpyp_cl = iwpthlp_cl
case ( "rtpthlp" )   ! rtpthlp clipping budget term
   ixpyp_cl = irtpthlp_cl
case default   ! scalars (or upwp/vpwp) are involved
   ixpyp_cl = 0
end select
 

if ( l_stats_samp ) then
   if ( l_first_clip_ts ) then
      call stat_begin_update( ixpyp_cl, real( xpyp / dt ), zm )
   else
      call stat_modify( ixpyp_cl, real( -xpyp / dt ), zm )
   endif
endif 

! Clipping for xpyp at an upper limit corresponding with 
! a correlation between x and y of 0.99.
where ( xpyp >  0.99 * sqrt( xp2 * yp2 ) ) & 
   xpyp =  0.99 * sqrt( xp2 * yp2 )

! Clipping for xpyp at a lower limit corresponding with 
! a correlation between x and y of -0.99.
where ( xpyp < -0.99 * sqrt( xp2 * yp2 ) ) & 
   xpyp = -0.99 * sqrt( xp2 * yp2 )

if ( l_stats_samp ) then
   if ( l_last_clip_ts ) then
      call stat_end_update( ixpyp_cl, real( xpyp / dt ), zm )
   else
      call stat_modify( ixpyp_cl, real( xpyp / dt ), zm )
   endif
endif


end subroutine covariance_clip

!===============================================================================
subroutine variance_clip( solve_type, dt, threshold, xp2 )

! Description:
! Clipping the value of variance x'^2 based on a minimum threshold value.  The 
! threshold value must be greater than or equal to 0.
!
! The values of x'^2 are found on the momentum levels.
!
! The following variances are found in the code:  
!
! r_t'^2, th_l'^2, u'^2, v'^2, sclr'^2, (computed in diag_var.F);
! w'^2 (computed in wp23.F).

! References:
!-----------------------------------------------------------------------

use grid_class, only: & 
    gr ! Variable(s)

use stats_precision, only: & 
    time_precision ! Variable(s)
 
use stats_type, only: & 
    stat_begin_update,  & ! Procedure(s)
    stat_end_update

use stats_variables, only: & 
    zm,  & ! Variable(s)
    iwp2_cl, & 
    irtp2_cl, & 
    ithlp2_cl, & 
    iup2_cl, & 
    ivp2_cl, & 
    l_stats_samp
 

implicit none

! Input Variables
character(len=*), intent(in) :: & 
  solve_type  ! Variable being solved; used for STATS.

real(kind=time_precision), intent(in) :: & 
  dt          ! Model timestep; used here for STATS     [s]

real, intent(in) :: & 
  threshold   ! Minimum value of x'^2                   [{x units}^2]

! Output Variable
real, dimension(gr%nnzp), intent(inout) :: & 
  xp2         ! Variance of x, x'^2 (momentum levels)   [{x units}^2]

! Local Variables
integer :: k   ! Array index

 
integer :: & 
  ixp2_cl


select case ( trim( solve_type ) )
case ( "wp2" )   ! wp2 clipping budget term
   ixp2_cl = iwp2_cl
case ( "rtp2" )   ! rtp2 clipping budget term
   ixp2_cl = irtp2_cl
case ( "thlp2" )   ! thlp2 clipping budget term
   ixp2_cl = ithlp2_cl
case ( "up2" )   ! up2 clipping budget term
   ixp2_cl = iup2_cl
case ( "vp2" )   ! vp2 clipping budget term
   ixp2_cl = ivp2_cl
case default   ! scalars are involved
   ixp2_cl = 0
end select


if ( l_stats_samp ) then 
   call stat_begin_update( ixp2_cl, real( xp2 / dt ), zm )
endif
 
! Limit the value of x'^2 at threshold.
do k = 2, gr%nnzp, 1
   if ( xp2(k) < threshold ) then
      xp2(k) = threshold
   endif
enddo

if ( l_stats_samp ) then
   call stat_end_update( ixp2_cl, real( xp2 / dt ), zm )
endif


end subroutine variance_clip

!===============================================================================
subroutine skewness_clip( dt, wp2_zt, wp3 )

! Description:
! Clipping the value of w'^3 based on the skewness of w, Sk_w.
!
! The skewness of w is:
!
! Sk_w = w'^3 / (w'^2)^(3/2).
!
! The value of Sk_w is limited to a range between an upper limit and a lower 
! limit.  The values of the limits depend on whether the level altitude is 
! within 100 meters of the surface.
!
! For altitudes within 100 meters of the surface:
!
! -0.2*sqrt(2) <= Sk_w <= 0.2*sqrt(2);
!
! while for all other altitudes:
!
! -4.5 <= Sk_w <= 4.5.
!
! Therefore, there is an upper limit on w'^3, such that:
!
! w'^3  <=  threshold_magnitude * (w'^2)^(3/2);
!
! and a lower limit on w'^3, such that:
!
! w'^3  >= -threshold_magnitude * (w'^2)^(3/2).
!
! The values of w'^3 are found on the thermodynamic levels, while the values of
! w'^2 are found on the momentum levels.  Therefore, the values of w'^2 are
! interpolated to the thermodynamic levels before being used to calculate the
! upper and lower limits for w'^3.

! References:
!-----------------------------------------------------------------------

use grid_class, only: & 
    gr ! Variable(s)

use stats_precision, only: & 
    time_precision ! Variable(s)

use stats_type, only: &
    stat_begin_update,  & ! Procedure(s)
    stat_end_update

use stats_variables, only: & 
    zt,  & ! Variable(s)
    iwp3_cl, & 
    l_stats_samp

implicit none

! Input Variables
real(kind=time_precision), intent(in) :: & 
  dt               ! Model timestep; used here for STATS        [s]

real, dimension(gr%nnzp), intent(in) :: &
  wp2_zt           ! w'^2 interpolated to thermodyamic levels   [m^2/s^2]

! Input/Output Variables
real, dimension(gr%nnzp), intent(inout) :: &
  wp3              ! w'^3 (thermodynamic levels)                [m^3/s^3]

! Local Variables
real, dimension(gr%nnzp) :: &
  wp3_upper_lim, & ! Keeps Sk_w from becoming > upper_limit     [m^3/s^3]
  wp3_lower_lim    ! Keeps Sk_w from becoming < lower_limit     [m^3/s^3]

integer :: k       ! Vertical array index.


if ( l_stats_samp ) then 
   call stat_begin_update( iwp3_cl, real( wp3 / dt ), zt )
endif

! Compute the upper and lower limits of w'^3 at every level,
! based on the skewness of w, Sk_w, such that:
! Sk_w = w'^3 / (w'^2)^(3/2);
! -4.5 <= Sk_w <= 4.5; 
! or, if the level altitude is within 100 meters of the surface,
! -0.2*sqrt(2) <= Sk_w <= 0.2*sqrt(2).

! The normal magnitude limit of skewness of w in the CLUBB code is 4.5.
! However, according to Andre et al. (1976b & 1978), wp3 should not exceed 
! [2*(wp2^3)]^(1/2) at any level.  However, this term should be multiplied 
! by 0.2 close to the surface to include surface effects.  There already is 
! a wp3 clipping term in place for all other altitudes, but this term will be
! included for the surface layer only.  Therefore, the lowest level wp3 should 
! not exceed 0.2 * sqrt(2) * wp2^(3/2).  Brian Griffin.  12/18/05.

do k = 1, gr%nnzp, 1
   if ( gr%zt(k) <= 100.0 ) then ! Clip for 100 m. above ground.
      wp3_upper_lim(k) =  0.2 * sqrt(2.0) * wp2_zt(k)**(3.0/2.0)
      wp3_lower_lim(k) = -0.2 * sqrt(2.0) * wp2_zt(k)**(3.0/2.0)
   else                          ! Clip skewness consistently with a.
      wp3_upper_lim(k) =  4.5 * wp2_zt(k)**(3.0/2.0)
      wp3_lower_lim(k) = -4.5 * wp2_zt(k)**(3.0/2.0)
   endif
enddo

! Clipping for w'^3 at an upper limit corresponding with 
! the appropriate value of Sk_w.
where ( wp3 > wp3_upper_lim ) &
   wp3 = wp3_upper_lim

! Clipping for w'^3 at a lower limit corresponding with 
! the appropriate value of Sk_w.
where ( wp3 < wp3_lower_lim ) &
   wp3 = wp3_lower_lim

if ( l_stats_samp ) then
   call stat_end_update( iwp3_cl, real( wp3 / dt ), zt )
endif


end subroutine skewness_clip

!===============================================================================

end module explicit_clip
