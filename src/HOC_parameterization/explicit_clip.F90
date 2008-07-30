!-----------------------------------------------------------------------
! $Id: explicit_clip.F90,v 1.4 2008-07-30 19:17:35 dschanen Exp $
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

!       Description:
!       Clipping the value of covariance x'y' based on the correlation 
!       between x and y.
!
!       The correlation between variables x and y is:
!
!       corr_(x,y) = x'y' / [ sqrt(x'^2) * sqrt(y'^2) ];
!
!       where x'^2 is the variance of x, y'^2 is the variance of y, and
!       x'y' is the covariance of x and y.
!
!       The correlation of two variables must always have a value 
!       between -1 and 1, such that:
!
!       -1 <= corr_(x,y) <= 1.
!
!       Therefore, there is an upper limit on x'y', such that:
!
!       x'y' <=  [ sqrt(x'^2) * sqrt(y'^2) ];
!
!       and a lower limit on x'y', such that:
!
!       x'y' >= -[ sqrt(x'^2) * sqrt(y'^2) ].
!
!       The values of x'y', x'^2, and y'^2 are all found on momentum 
!       levels.
!
!       The value of x'y' may need to be clipped whenever x'y', x'^2, 
!       or y'^2 is updated.
!
!       The following covariances are found in the code:  
!
!       w'r_t', w'th_l', w'sclr', (computed in mixing.F);
!       r_t'th_l', sclr'r_t', sclr'th_l', (computed in diag_var.F);
!       u'w', v'w', w'edsclr' (computed in compute_um_edsclrm_mod.F).

!       References:
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
    lstats_samp
 

implicit none

! Input Variables
character(len=*), intent(in) :: & 
  solve_type       ! Variable being solved; used for STATS.

logical, intent(in) :: & 
  l_first_clip_ts,    & ! First instance of clipping in a timestep.
  l_last_clip_ts     ! Last instance of clipping in a timestep.

real(kind=time_precision), intent(in) ::  & 
  dt     ! Model timestep; used here for STATS           [s]

real, dimension(gr%nnzp), intent(in) :: & 
  xp2,    & ! Variance of x, x'^2 (momentum levels)         [{x units}^2]
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
 

 
if ( lstats_samp ) then
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

 
if ( lstats_samp ) then
   if ( l_last_clip_ts ) then
      call stat_end_update( ixpyp_cl, real( xpyp / dt ), zm )
   else
      call stat_modify( ixpyp_cl, real( xpyp / dt ), zm )
   endif
endif
 


end subroutine covariance_clip

!===============================================================================
subroutine variance_clip( solve_type, dt, threshold, xp2 )

!       Description:
!       Clipping the value of variance x'^2 based on a minimum threshold
!       value.  The threshold value must be greater than or equal to 0.
!
!       The values of x'^2 are found on the momentum levels.
!
!       The following variances are found in the code:  
!
!       r_t'^2, th_l'^2, u'^2, v'^2, sclr'^2, (computed in diag_var.F);
!       w'^2 (computed in wp23.F).

!       References:
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
    lstats_samp
 

implicit none

! Input Variables
character(len=*), intent(in) :: & 
  solve_type  ! Variable being solved; used for STATS.

real(kind=time_precision), intent(in) :: & 
  dt          ! Model timestep; used here for STATS     [s]

real, dimension(gr%nnzp), intent(in) :: & 
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
 


 
if ( lstats_samp ) then
   call stat_begin_update( ixp2_cl, real( xp2 / dt ), zm )
endif
 

! Limit the value of x'^2 at threshold.
do k = 2, gr%nnzp, 1
   if ( xp2(k) < threshold(k) ) then
      xp2(k) = threshold(k)
   endif
enddo

 
if ( lstats_samp ) then
   call stat_end_update( ixp2_cl, real( xp2 / dt ), zm )
endif
 


end subroutine variance_clip

!===============================================================================
subroutine skewness_clip

!       Description:

!       References:
!-----------------------------------------------------------------------

implicit none


end subroutine skewness_clip

!===============================================================================

end module explicit_clip
