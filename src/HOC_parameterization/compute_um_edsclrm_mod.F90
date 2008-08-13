!------------------------------------------------------------------------
! $Id: compute_um_edsclrm_mod.F90,v 1.11 2008-08-13 16:52:10 griffinb Exp $
!------------------------------------------------------------------------
module compute_um_edsclrm_mod

implicit none

private ! Set Default Scope

public :: compute_um_edsclrm, compute_uv_tndcy

contains

!------------------------------------------------------------------------
subroutine compute_um_edsclrm( solve_type, xpwp_sfc, &
                               xm_tndcy,  & 
                               Kh_zm, dt,   &
                               xm, xpwp,  &
                               err_code )

! Description:
! Solves the horizontal wind or eddy-scalar time-tendency equation, and 
! diagnoses the turbulent flux.  A Crank-Nicholson time-stepping algorithm is
! used in solving the turbulent advection term and in diagnosing the turbulent
! flux.
!
! The rate of change of an eddy-scalar quantity, xm, is:
!
! d(xm)/dt = - w * d(xm)/dz - d(x'w')/dz + xm_forcings.
!
!
! The Turbulent Advection Term
! ----------------------------
!
! The above equation contains a turbulent advection term:
!
! - d(x'w')/dz;
!
! where the momentum flux, x'w', is closed using a down gradient approach:
!
! x'w' = - K_zm * d(xm)/dz.
!
! The turbulent advection term becomes:
!
! + d [ K_zm * d(xm)/dz ] / dz;
!
! which is the same as an eddy-diffusion term.  Thus, the turbulent advection
! term is treated and solved in the same way that an eddy-diffusion term would 
! be solved.  The term is discretized as follows:
!
! The values of xm are found on the thermodynamic levels, while the values of 
! K_zm are found on the momentum levels.  The derivatives (d/dz) of xm are taken
! over the intermediate momentum levels.  At the intermediate momentum levels, 
! d(xm)/dz is multiplied by K_zm.  Then, the derivative of the whole 
! mathematical expression is taken over the central thermodynamic level, which 
! yields the desired result.
!
! ---xm(kp1)----------------------------------------------- t(k+1)
!
! =============d(xm)/dz===K_zm(k)========================== m(k)
!
! ---xm(k)---------------------------d[K_zm*d(xm)/dz]/dz--- t(k)
!
! =============d(xm)/dz===K_zm(km1)======================== m(k-1)
!
! ---xm(km1)----------------------------------------------- t(k-1)
!
! The vertical indices t(k+1), m(k), t(k), m(k-1), and t(k-1) correspond with 
! altitudes zt(k+1), zm(k), zt(k), zm(k-1), and zt(k-1), respectively.  The 
! letter "t" is used for thermodynamic levels and the letter "m" is used for 
! momentum levels.
!
! dzt(k)   = 1 / ( zm(k) - zm(k-1) )
! dzm(k)   = 1 / ( zt(k+1) - zt(k) )
! dzm(k-1) = 1 / ( zt(k) - zt(k-1) )
!
! The vertically discretized form of the term is written out as:
!
! + dzt(k) * [   K_zm(k) * dzm(k) * ( xm(k+1) - xm(k) ) 
!              - K_zm(k-1) * dzm(k-1) * ( xm(k) - xm(k-1) ) ].
!
! For this equation, a Crank-Nicholson (semi-implicit) diffusion scheme is used
! to solve the d [ K_zm * d(xm)/dz ] / dz eddy-diffusion term.  The discretized 
! implicit form of the term is written out as:
!
! + (1/2)*dzt(k) * [   K_zm(k) * dzm(k) * ( xm(k+1,<t+1>) - xm(k,<t+1>) )
!                    - K_zm(k-1) * dzm(k-1) * ( xm(k,<t+1>) - xm(k-1,<t+1>) ) ].
!
! Note:  When the implicit term is brought over to the left-hand side, the 
!        sign is reversed and the leading "+" in front of the term is changed 
!        to a "-".
!
! The discretized explicit form of the term is written out as:
!
! + (1/2)*dzt(k) * [   K_zm(k) * dzm(k) * ( xm(k+1,<t>) - xm(k,<t>) )
!                    - K_zm(k-1) * dzm(k-1) * ( xm(k,<t>) - xm(k-1,<t>) ) ].
!
! Timestep index (t) stands for the index of the current timestep, while 
! timestep index (t+1) stands for the index of the next timestep, which is being
! advanced to in solving the d(xm)/dt equation.
!
!
! Boundary Conditions:
!
! An eddy-scalar quantity is not allowed to flux out the upper boundary.  Thus,
! a zero-flux boundary condition is used for the upper boundary in the eddy
! diffusion equation.
!
! The lower boundary condition is much more complicated.  It is neither a 
! zero-flux nor a fixed-point boundary condition.  Rather, it is a fixed-flux 
! boundary condition.  This term is a turbulent advection term, but with the
! eddy-scalars, the only value of x'w' relevant in solving the d(xm)/dt equation
! is the value of x'w' at the surface (written as x'w'|_sfc).
!
! Since x'w' = - K_zm * d(xm)/dz,
!
! x'w'|_sfc = - K_zm(1) * dzm(1) * ( xm(2) - xm(1) )
!
! The lower boundary condition, which in this case is applied to the d(xm)/dt 
! equation at level 2, is discretized as follows:
!
! ---xm(3)------------------------------------------------- t(3)
!
! =============d(xm)/dz===K_zm(2)========================== m(2)
!
! ---xm(2)---------------------------d[K_zm*d(xm)/dz]/dz--- t(2)
!
! =============[ x'w'|_sfc = - K_zm(1) * d(xm)/dz ]======== m(1) (surface)
!
! ---xm(1)------------------------------------------------- t(1)
! 
! The vertically discretized form of the term is written out as:
!
! + dzt(2) * [ K_zm(2) * dzm(2) * ( xm(3) - xm(2) ) + x'w'|_sfc ];
!
! which can be re-written as:
!
! + dzt(2) * [ K_zm(2) * dzm(2) * ( xm(3) - xm(2) ) ]  +  dzt(2) * x'w'|_sfc ].
!
! For this equation, a Crank-Nicholson (semi-implicit) diffusion scheme is used
! to solve the d [ K_zm * d(xm)/dz ] / dz eddy-diffusion term.  The discretized
! implicit form of the term is written out as:
!
! + (1/2)*dzt(2) * [ K_zm(2) * dzm(2) * ( xm(3,<t+1>) - xm(2,<t+1>) ) ].
!
! Note:  When the implicit term is brought over to the left-hand side, the
!        sign is reversed and the leading "+" in front of the term is changed
!        to a "-".
!
! The discretized explicit form of the term is written out as:
!
! + (1/2)*dzt(2) * [ K_zm(2) * dzm(2) * ( xm(3,<t>) - xm(2,<t>) ) ]
!    + dzt(2) * x'w'|_sfc.
!
! Note:  The x'w'|_sfc portion of the term written above has been pulled away 
!        from the rest of the explicit form written above because the (1/2) 
!        factor due to Crank-Nicholson time_stepping does not apply to it, as 
!        there isn't an implicit portion for x'w'|_sfc.
!
! Timestep index (t) stands for the index of the current timestep, while
! timestep index (t+1) stands for the index of the next timestep, which is being
! advanced to in solving the d(xm)/dt equation.
!
! The lower boundary condition for the implicit and explicit portions of the 
! term, without the x'w'|_sfc portion of the term, can easily be invoked by 
! using the zero-flux boundary conditions found in the generalized diffusion 
! function (function diffusion_zt_lhs) used for many other equations in this 
! model.  The dzt(2) * x'w'|_sfc term needs to be added onto the explicit term
! after this function has been called.  However, all other equations in this 
! model that use zero-flux diffusion have level 1 as the level to which the
! lower boundary condition is applied.  Thus, an adjuster will have to be used
! at level 2 to call diffusion_zt_lhs with level 1 as the input level (the last
! variable being passed in during the function call).  However, the other
! variables passed in (K_zm, gr%dzt, and gr%dzm variables) will have to be 
! passed in as solving for level 2.
!
! The value of xm(1) will just be set equal to itself and not effect the rest
! of the xm array during solving.  However, the value of xm(1) can be found
! after solving, since x'w'|_sfc = - K_zm(1) * dzm(1) * ( xm(2) - xm(1) ):
!
! xm(1)  =  xm(2)  +  x'w'|_sfc / ( K_zm(1) * dzm(1) ).

! References:
! Eqn. 8 & 9 on p. 3545 of 
! ``A PDF-Based Model for Boundary Layer Clouds. Part I:
!   Method and Model Description'' Golaz, et al. (2002)
! JAS, Vol. 59, pp. 3540--3551.
!-----------------------------------------------------------------------

use grid_class, only: & 
    gr ! Variable(s)

use lapack_wrap, only:  & 
    tridag_solve ! Procedure(s)
 
use stats_variables, only: & 
    ium_bt,  & ! Variable(s)
    ium_ta, & 
    ivm_bt, & 
    ivm_ta, & 
    zt, & 
    ztscr01, & 
    ztscr02, & 
    ztscr03, & 
    l_stats_samp

use stats_type, only: stat_update_var_pt, stat_modify_pt,  & 
                stat_begin_update, stat_end_update_pt

use constants, only:  & 
    fstderr ! Variable(s)
use stats_precision, only:  & 
    time_precision ! Variable(s)
use error_code, only:  & 
    lapack_error,  & ! Procedure(s)
    clubb_at_debug_level

implicit none

! Input Variables
character(len=*), intent(in) :: & 
  solve_type ! Desc. of what is being solved for

real(kind=time_precision), intent(in) ::  & 
  dt           ! Timestep                             [s]

real, intent(in) ::  & 
  xpwp_sfc     ! sfc flux                             [units vary]

real, dimension(gr%nnzp), intent(in) ::  & 
  xm_tndcy,  & ! x tendency                              [units vary]
  Kh_zm       ! Diffusion coefficient on momentum grid  [m^2/s]

! Input/Output Variables
real, dimension(gr%nnzp), intent(inout) ::  & 
  xm ! Prognostic array on the thermodynamic grid     [units vary]

! Output Variables
real, dimension(gr%nnzp), intent(inout) ::  & 
  xpwp   ! Momentum flux                              [units vary]

integer, intent(out) :: & 
  err_code ! clubb_singular_matrix when matrix is singular

! Local Variables
real,  dimension(gr%nnzp) :: a, b, c, rhs

real :: atmp, ctmp

integer :: k, kp1, km1 ! Indices

 
integer :: ixm_ta, ixm_bt

select case ( trim( solve_type ) )
case ( "um" )
  ixm_bt = ium_bt
  ixm_ta = ium_ta
case ( "vm" )
  ixm_bt = ivm_bt
  ixm_ta = ivm_ta
case default  ! Eddy scalars
  ixm_bt = 0
  ixm_ta = 0
end select

if ( l_stats_samp ) then
 
  ! xm total time tendency ( 1st calculation)
  call stat_begin_update( ixm_bt, real( xm / dt ), zt )
end if
 

! Prepare tridiagonal system

! zt(1) is below ground, we don't have to worry about it

a(1)   = 0.
b(1)   = real( 1./dt )
c(1)   = 0.
rhs(1) = real( 1./dt )

! zt(2) is the first active model layer. We need to impose the 
! surface momentum flux xpwp_sfc.

atmp = 0.
ctmp = -0.5 * Kh_zm(2) * gr%dzm(2) * gr%dzt(2)

a(2)   = atmp
c(2)   = ctmp
b(2)   = real( - ctmp + 1./dt )
rhs(2) = real( ( ctmp + 1./dt ) * xm(2) & 
         - ctmp * xm(3) & 
         + xm_tndcy(2) & 
         + xpwp_sfc * gr%dzt(2) )

   if ( l_stats_samp .and. ixm_ta > 0 ) then
 
     ztscr01(1) = 0.0
     ztscr02(1) = 0.0
     ztscr03(1) = 0.0

     ztscr01(2) = -atmp
     ztscr02(2) = ctmp
     ztscr03(2) = -ctmp
   end if
! Loop from level 3 to gr%nnzp

do k=3, gr%nnzp-1, 1

   atmp = -0.5 * Kh_zm(k-1) * gr%dzt(k) * gr%dzm(k-1)
   ctmp = -0.5 * Kh_zm(k) * gr%dzt(k) * gr%dzm(k)

   a(k)   = atmp
   c(k)   = ctmp
   b(k)   = real( - atmp - ctmp + 1./dt )
   rhs(k) = real( - atmp * xm(k-1) & 
            + ( atmp + ctmp + 1./dt ) * xm(k) & 
            - ctmp * xm(k+1) & 
            + xm_tndcy(k) )
   if ( l_stats_samp .and. ixm_ta > 0 ) then
 
     ztscr01(k) = -atmp
     ztscr02(k) =  atmp + ctmp
     ztscr03(k) = -ctmp
   end if
end do

! Level gr%nnzp. We impose zero flux from model top

! Vince Larson simplified the upper BC.  It still imposes
!    zero flux. 7 Jul 2007 
!        atmp = -0.5 * Kh_zm(gr%nnzp-1) * gr%dzm(gr%nnzp) 
!     .              * gr%dzt(gr%nnzp-1)
!        ctmp = 0.

!        a(gr%nnzp)   = atmp
!        c(gr%nnzp)   = ctmp
!        b(gr%nnzp)   = - atmp + 1./dt
!        rhs(gr%nnzp) = - atmp * xm(gr%nnzp-1)
!     .               + ( atmp + 1./dt ) * xm(gr%nnzp)
atmp = -1.
ctmp = 0.

a(gr%nnzp)   =  atmp
c(gr%nnzp)   =  ctmp
b(gr%nnzp)   =  -atmp
rhs(gr%nnzp) =  0.
! End of Vince Larson's change


 
   ! Zero flux
   ! This new code should make the budget balance at nnzp
   ! -dschanen 18 Jul 2008
   if ( l_stats_samp .and. ixm_ta > 0 ) then
     ztscr01(gr%nnzp) = -atmp
     ztscr02(gr%nnzp) =  ctmp
     ztscr03(gr%nnzp) =  atmp
   end if

   
!    Caused problems with DYCOMS II RF02
!     .               + xm_tndcy(gr%nnzp)

!    Attempted to compensate for the DYCOMS problem using the code below.
!    Doesn't actually seem to make a difference
!       atmp = -0.5 * Kh_zm(gr%nnzp-1) * gr%dzt(gr%nnzp) * gr%dzm(gr%nnzp-1)

!       a(gr%nnzp)   = atmp
!       b(gr%nnzp)   = - c(gr%nnzp-1) + 1./dt
!       c(gr%nnzp)   = UNDEFINED
!       rhs(gr%nnzp) = - atmp * x(gr%nnzp-1)
!    .                 + ( atmp + 1./dt ) * x(gr%nnzp)
!    .                 + xm_tndcy(gr%nnzp)

! Store momentum flux (explicit component)
xpwp(1) = xpwp_sfc
do k=2,gr%nnzp-1
  xpwp(k) = -0.5 * Kh_zm(k) * gr%dzm(k) * ( xm(k+1) - xm(k) ) 
end do
xpwp(gr%nnzp) = 0.

 
   ! Turbulent transport (explicit component)
   if ( l_stats_samp .and. ixm_ta > 0 ) then
     do k=1,gr%nnzp,1
       km1 = max( k-1, 1 )
       kp1 = min( k+1, gr%nnzp )

       call stat_modify_pt( ixm_ta, k, & 
         ztscr01(k) * xm(km1) & 
       + ztscr02(k) * xm(k) & 
       + ztscr03(k) * xm(kp1), zt )

     end do
   end if

! Solve tridiagonal system
call tridag_solve( solve_type, gr%nnzp, 1, c, b, a,  & 
                   rhs, xm, err_code )

! Error report
! Joshua Fasching February 2008
if ( lapack_error( err_code ) .and.  & 
     clubb_at_debug_level( 1 ) ) then
               
    write(fstderr,*) "Error in compute_um_edsclrm"
   
    write(fstderr,*) "Intent(in)"
   
    write(fstderr,*) "dt = ", dt
    write(fstderr,*) "xpwp_sfc = ", xpwp_sfc
    write(fstderr,*) "xm_tndcy = ", xm_tndcy
    write(fstderr,*) "Kh_zm = ", Kh_zm
   
    write(fstderr,*) "Intent(inout)"
   
    write(fstderr,*) "xm = ", xm
   
    write(fstderr,*) "Intent(out)"
   
    write(fstderr,*) "xpwp = ", xpwp
             
  return
end if

! Second part of momentum (implicit component)
do k=2, gr%nnzp-1, 1
  xpwp(k) = xpwp(k)  & 
             - 0.5 * Kh_zm(k) * gr%dzm(k) * ( xm(k+1) - xm(k) )
end do

if ( l_stats_samp ) then
 
  do k = 1, gr%nnzp, 1
    km1 = max( k-1, 1 )
    kp1 = min( k+1, gr%nnzp )

    ! x time tendency (2nd calculation)
    call stat_end_update_pt( ixm_bt, k, real( xm(k) / dt), zt )

    ! x turbulent transport (implicit component)
    call stat_update_var_pt( ixm_ta, k, & 
        ztscr01(k) * xm(km1) & 
      + ztscr02(k) * xm(k) & 
      + ztscr03(k) * xm(kp1), zt )
    
   end do

 end if

return

end subroutine compute_um_edsclrm

!-----------------------------------------------------------------------
subroutine compute_uv_tndcy( solve_type, xm, wm_zt, fcor, perp_wind_m, perp_wind_g, implemented, & 
                             xmt )
!
!       Description: Computes the tendency for the u/v wind components.
!
!
!-----------------------------------------------------------------------
  use grid_class, only: & 
    gr

  use grid_class, only: & 
    zt2zm, & 
    ddzm         

 
  use stats_type, only: & 
    stat_update_var

  use stats_variables, only:      & 
    ium_ma,  & ! Variable(s)
    ium_gf, & 
    ium_cf, & 
    ivm_ma, & 
    ivm_gf, & 
    ivm_cf, & 
    zt, & 
    l_stats_samp

  implicit none

  ! Input Variables
  character(len=*), intent(in) :: solve_type

  real, dimension(gr%nnzp), intent(in) ::  & 
    xm,  & ! u/v wind                                       [m/s]   
    wm_zt ! wm on thermodynaming grid                     [m/s]

  real, intent(in) ::  & 
    fcor ! Coriolis forcing                             [s^-1]

  real, dimension(gr%nnzp), intent(in) :: & 
    perp_wind_m,  & ! Perpendicular component of the mean wind (e.g. v, for the u-eqn) [m/s]
    perp_wind_g  ! Perpendicular component of the geostropic wind (e.g. vg)         [m/s]

  logical, intent(in) :: & 
    implemented

  ! Output Variables
  real, dimension(gr%nnzp), intent(out) :: xmt ! xm tendency  [m/s^2]

  ! Local Variables
  integer :: & 
    ixm_ma, & 
    ixm_gf, & 
    ixm_cf

  real, dimension(gr%nnzp) :: & 
    xm_ma, & 
    xm_gf, & 
    xm_cf


if (.not. implemented) then
  select case (trim(solve_type))
  case("um")
 
    ixm_ma = ium_ma
    ixm_gf = ium_gf
    ixm_cf = ium_cf

    xm_ma = - wm_zt * ddzm( zt2zm( xm ) )

    xm_gf = - fcor * perp_wind_g

    xm_cf = fcor * perp_wind_m

  case("vm") 
 
    ixm_ma = ivm_ma
    ixm_gf = ivm_gf
    ixm_cf = ivm_cf

    xm_ma = - wm_zt * ddzm( zt2zm( xm ) )

    xm_gf = fcor * perp_wind_g

    xm_cf = -fcor * perp_wind_m
  case default
    ixm_ma = 0
    ixm_gf = 0
    ixm_cf = 0
    xm_ma = 0.
    xm_gf = 0.
    xm_cf = 0.
  end select

  xmt = xm_ma + xm_gf + xm_cf 

  if ( l_stats_samp ) then
 
    call stat_update_var( ixm_ma, xm_ma, zt )

    call stat_update_var( ixm_gf, xm_gf, zt )
  
    call stat_update_var( ixm_cf, xm_cf, zt )
  endif
else
  xmt = 0.0
endif

end subroutine compute_uv_tndcy

!===============================================================================
subroutine compute_um_edsclrm_lhs( solve_type, dt, wm_zt, Kh_zm, lhs )

! Description:
! Calculate the implicit portion of the horizontal wind or eddy-scalar 
! time-tendency equation.  See the description in subroutine compute_um_edsclrm
! for more details.

! References:
!-----------------------------------------------------------------------

use grid_class, only:  & 
    gr  ! Variable(s)

use stats_precision, only:  & 
    time_precision ! Variable(s)

use diffusion, only:  & 
    diffusion_zt_lhs ! Procedure(s)

use mean_adv, only: & 
    term_ma_zt_lhs  ! Procedures

use stats_variables, only: &
    ium_ma,  & ! Variable(s)
    ium_ta,  &
    ivm_ma,  &
    ivm_ta,  &
    ztscr01, &
    ztscr02, &
    ztscr03, &
    ztscr04, &
    ztscr05, &
    ztscr06, &
    l_stats_samp

implicit none

! Constant parameters
integer, parameter :: & 
  kp1_tdiag = 1,    & ! Thermodynamic superdiagonal index.
  k_tdiag   = 2,    & ! Thermodynamic main diagonal index.
  km1_tdiag = 3       ! Thermodynamic subdiagonal index.

! Input Variables
character(len=*), intent(in) :: &
  solve_type ! Desc. of what is being solved for

real(kind=time_precision), intent(in) :: & 
  dt         ! Model timestep                           [s]

real, dimension(gr%nnzp), intent(in) :: &
  wm_zt,   & ! w wind component on thermodynamic levels [m/s]
  Kh_zm      ! Eddy diffusivity on momentum levels      [m^2/s]

! Output Variable
real, dimension(3,gr%nnzp), intent(out) :: &
  lhs        ! Implicit contributions to the term

! Local Variables
integer :: k, km1  ! Array indices
integer :: diff_k_in

real, dimension(3) :: tmp
integer :: ixm_ma, ixm_ta

select case ( trim( solve_type ) )
   case ( "um" )
      ixm_ma = ium_ma
      ixm_ta = ium_ta
   case ( "vm" )
      ixm_ma = ivm_ma
      ixm_ta = ivm_ta
   case default  ! Eddy scalars
      ixm_ma = 0
      ixm_ta = 0
end select

! Initialize the LHS array.
lhs = 0.0

do k = 2, gr%nnzp-1, 1

   ! Define index
   km1 = max( k-1, 1 )

   ! LHS mean advection term.
   lhs(kp1_tdiag:km1_tdiag,k)  &
   = lhs(kp1_tdiag:km1_tdiag,k)  &
   + term_ma_zt_lhs( wm_zt(k), gr%dzt(k), k )

   ! LHS turbulent advection term (solved as an eddy-diffusion term).
   if ( k == 2 ) then
      ! The lower boundary condition needs to be applied here at level 2.
      diff_k_in = 1
   else
      diff_k_in = k
   endif
   lhs(kp1_tdiag:km1_tdiag,k)  &
   = lhs(kp1_tdiag:km1_tdiag,k)  &
   + (1.0/2.0)  &
   * diffusion_zt_lhs( Kh_zm(k), Kh_zm(km1), 0.0,  & 
                       gr%dzm(km1), gr%dzm(k), gr%dzt(k), diff_k_in )

   ! LHS time tendency.
   lhs(k_tdiag,k)  &
   = real( lhs(k_tdiag,k) + ( 1.0 / dt ) )

   if ( l_stats_samp ) then

      ! Statistics:  implicit contributions for um or vm.

      if ( ixm_ma > 0 ) then
         tmp(1:3) &
         = term_ma_zt_lhs( wm_zt(k), gr%dzt(k), k )
         ztscr01(k) = -tmp(3)
         ztscr02(k) = -tmp(2)
         ztscr03(k) = -tmp(1)
      endif

      if ( ixm_ta > 0 ) then
         tmp(1:3)  &
         = (1.0/2.0)  &
         * diffusion_zt_lhs( Kh_zm(k), Kh_zm(km1), 0.0,  &
                             gr%dzm(km1), gr%dzm(k), gr%dzt(k), diff_k_in )
         ztscr04(k) = -tmp(3)
         ztscr05(k) = -tmp(2)
         ztscr06(k) = -tmp(1)
      endif

   endif  ! lstats_samp

enddo


! Boundary Conditions
! Lower Boundary
lhs(:,1) = 0.0
lhs(k_tdiag,1) = real( 1.0 / dt )

! Upper Boundary
lhs(:,gr%nnzp) = 0.0
lhs(k_tdiag,gr%nnzp) = real( 1.0 / dt )


return
end subroutine compute_um_edsclrm_lhs

!===============================================================================
subroutine compute_um_edsclrm_rhs( solve_type, dt, Kh_zm, xm,  &
                                   xm_forcing, xpwp_sfc, rhs )

! Description:
! Calculate the explicit portion of the horizontal wind or eddy-scalar 
! time-tendency equation.  See the description in subroutine compute_um_edsclrm
! for more details.

! References:
!-----------------------------------------------------------------------

use grid_class, only:  & 
    gr  ! Variable(s)

use stats_precision, only:  & 
    time_precision ! Variable(s)

use diffusion, only:  & 
    diffusion_zt_lhs ! Procedure(s)

use stats_variables, only: &
    ium_ta,  & ! Variable(s)
    ivm_ta,  &
    zt,      &
    l_stats_samp

use stats_type, only: &
    stat_begin_update_pt,  & ! Procedure(s)
    stat_modify_pt

implicit none

! Input Variables
character(len=*), intent(in) :: &
  solve_type ! Desc. of what is being solved for

real(kind=time_precision), intent(in) :: & 
  dt           ! Model timestep                                   [s]

real, dimension(gr%nnzp), intent(in) :: &
  Kh_zm,     & ! Eddy diffusivity on momentum levels              [m^2/s]
  xm,        & ! Eddy-scalar variable, xm (thermodynamic levels)  [units vary]
  xm_forcing   ! The explicit time-tendency acting on xm          [units vary]

real, intent(in) :: &
  xpwp_sfc     ! x'w' at the surface                              [units vary]

! Output Variable
real, dimension(gr%nnzp), intent(out) :: &
  rhs

! Local Variables
integer :: k, kp1, km1  ! Array indices
integer :: diff_k_in

! For use in Crank-Nicholson eddy diffusion.
real, dimension(3) :: rhs_diff

integer :: ixm_ta

select case ( trim( solve_type ) )
   case ( "um" )
      ixm_ta = ium_ta
   case ( "vm" )
      ixm_ta = ivm_ta
   case default  ! Eddy scalars
      ixm_ta = 0
end select


! Initialize the RHS vector.
rhs = 0.0

do k = 2, gr%nnzp-1, 1

   ! Define indices
   km1 = max( k-1, 1 )
   kp1 = min( k+1, gr%nnzp )

   ! RHS turbulent advection term (solved as an eddy-diffusion term).
   if ( k == 2 ) then
      ! The lower boundary condition needs to be applied here at level 2.
      diff_k_in = 1
   else
      diff_k_in = k
   endif
   rhs_diff(1:3)  & 
   = (1.0/2.0)  &
   * diffusion_zt_lhs( Kh_zm(k), Kh_zm(km1), 0.0,  &
                       gr%dzm(km1), gr%dzm(k), gr%dzt(k), diff_k_in )
   rhs(k)   =   rhs(k) & 
              - rhs_diff(3) * xm(km1) &
              - rhs_diff(2) * xm(k)   &
              - rhs_diff(1) * xm(kp1)

   ! RHS forcings.
   rhs(k) = rhs(k) + xm_forcing(k)

   ! RHS time tendency
   rhs(k) = real( rhs(k) + ( 1.0 / dt ) * xm(k) )

   if ( l_stats_samp ) then

      ! Statistics:  explicit contributions for um or vm.

      if ( ixm_ta > 0 ) then
         call stat_begin_update_pt( ixm_ta, k, & 
                rhs_diff(3) * xm(km1) &
              + rhs_diff(2) * xm(k)   &
              + rhs_diff(1) * xm(kp1), zt )
      endif

   endif  ! lstats_samp

enddo


! Boundary Condition
! Lower Boundary
rhs(2) = rhs(2) + gr%dzt(2) * xpwp_sfc
rhs(1) = real( ( 1.0 / dt ) * xm(1) )

if ( l_stats_samp ) then

   ! Statistics:  explicit contributions for um or vm.

   if ( ixm_ta > 0 ) then
      call stat_modify_pt( ixm_ta, 2, &
           + gr%dzt(2) * xpwp_sfc, zt )
   endif

endif  ! lstats_samp


! Upper Boundary
rhs(gr%nnzp) = real( ( 1.0 / dt ) * xm(gr%nnzp) )


return
end subroutine compute_um_edsclrm_rhs

!===============================================================================

end module compute_um_edsclrm_mod
