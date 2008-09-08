!------------------------------------------------------------------------
! $Id$
!===============================================================================
module advance_windm_edsclrm_module

implicit none

private ! Set Default Scope

public :: advance_windm_edsclrm

private :: compute_um_edsclrm_solve, &
           compute_uv_tndcy,         &
           compute_um_edsclrm_lhs,   &
           compute_um_edsclrm_rhs

contains

!===============================================================================
subroutine advance_windm_edsclrm( dt, wm_zt, Kh_zm, ug, vg, um_ref, vm_ref,  &
                                  wp2, up2, vp2, upwp_sfc, vpwp_sfc, fcor,  &
                                  l_implemented, um, vm, edsclrm,  &
                                  upwp, vpwp, wpedsclrp, err_code )

! Description:
! Solves for both mean horizontal wind components, um and vm, and for mean 
! eddy-scalars.

! References:
! Eqn. 8 & 9 on p. 3545 of 
! ``A PDF-Based Model for Boundary Layer Clouds. Part I:
!   Method and Model Description'' Golaz, et al. (2002)
! JAS, Vol. 59, pp. 3540--3551.
!-----------------------------------------------------------------------

use grid_class, only:  &
    gr  ! Variables(s)

use parameters_tunable, only:  &
    ts_nudge,  & ! Variable(s)
    sclr_dim

use model_flags, only:  &
    l_uv_nudge  ! Variable(s)

use stats_precision, only:  &
    time_precision  ! Variable(s)

use clip_explicit, only:  &
    clip_covariance  ! Procedure(s)

use error_code, only:  & 
    lapack_error,  & ! Procedure(s)
    clubb_at_least_debug_level

use constants, only:  & 
    fstderr  ! Constant

implicit none

! Input Variables
real(kind=time_precision), intent(in) ::  &
  dt             ! Model timestep                                [s]

real, dimension(gr%nnzp), intent(in) ::  &
  wm_zt,       & ! w wind component on thermodynamic levels      [m/s]
  Kh_zm,       & ! Eddy diffusivity on momentum levels           [m^2/s]
  ug,          & ! u (west-to-east) geostrophic wind component   [m/s]
  vg,          & ! v (south-to-north) geostrophic wind component [m/s]
  um_ref,      & ! Reference u wind component for nudging        [m/s]
  vm_ref,      & ! Reference v wind component for nudging        [m/s]
  wp2,         & ! w'^2 (momentum levels)                        [m^2/s^2]
  up2,         & ! u'^2 (momentum levels)                        [m^2/s^2]
  vp2            ! v'^2 (momentum levels)                        [m^2/s^2]

real, intent(in) ::  &
  upwp_sfc,    & ! u'w' at the surface (momentum level 1)        [m^2/s^2]
  vpwp_sfc,    & ! v'w' at the surface (momentum level 1)        [m^2/s^2]
  fcor           ! Coriolis parameter                            [s^-1]

logical, intent(in) ::  &
  l_implemented  ! Flag for CLUBB being implemented in a larger model.

! Input/Output Variables
real, dimension(gr%nnzp), intent(inout) ::  &
  um,          & ! mean u (west-to-east) wind component          [m/s]
  vm             ! mean v (south-to-north) wind component        [m/s]

! Input/Output Variable for eddy-scalars
real, dimension(gr%nnzp,sclr_dim), intent(inout) ::  &
  edsclrm        ! mean eddy scalar quantity                     [units vary]

! Output Variables
real, dimension(gr%nnzp), intent(out) ::  &
  upwp,        & ! u'w' (momentum levels)                        [m^2/s^2]
  vpwp           ! v'w' (momentum levels)                        [m^2/s^2]

! Output Variable for eddy-scalars
real, dimension(gr%nnzp,sclr_dim), intent(out) ::  &
  wpedsclrp      ! w'edsclr' (momentum levels)                   [units vary]

integer, intent(out) :: &
  err_code       ! clubb_singular_matrix when matrix is singular

! Local Variables
real, dimension(gr%nnzp) ::  &
  um_tndcy,    & ! u wind component tendency                     [m/s^2]
  vm_tndcy       ! v wind component tendency                     [m/s^2]

real, dimension(gr%nnzp,sclr_dim) ::  &
  edsclrm_tndcy  ! eddy-scalar tendency                          [units vary]

integer :: i     ! Array index


!----------------------------------------------------------------
! Update zonal (west-to-east) component of mean wind, um
!----------------------------------------------------------------

call compute_uv_tndcy( "um", fcor, vm, vg, &
                       l_implemented, um_tndcy )

call compute_um_edsclrm_solve( "um", dt, upwp_sfc, um_tndcy,  &
                               wm_zt, Kh_zm, l_implemented,  &
                               um, upwp, err_code )


!----------------------------------------------------------------
! Update meridional (south-to-north) component of mean wind, vm
!----------------------------------------------------------------

call compute_uv_tndcy( "vm", fcor, um, ug, &
                       l_implemented, vm_tndcy )

call compute_um_edsclrm_solve( "vm", dt, vpwp_sfc, vm_tndcy,  &
                               wm_zt, Kh_zm, l_implemented,  &
                               vm, vpwp, err_code )
       

! Adjust um and vm if nudging is turned on.
if ( l_uv_nudge ) then
  um(1:gr%nnzp) = real( um(1:gr%nnzp)  & 
        - ((um(1:gr%nnzp) - um_ref(1:gr%nnzp)) * (dt/ts_nudge)) )
  vm(1:gr%nnzp) = real( vm(1:gr%nnzp)  & 
        - ((vm(1:gr%nnzp) - vm_ref(1:gr%nnzp)) * (dt/ts_nudge)) )
endif


! Clipping for u'w'
!
! Clipping u'w' at each vertical level, based on the
! correlation of u and w at each vertical level, such that:
! corr_(u,w) = u'w' / [ sqrt(u'^2) * sqrt(w'^2) ];
! -1 <= corr_(u,w) <= 1.
! Since u'^2, w'^2, and u'w' are updated in different
! places from each other, clipping for u'w' has to be done
! three times.  This is the third instance of u'w' clipping.
call clip_covariance( "upwp", .false.,      & ! intent(in)
                      .true., dt, wp2, up2, & ! intent(in)
                      upwp )                  ! intent(inout)

! Clipping for v'w'
!
! Clipping v'w' at each vertical level, based on the
! correlation of v and w at each vertical level, such that:
! corr_(v,w) = v'w' / [ sqrt(v'^2) * sqrt(w'^2) ];
! -1 <= corr_(v,w) <= 1.
! Since v'^2, w'^2, and v'w' are updated in different
! places from each other, clipping for v'w' has to be done
! three times.  This is the third instance of v'w' clipping.
call clip_covariance( "vpwp", .false.,      & ! intent(in)
                      .true., dt, wp2, vp2, & ! intent(in)
                      vpwp )                  ! intent(inout)


!----------------------------------------------------------------
! Compute Eddy-diff. Passive Scalars
!----------------------------------------------------------------

if ( sclr_dim > 0 ) then

   do i = 1, sclr_dim

      edsclrm_tndcy(1:gr%nnzp,i) = 0.0

      call compute_um_edsclrm_solve( "edsclr", dt, wpedsclrp(1,i),  &
                                     edsclrm_tndcy(:,i), wm_zt, Kh_zm,  &
                                     l_implemented, edsclrm(:,i),  &
                                     wpedsclrp(:,i), err_code )

   enddo

endif ! sclr_dim > 0


! Error report
! Joshua Fasching February 2008
if ( lapack_error( err_code ) .and.  &
     clubb_at_least_debug_level( 1 ) ) then

    write(fstderr,*) "Error in compute_um_edsclrm"

    write(fstderr,*) "Intent(in)"

    write(fstderr,*) "dt = ", dt
    write(fstderr,*) "wm_zt = ", wm_zt
    write(fstderr,*) "Kh_zm = ", Kh_zm
    write(fstderr,*) "ug = ", ug
    write(fstderr,*) "vg = ", vg
    write(fstderr,*) "um_ref = ", um_ref
    write(fstderr,*) "vm_ref = ", vm_ref
    write(fstderr,*) "wp2 = ", wp2
    write(fstderr,*) "up2 = ", up2
    write(fstderr,*) "vp2 = ", vp2
    write(fstderr,*) "upwp_sfc = ", upwp_sfc
    write(fstderr,*) "vpwp_sfc = ", vpwp_sfc
    write(fstderr,*) "fcor = ", fcor
    write(fstderr,*) "l_implemented = ", l_implemented

    write(fstderr,*) "Intent(inout)"

    write(fstderr,*) "um = ", um
    write(fstderr,*) "vm = ", vm
    write(fstderr,*) "edsclrm = ", edsclrm

    write(fstderr,*) "Intent(out)"

    write(fstderr,*) "upwp = ", upwp
    write(fstderr,*) "vpwp = ", vpwp
    write(fstderr,*) "wpedsclrp = ", wpedsclrp

    return

endif


return
end subroutine advance_windm_edsclrm

!===============================================================================
subroutine compute_um_edsclrm_solve( solve_type, dt, xpwp_sfc, xm_tndcy,  &
                                     wm_zt, Kh_zm, l_implemented,  &
                                     xm, xpwp, err_code )

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
! The vertically discretized form of the turbulent advection term (treated as an
! eddy diffusion term) is written out as:
!
! + dzt(k) * [   K_zm(k) * dzm(k) * ( xm(k+1) - xm(k) ) 
!              - K_zm(k-1) * dzm(k-1) * ( xm(k) - xm(k-1) ) ].
!
! For this equation, a Crank-Nicholson (semi-implicit) diffusion scheme is used
! to solve the d [ K_zm * d(xm)/dz ] / dz eddy-diffusion term.  The discretized 
! implicit portion of the term is written out as:
!
! + (1/2)*dzt(k) * [   K_zm(k) * dzm(k) * ( xm(k+1,<t+1>) - xm(k,<t+1>) )
!                    - K_zm(k-1) * dzm(k-1) * ( xm(k,<t+1>) - xm(k-1,<t+1>) ) ].
!
! Note:  When the implicit term is brought over to the left-hand side, the 
!        sign is reversed and the leading "+" in front of the term is changed 
!        to a "-".
!
! The discretized explicit portion of the term is written out as:
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
! is the value of x'w' at the surface (the first momentum level), which is 
! written as x'w'|_sfc.
!
! Since x'w' = - K_zm * d(xm)/dz, a substitution can be made in the xm eddy 
! diffusion equation for - K_zm * d(xm/dz) at the surface, such that:
!
! - K_zm(1) * dzm(1) * ( xm(2) - xm(1) ) = x'w'|_sfc;
!
! where x'w'|_sfc is the surface momentum flux that is computed elsewhere in the
! model code.
!
! The lower boundary condition, which in this case needs to be applied to the 
! d(xm)/dt equation at level 2, is discretized as follows:
!
! ---xm(3)------------------------------------------------- t(3)
!
! =============[ + K_zm(2) * d(xm)/dz ]==================== m(2)
!
! ---xm(2)---------------------------d[K_zm*d(xm)/dz]/dz--- t(2)
!
! =============[ - K_zm(1) * d(xm)/dz = x'w'|_sfc ]======== m(1) (surface)
!
! ---xm(1)------------------------------------------------- t(1)
! 
! The vertically discretized form of the turbulent advection term (treated as an
! eddy diffusion term), with the x'w'|_sfc substitution in place, is written 
! out as:
!
! + dzt(2) * [ K_zm(2) * dzm(2) * ( xm(3) - xm(2) ) + x'w'|_sfc ];
!
! which can be re-written as:
!
! + dzt(2) * [ K_zm(2) * dzm(2) * ( xm(3) - xm(2) ) ]  +  dzt(2) * x'w'|_sfc.
!
! For this equation, a Crank-Nicholson (semi-implicit) diffusion scheme is used
! to solve the d [ K_zm * d(xm)/dz ] / dz eddy-diffusion term.  The discretized
! implicit portion of the term is written out as:
!
! + (1/2)*dzt(2) * [ K_zm(2) * dzm(2) * ( xm(3,<t+1>) - xm(2,<t+1>) ) ].
!
! Note:  When the implicit term is brought over to the left-hand side, the
!        sign is reversed and the leading "+" in front of the term is changed
!        to a "-".
!
! The discretized explicit portion of the term is written out as:
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
! The value of xm(1) is located below the model surface and is not treated as a
! time-tendency equation, but rather is solved for based on the result of xm(2) 
! after xm(2) has been advanced to the next timestep.  Since, 
! x'w'|_sfc = - K_zm(1) * dzm(1) * ( xm(2) - xm(1) ):
!
! xm(1)  =  xm(2)  +  x'w'|_sfc / ( K_zm(1) * dzm(1) ).
!
! This equation can be broken down into implicit and explicit components.  The
! implicit portion of this equation is:
!
! [+1] * xm(1,<t+1>) + [-1] * xm(2,<t+1>);
!
! while the explicit portion of this equation is:
!
! + x'w'|_sfc / ( K_zm(1) * dzm(1) ).
!
!
! Conservation Properties:
!
! When a fixed-flux lower boundary condition is used (combined with a zero-flux 
! upper boundary condition), this technique of discretizing the turbulent
! advection term (treated as an eddy-diffusion term) leads to conservative 
! differencing.  The column totals for each column in the left-hand side matrix 
! (for the turbulent advection term) should be equal to 0, while the column 
! total for the right-hand side vector (for the turbulent advection term) should
! be equal to the surface flux.  This ensures that the total amount of quantity 
! xm over the entire vertical domain is only changed by the surface flux 
! (neglecting any forcing terms).  The total amount of change is equal to the 
! surface flux.
!
! To see that this conservation law is satisfied by the left-hand side matrix,
! compute the turbulent advection (treated as eddy diffusion) of xm and 
! integrate vertically.  In discretized matrix notation (where "i" stands for 
! the matrix column and "j" stands for the matrix row):
!
!  0 = Sum_j Sum_i ( 1/dzt )_i ( 0.5 * dzt * (K_zm*dzm) )_ij (xm<t+1>)_j.
!
! The left-hand side matrix, ( 0.5 * dzt * (K_zm*dzm) )_ij, is partially written
! below.  The sum over i in the above equation removes dzt everywhere from the 
! matrix below.  The sum over j leaves the column totals that are desired, which
! are 0.
!
! Left-hand side matrix contributions from the turbulent advection term
! (treated as an eddy-diffusion term using a Crank-Nicholson timestep); 
! first five vertical levels:
!
!     ------------------------------------------------------------------------------->
!k=1 |  0             0                        0                          0
!    |
!k=2 |  0   +0.5*dzt(k)*           -0.5*dzt(k)*                           0
!    |        K_zm(k)*dzm(k)         K_zm(k)*dzm(k)
!    | 
!k=3 |  0   -0.5*dzt(k)*           +0.5*dzt(k)*               -0.5*dzt(k)*
!    |        K_zm(k-1)*dzm(k-1)     [ K_zm(k)*dzm(k)           K_zm(k)*dzm(k)
!    |                                +K_zm(k-1)*dzm(k-1) ]
!    |
!k=4 |  0             0            -0.5*dzt(k)*               +0.5*dzt(k)*
!    |                               K_zm(k-1)*dzm(k-1)         [ K_zm(k)*dzm(k)
!    |                                                           +K_zm(k-1)*dzm(k-1) ]
!    |
!k=5 |  0             0                        0              -0.5*dzt(k)*
!    |                                                          K_zm(k-1)*dzm(k-1)
!   \ /
!
! Note:  The superdiagonal term from level 4 and both the main diagonal and 
!        superdiagonal terms from level 5 are not shown on this diagram.
!
! To see that the above conservation law is satisfied by the right-hand side 
! vector, compute the turbulent advection (treated as eddy diffusion) of xm and
! integrate vertically.  In discretized matrix notation (where "i" stands for 
! the matrix column and "j" stands for the matrix row):
!
!  x'w'|_sfc = Sum_j Sum_i ( 1/dzt )_i ( rhs_vector )_j.
!
! The right-hand side vector, ( rhs_vector )_j, is partially written below.
! The sum over i in the above equation removes dzt everywhere from the 
! vector below.  The sum over j leaves the column total that is desired, 
! which is x'w'|_sfc.
!
! Right-hand side vector contributions from the turbulent advection term
! (treated as an eddy-diffusion term using a Crank-Nicholson timestep);
! first five vertical levels:
!
!     --------------------------------------------
!k=1 |                      0                     |
!    |                                            |
!    |                                            |
!k=2 | +0.5*dzt(k)*                               |
!    |       [ K_zm(k)*dzm(k)*                    |
!    |                  (xm(k+1,<t>)-xm(k,<t>)) ] |
!    | +dzt(k) * x'w'|_sfc                        |
!    |                                            |
!k=3 | +0.5*dzt(k)*                               |
!    |       [ K_zm(k)*dzm(k)*                    |
!    |                  (xm(k+1,<t>)-xm(k,<t>))   |
!    |        -K_zm(k-1)*dzm(k-1)*                |
!    |                  (xm(k,<t>)-xm(k-1,<t>)) ] |
!    |                                            |
!k=4 | +0.5*dzt(k)*                               |
!    |       [ K_zm(k)*dzm(k)*                    |
!    |                  (xm(k+1,<t>)-xm(k,<t>))   |
!    |        -K_zm(k-1)*dzm(k-1)*                |
!    |                  (xm(k,<t>)-xm(k-1,<t>)) ] |
!    |                                            |
!k=5 | +0.5*dzt(k)*                               |
!    |       [ K_zm(k)*dzm(k)*                    |
!    |                  (xm(k+1,<t>)-xm(k,<t>))   |
!    |        -K_zm(k-1)*dzm(k-1)*                |
!    |                  (xm(k,<t>)-xm(k-1,<t>)) ] |
!   \ /                                          \ /
!
! Note:  Only the contributions by the turbulent advection term are shown for
!        both the left-hand side matrix and the right-hand side vector.  There
!        are more terms in the equation, and thus more factors to be added to
!        both the left-hand side matrix (such as time tendency and mean 
!        advection) and the right-hand side vector (such as xm forcings).  The
!        left-hand side matrix is set-up so that a singular matrix is not 
!        encountered.

! References:
! Eqn. 8 & 9 on p. 3545 of 
! ``A PDF-Based Model for Boundary Layer Clouds. Part I:
!   Method and Model Description'' Golaz, et al. (2002)
! JAS, Vol. 59, pp. 3540--3551.
!-----------------------------------------------------------------------

use grid_class, only: & 
    gr ! Variable(s)

use lapack_wrap, only:  & 
  tridag_solve, & ! Procedure(s)
  tridag_solvex 

use stats_variables, only: & 
    ium_bt,  & ! Variable(s)
    ium_ma,  &
    ium_ta,  & 
    ium_matrix_condt_num,  & 
    ivm_bt,  & 
    ivm_ma,  &
    ivm_ta,  & 
    ivm_matrix_condt_num,  & 
    zt,      & 
    sfc,     & 
    ztscr01, & 
    ztscr02, & 
    ztscr03, & 
    ztscr04, & 
    ztscr05, & 
    ztscr06, & 
    l_stats_samp

use stats_type, only:  &
    stat_begin_update,  & ! Procedure(s)
    stat_end_update,  &
    stat_end_update_pt,  &
    stat_update_var_pt

use constants, only:  & 
    fstderr ! Variable(s)

use stats_precision, only:  & 
    time_precision ! Variable(s)

implicit none

! Constant parameters
integer, parameter :: &
  kp1_tdiag = 1,    & ! Thermodynamic superdiagonal index.
  k_tdiag   = 2,    & ! Thermodynamic main diagonal index.
  km1_tdiag = 3       ! Thermodynamic subdiagonal index.

! Input Variables
character(len=*), intent(in) :: & 
  solve_type     ! Desc. of what is being solved for

real(kind=time_precision), intent(in) ::  & 
  dt             ! Model timestep                            [s]

real, intent(in) ::  & 
  xpwp_sfc       ! x'w' at the surface                       [units vary]

real, dimension(gr%nnzp), intent(in) ::  & 
  xm_tndcy,    & ! The explicit time-tendency acting on xm [units vary]
  wm_zt,       & ! w wind component on thermodynamic levels  [m/s]
  Kh_zm          ! Eddy diffusivity on momentum levels       [m^2/s]

logical, intent(in) ::  &
  l_implemented  ! Flag for CLUBB being implemented in a larger model.

! Input/Output Variables
real, dimension(gr%nnzp), intent(inout) ::  & 
  xm     ! xm array on the thermodynamic grid              [units vary]

! Output Variables
real, dimension(gr%nnzp), intent(inout) ::  & 
  xpwp   ! Momentum flux                                   [units vary]

integer, intent(out) :: & 
  err_code ! clubb_singular_matrix when matrix is singular

! Local Variables
real, dimension(3,gr%nnzp) :: &
  lhs    ! Implicit contributions to xm (tridiagonal matrix)

real, dimension(gr%nnzp) :: &
  rhs    ! Right-hand side (explicit) contributions.

integer :: k, kp1, km1 ! Indices

integer :: ixm_bt, ixm_ma, ixm_ta, ixm_matrix_condt_num

real :: rcond ! Estimate of the reciprocal of the condition number on the LHS matrix


select case ( trim( solve_type ) )
case ( "um" )
  ixm_bt = ium_bt
  ixm_ma = ium_ma
  ixm_ta = ium_ta

  ! This is a diagnostic from inverting the matrix, not a budget
  ixm_matrix_condt_num = ium_matrix_condt_num
case ( "vm" )
  ixm_bt = ivm_bt
  ixm_ma = ivm_ma
  ixm_ta = ivm_ta

  ! This is a diagnostic from inverting the matrix, not a budget
  ixm_matrix_condt_num = ivm_matrix_condt_num
case default  ! Eddy scalars
  ixm_bt = 0
  ixm_ma = 0
  ixm_ta = 0

  ! This is a diagnostic from inverting the matrix, not a budget
  ixm_matrix_condt_num = 0

end select

if ( l_stats_samp ) then
 
   ! xm total time tendency (1st calculation)
   call stat_begin_update( ixm_bt, real( xm / dt ), zt )

endif
 

! Prepare tridiagonal system

! Compute the implicit portion of the xm equation.
! Build the left-hand side matrix.
call compute_um_edsclrm_lhs( solve_type, dt, wm_zt, Kh_zm,  &
                             l_implemented, lhs )

! Compute the explicit portion of the xm equation.
! Build the right-hand side vector.
call compute_um_edsclrm_rhs( solve_type, dt, Kh_zm, xm,  &
                             xm_tndcy, xpwp_sfc, rhs )


! Store momentum flux (explicit component)
! The surface flux, x'w'(1) = x'w'|_sfc, is set elsewhere in the model.
xpwp(1) = xpwp_sfc
! Solve for x'w' at all intermediate model levels.
! A Crank-Nicholson timestep is used.
do k = 2, gr%nnzp-1, 1
   xpwp(k) = -0.5 * Kh_zm(k) * gr%dzm(k) * ( xm(k+1) - xm(k) ) 
enddo
! A zero-flux boundary condition at the top of the model, d(xm)/dz = 0, means 
! that u'w' at the top model level is 0, since x'w' = - K_zm * d(xm)/dz.
xpwp(gr%nnzp) = 0.


! Solve tridiagonal system for xm.
if ( l_stats_samp .and. ixm_matrix_condt_num > 0 ) then
  call tridag_solvex & 
       ( solve_type, gr%nnzp, 1, &                                  ! Intent(in) 
         lhs(kp1_tdiag,:), lhs(k_tdiag,:), lhs(km1_tdiag,:), rhs, & ! Intent(inout)
         xm, rcond, err_code )                                      ! Intent(out)

  ! Est. of the condition number of the variance LHS matrix 
  call stat_update_var_pt( ixm_matrix_condt_num, 1, 1.0 / rcond, &  ! Intent(in)
                           sfc )                      ! Intent(inout)
else
 
  call tridag_solve( solve_type, gr%nnzp, 1, & ! In
                     lhs(kp1_tdiag,:),  lhs(k_tdiag,:), lhs(km1_tdiag,:), rhs, & ! Inout
                     xm, err_code ) ! Out
end if


! Second part of momentum (implicit component)
! Solve for x'w' at all intermediate model levels.
! A Crank-Nicholson timestep is used.
do k = 2, gr%nnzp-1, 1
   xpwp(k) = xpwp(k)  &
              - 0.5 * Kh_zm(k) * gr%dzm(k) * ( xm(k+1) - xm(k) )
enddo


if ( l_stats_samp ) then

   ! Finalize implicit contributions for xm

   do k = 2, gr%nnzp-1, 1

      km1 = max( k-1, 1 )
      kp1 = min( k+1, gr%nnzp )

      ! xm mean advection
      call stat_update_var_pt( ixm_ma, k, &
             ztscr01(k) * xm(km1) &
           + ztscr02(k) * xm(k) &
           + ztscr03(k) * xm(kp1), zt )

      ! xm turbulent transport (implicit component)
      call stat_end_update_pt( ixm_ta, k, &
             ztscr04(k) * xm(km1) &
           + ztscr05(k) * xm(k) &
           + ztscr06(k) * xm(kp1), zt )

    enddo

    ! Upper boundary condition:
    ! xm turbulent transport (implicit component)
    k   = gr%nnzp
    km1 = max( k-1, 1 )
    call stat_end_update_pt( ixm_ta, k, &
           ztscr04(k) * xm(km1) &
         + ztscr05(k) * xm(k), zt )

endif

if ( l_stats_samp ) then

   ! xm total time tendency (2nd calculation)
   call stat_end_update( ixm_bt, real( xm / dt ), zt )

endif


return
end subroutine compute_um_edsclrm_solve

!===============================================================================
subroutine compute_uv_tndcy( solve_type, fcor, perp_wind_m, perp_wind_g,  &
                             l_implemented, xm_tndcy )
!
! Description:
! Computes the explicit tendency for the um and vm wind components.
!
! The only explicit tendency that is involved in the d(um)/dt or d(vm)/dt 
! equations is the Coriolis tendency.
!
! The d(um)/dt equation contains the term:
!
! - f * ( v_g - vm );
!
! where f is the Coriolis parameter and v_g is the v component of the 
! geostrophic wind.
!
! Likewise, the d(vm)/dt equation contains the term:
!
! + f * ( u_g - um );
!
! where u_g is the u component of the geostrophic wind.
!
! This term is treated completely explicitly.  The values of um, vm, u_g, and 
! v_g are all found on the thermodynamic levels.
!
! References:
!-----------------------------------------------------------------------

use grid_class, only: & 
    gr

use grid_class, only: & 
    zt2zm, & 
    ddzm         

use stats_type, only: & 
    stat_update_var

use stats_variables, only:      &
    ium_gf, & 
    ium_cf, & 
    ivm_gf, & 
    ivm_cf, & 
    zt, & 
    l_stats_samp

implicit none

! Input Variables
character(len=*), intent(in) ::  &
  solve_type      ! Description of what is being solved for

real, intent(in) ::  & 
  fcor            ! Coriolis parameter     [s^-1]

real, dimension(gr%nnzp), intent(in) :: & 
  perp_wind_m,  & ! Perpendicular component of the mean wind (e.g. v, for the u-eqn) [m/s]
  perp_wind_g     ! Perpendicular component of the geostropic wind (e.g. vg)         [m/s]

logical, intent(in) :: & 
  l_implemented   ! Flag for CLUBB being implemented in a larger model.

! Output Variables
real, dimension(gr%nnzp), intent(out) ::  &
  xm_tndcy        ! xm tendency            [m/s^2]

! Local Variables
integer :: & 
  ixm_gf, & 
  ixm_cf

real, dimension(gr%nnzp) :: & 
  xm_gf, & 
  xm_cf


if ( .not. l_implemented ) then
  ! Only compute the Coriolis term if the model is running on it's own, 
  ! and is not part of a larger, host model.

  select case (trim(solve_type))

  case("um")
 
    ixm_gf = ium_gf
    ixm_cf = ium_cf

    xm_gf = - fcor * perp_wind_g

    xm_cf = fcor * perp_wind_m

  case("vm") 
 
    ixm_gf = ivm_gf
    ixm_cf = ivm_cf

    xm_gf = fcor * perp_wind_g

    xm_cf = -fcor * perp_wind_m

  case default

    ixm_gf = 0
    ixm_cf = 0

    xm_gf = 0.

    xm_cf = 0.

  end select

  xm_tndcy = xm_gf + xm_cf 

  if ( l_stats_samp ) then
 
    call stat_update_var( ixm_gf, xm_gf, zt )
  
    call stat_update_var( ixm_cf, xm_cf, zt )

  endif

else   ! implemented in a host model.

  xm_tndcy = 0.0

endif


return
end subroutine compute_uv_tndcy

!===============================================================================
subroutine compute_um_edsclrm_lhs( solve_type, dt, wm_zt, Kh_zm,  &
                                   l_implemented, lhs )

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
  solve_type    ! Description of what is being solved for

real(kind=time_precision), intent(in) :: & 
  dt            ! Model timestep                           [s]

real, dimension(gr%nnzp), intent(in) :: &
  wm_zt,      & ! w wind component on thermodynamic levels [m/s]
  Kh_zm         ! Eddy diffusivity on momentum levels      [m^2/s]

logical, intent(in) ::  & 
  l_implemented ! Flag for CLUBB being implemented in a larger model.

! Output Variable
real, dimension(3,gr%nnzp), intent(out) :: &
  lhs           ! Implicit contributions to xm (tridiagonal matrix)

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
   if ( .not. l_implemented ) then

      lhs(kp1_tdiag:km1_tdiag,k)  &
      = lhs(kp1_tdiag:km1_tdiag,k)  &
      + term_ma_zt_lhs( wm_zt(k), gr%dzt(k), k )

   else

      lhs(kp1_tdiag:km1_tdiag,k)  &
      = lhs(kp1_tdiag:km1_tdiag,k) + 0.0

   endif

   ! LHS turbulent advection term (solved as an eddy-diffusion term).
   if ( k == 2 ) then
      ! The lower boundary condition needs to be applied here at level 2.
      ! The lower boundary condition is a "fixed flux" boundary condition.  
      ! The coding is the same as for a zero-flux boundary condition, but with 
      ! an extra term added on the right-hand side at the boundary level.  For
      ! the rest of the model code, a zero-flux boundary condition is applied
      ! at level 1, and thus subroutine diffusion_zt_lhs is set-up to do that.
      ! In order to apply the same boundary condition code here at level 2, and
      ! adjuster needs to be used to tell diffusion_zt_lhs to use the code at 
      ! level 2 that it normally uses at level 1.
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
         if ( .not. l_implemented ) then
            tmp(1:3) &
            = term_ma_zt_lhs( wm_zt(k), gr%dzt(k), k )
            ztscr01(k) = -tmp(3)
            ztscr02(k) = -tmp(2)
            ztscr03(k) = -tmp(1)
         else
            ztscr01(k) = 0.0
            ztscr02(k) = 0.0
            ztscr03(k) = 0.0
         endif
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

! The lower boundary condition is a fixed-flux boundary condition, which gets 
! added into the time-tendency equation at level 2.
! The value of xm(1) is solved for based on the result of xm(2) after xm(2) has
! been advanced to the next timestep.  The new result for xm(1) is based on the
! relationship:  x'w'|_sfc = - K_zm(1) * dzm(1) * ( xm(2) - xm(1) ); where 
! x'w'|_sfc is the given value of the surface flux, and K_zm(1), dzm(1), and 
! xm(2) are known values.  The implicit portion of the equation is:
! [+1] * xm(1,<t+1>)  +  [-1] * xm(2,<t+1>); while the explicit portion of the 
! equation is:  + x'w'|_sfc / ( K_zm(1) * dzm(1) ).
lhs(k_tdiag,1)   = 1.0
lhs(kp1_tdiag,1) = -1.0

! Upper Boundary

! A zero-flux boundary condition is used at the upper boundary, meaning that xm 
! is not allowed to exit the model through the upper boundary.  This boundary
! condition is invoked by calling diffusion_zt_lhs at the uppermost level.
k   = gr%nnzp
km1 = max( k-1, 1 )

! LHS turbulent advection term (solved as an eddy-diffusion term) 
! at the upper boundary.
lhs(kp1_tdiag:km1_tdiag,k)  &
= lhs(kp1_tdiag:km1_tdiag,k)  &
+ (1.0/2.0)  &
* diffusion_zt_lhs( Kh_zm(k), Kh_zm(km1), 0.0,  & 
                    gr%dzm(km1), gr%dzm(k), gr%dzt(k), k )

! LHS time tendency term at the upper boundary.
lhs(k_tdiag,k) = real( lhs(k_tdiag,k) + ( 1.0 / dt ) )

if ( l_stats_samp ) then

   ! Statistics:  implicit contributions for um or vm.

   if ( ixm_ta > 0 ) then
      tmp(1:3)  &
      = (1.0/2.0)  &
      * diffusion_zt_lhs( Kh_zm(k), Kh_zm(km1), 0.0,  &
                          gr%dzm(km1), gr%dzm(k), gr%dzt(k), k )
      ztscr04(k) = -tmp(3)
      ztscr05(k) = -tmp(2)
      ztscr06(k) = -tmp(1)
   endif

endif  ! lstats_samp


return
end subroutine compute_um_edsclrm_lhs

!===============================================================================
subroutine compute_um_edsclrm_rhs( solve_type, dt, Kh_zm, xm,  &
                                   xm_tndcy, xpwp_sfc, rhs )

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
  solve_type ! Description of what is being solved for

real(kind=time_precision), intent(in) :: & 
  dt           ! Model timestep                                   [s]

real, dimension(gr%nnzp), intent(in) :: &
  Kh_zm,     & ! Eddy diffusivity on momentum levels              [m^2/s]
  xm,        & ! Eddy-scalar variable, xm (thermodynamic levels)  [units vary]
  xm_tndcy     ! The explicit time-tendency acting on xm          [units vary]

real, intent(in) :: &
  xpwp_sfc     ! x'w' at the surface                              [units vary]

! Output Variable
real, dimension(gr%nnzp), intent(out) :: &
  rhs          ! Right-hand side (explicit) contributions.

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
      ! The lower boundary condition is a "fixed flux" boundary condition.  
      ! The coding is the same as for a zero-flux boundary condition, but with 
      ! an extra term added on the right-hand side at the boundary level.  For
      ! the rest of the model code, a zero-flux boundary condition is applied
      ! at level 1, and thus subroutine diffusion_zt_lhs is set-up to do that.
      ! In order to apply the same boundary condition code here at level 2, and
      ! adjuster needs to be used to tell diffusion_zt_lhs to use the code at 
      ! level 2 that it normally uses at level 1.
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
   rhs(k) = rhs(k) + xm_tndcy(k)

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


! Boundary Conditions

! Lower Boundary

! The lower boundary condition is a fixed-flux boundary condition, which gets 
! added into the time-tendency equation at level 2.
rhs(2) = rhs(2) + gr%dzt(2) * xpwp_sfc
! The value of xm(1) is solved for based on the result of xm(2) after xm(2) has
! been advanced to the next timestep.  The new result for xm(1) is based on the
! relationship:  x'w'|_sfc = - K_zm(1) * dzm(1) * ( xm(2) - xm(1) ); where 
! x'w'|_sfc is the given value of the surface flux, and K_zm(1), dzm(1), and 
! xm(2) are known values.  The implicit portion of the equation is:
! [+1] * xm(1,<t+1>)  +  [-1] * xm(2,<t+1>); while the explicit portion of the 
! equation is:  + x'w'|_sfc / ( K_zm(1) * dzm(1) ).
rhs(1) = + xpwp_sfc / ( Kh_zm(1) * gr%dzm(1) ) 

if ( l_stats_samp ) then

   ! Statistics:  explicit contributions for um or vm.

   if ( ixm_ta > 0 ) then
      call stat_modify_pt( ixm_ta, 2, &
           + gr%dzt(2) * xpwp_sfc, zt )
   endif

endif  ! lstats_samp


! Upper Boundary

! A zero-flux boundary condition is used at the upper boundary, meaning that xm 
! is not allowed to exit the model through the upper boundary.  This boundary
! condition is invoked by calling diffusion_zt_lhs at the uppermost level.
k   = gr%nnzp
km1 = max( k-1, 1 )

! RHS turbulent advection term (solved as an eddy-diffusion term) 
! at the upper boundary.
rhs_diff(1:3)  &
= (1.0/2.0)  &
* diffusion_zt_lhs( Kh_zm(k), Kh_zm(km1), 0.0,  &
                    gr%dzm(km1), gr%dzm(k), gr%dzt(k), k )
rhs(k)   =   rhs(k) &
           - rhs_diff(3) * xm(km1) &
           - rhs_diff(2) * xm(k)

! RHS time tendency term at the upper boundary.
rhs(k) = real( rhs(k) + ( 1.0 / dt ) * xm(k) )

if ( l_stats_samp ) then

   ! Statistics:  explicit contributions for um or vm.

   if ( ixm_ta > 0 ) then
      call stat_begin_update_pt( ixm_ta, k, &
             rhs_diff(3) * xm(km1) &
           + rhs_diff(2) * xm(k), zt )
   endif

endif  ! lstats_samp


return
end subroutine compute_um_edsclrm_rhs

!===============================================================================

end module advance_windm_edsclrm_module
