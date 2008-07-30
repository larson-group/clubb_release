!------------------------------------------------------------------------
! $Id: compute_um_edsclrm_mod.F90,v 1.5 2008-07-30 19:17:35 dschanen Exp $
!------------------------------------------------------------------------
module compute_um_edsclrm_mod

implicit none

private ! Set Default Scope

public :: compute_um_edsclrm, compute_uv_tndcy

contains

!------------------------------------------------------------------------
subroutine compute_um_edsclrm( solve_type, xpwp_sfc, &
                               xm_tndcy,  & 
                               Khm, dt,   &
                               xm, xpwp,  &
                               err_code )
!       Description:
!         Prognoses a horizontal wind component or other 
!         eddy diffusivity variable.
!         Diagnoses the turbulent flux.
!         Uses a Crank-Nicholson time-stepping algorithm.

!       References:
!       Eqn. 8 & 9 on p. 3545 of 
!       ``A PDF-Based Model for Boundary Layer Clouds. Part I:
!         Method and Model Description'' Golaz, et al. (2002)
!       JAS, Vol. 59, pp. 3540--3551.
!------------------------------------------------------------------------

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
    lstats_samp

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
  Khm       ! Diffusion coefficient on momentum grid  [m^2/s]

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

 
if ( lstats_samp ) then
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
ctmp = -0.5 * Khm(2) * gr%dzm(2) * gr%dzt(2)

a(2)   = atmp
c(2)   = ctmp
b(2)   = real( - ctmp + 1./dt )
rhs(2) = real( ( ctmp + 1./dt ) * xm(2) & 
         - ctmp * xm(3) & 
         + xm_tndcy(2) & 
         + xpwp_sfc * gr%dzt(2) )

 
   if ( lstats_samp .and. ixm_ta > 0 ) then
     ztscr01(1) = 0.0
     ztscr02(1) = 0.0
     ztscr03(1) = 0.0

     ztscr01(2) = -atmp
     ztscr02(2) = ctmp
     ztscr03(2) = -ctmp
   end if
! Loop from level 3 to gr%nnzp

do k=3, gr%nnzp-1, 1

   atmp = -0.5 * Khm(k-1) * gr%dzt(k) * gr%dzm(k-1)
   ctmp = -0.5 * Khm(k) * gr%dzt(k) * gr%dzm(k)

   a(k)   = atmp
   c(k)   = ctmp
   b(k)   = real( - atmp - ctmp + 1./dt )
   rhs(k) = real( - atmp * xm(k-1) & 
            + ( atmp + ctmp + 1./dt ) * xm(k) & 
            - ctmp * xm(k+1) & 
            + xm_tndcy(k) )
 
   if ( lstats_samp .and. ixm_ta > 0 ) then
     ztscr01(k) = -atmp
     ztscr02(k) =  atmp + ctmp
     ztscr03(k) = -ctmp
   end if
end do

! Level gr%nnzp. We impose zero flux from model top

! Vince Larson simplified the upper BC.  It still imposes
!    zero flux. 7 Jul 2007 
!        atmp = -0.5 * Khm(gr%nnzp-1) * gr%dzm(gr%nnzp) 
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
   if ( lstats_samp .and. ixm_ta > 0 ) then
     ztscr01(gr%nnzp) = -atmp
     ztscr02(gr%nnzp) =  ctmp
     ztscr03(gr%nnzp) =  atmp
   end if

   
!    Caused problems with DYCOMS II RF02
!     .               + xm_tndcy(gr%nnzp)

!    Attempted to compensate for the DYCOMS problem using the code below.
!    Doesn't actually seem to make a difference
!       atmp = -0.5 * Khm(gr%nnzp-1) * gr%dzt(gr%nnzp) * gr%dzm(gr%nnzp-1)

!       a(gr%nnzp)   = atmp
!       b(gr%nnzp)   = - c(gr%nnzp-1) + 1./dt
!       c(gr%nnzp)   = UNDEFINED
!       rhs(gr%nnzp) = - atmp * x(gr%nnzp-1)
!    .                 + ( atmp + 1./dt ) * x(gr%nnzp)
!    .                 + xm_tndcy(gr%nnzp)

! Store momentum flux (explicit component)
xpwp(1) = xpwp_sfc
do k=2,gr%nnzp-1
  xpwp(k) = -0.5 * Khm(k) * gr%dzm(k) * ( xm(k+1) - xm(k) ) 
end do
xpwp(gr%nnzp) = 0.

 
   ! Turbulent transport (explicit component)
   if ( lstats_samp .and. ixm_ta > 0 ) then
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
    write(fstderr,*) "Khm = ", Khm
   
    write(fstderr,*) "Intent(inout)"
   
    write(fstderr,*) "xm = ", xm
   
    write(fstderr,*) "Intent(out)"
   
    write(fstderr,*) "xpwp = ", xpwp
             
  return
end if

! Second part of momentum (implicit component)
do k=2, gr%nnzp-1, 1
  xpwp(k) = xpwp(k)  & 
             - 0.5 * Khm(k) * gr%dzm(k) * ( xm(k+1) - xm(k) )
end do

 
if ( lstats_samp ) then
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
subroutine compute_uv_tndcy( solve_type, xm, wmt, fcor, perp_wind_m, perp_wind_g, implemented, & 
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
    lstats_samp

  implicit none

  ! Input Variables
  character(len=*), intent(in) :: solve_type

  real, dimension(gr%nnzp), intent(in) ::  & 
    xm,  & ! u/v wind                                       [m/s]   
    wmt ! wm on thermodynaming grid                     [m/s]

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

    xm_ma = - wmt * ddzm( zt2zm( xm ) )

    xm_gf = - fcor * perp_wind_g

    xm_cf = fcor * perp_wind_m

  case("vm") 
 
    ixm_ma = ivm_ma
    ixm_gf = ivm_gf
    ixm_cf = ivm_cf

    xm_ma = - wmt * ddzm( zt2zm( xm ) )

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

 
  if ( lstats_samp ) then
    call stat_update_var( ixm_ma, xm_ma, zt )

    call stat_update_var( ixm_gf, xm_gf, zt )
  
    call stat_update_var( ixm_cf, xm_cf, zt )
  endif
else
  xmt = 0.0
endif

end subroutine compute_uv_tndcy


end module compute_um_edsclrm_mod
