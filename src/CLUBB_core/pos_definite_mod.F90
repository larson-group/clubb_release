!$Id$
module pos_definite_mod

  implicit none

  public :: pos_definite_adj

  private ! Default Scope

  contains
!-----------------------------------------------------------------------
  subroutine pos_definite_adj & 
            ( dt, field_grid, field_np1,  & 
              flux_np1, field_n, field_pd, flux_pd )
!       Description:
!       Applies a  flux conservative positive definite scheme to a variable

!       There are two possible grids:       
!       (1) flux on zm  field on zt
!       then
!       flux_zt(k) = ( flux_zm(k) + flux_zm(k-1) ) / 2

!             CLUBB grid                  Smolarkiewicz grid      
!       m +-- flux  zm(k)  --+               flux        k + 1/2
!       t +-- field zt(k)  --+               field, fout k
!       m +-- flux  zm(k-1) --+              flux        k - 1/2
!       t +-- field zt(k-1) --+

!       (1) flux on zt field on zm
!       then
!       flux_zm(k) = ( flux_zt(k) + flux_zt(k+1) ) / 2

!             CLUBB grid                  Smolarkiewicz grid      
!       m +-- field  (k+1)  --+               
!       t +-- flux   (k+1)  --+               flux        k + 1/2
!       m +-- field  (k)    --+               field, fout k 
!       t +-- flux   (k)    --+               flux        k - 1/2


!       References:
!       ``A Positive Definite Advection Scheme Obtained by
!       Nonlinear Renormalization of the Advective Fluxes'' Smolarkiewicz (1989)
!       Monthly Weather Review, Vol. 117, pp. 2626--2632
!-----------------------------------------------------------------------

  use grid_class, only: & 
      gr, & ! Variable(s)
      ddzt,  & ! Function
      ddzm  ! Function

  use constants, only :  & 
      eps, & ! Variable(s)
      zero_threshold

  use stats_precision, only:  & 
      time_precision ! Variable(s)

  implicit none

  ! Constant parameters
!       real, parameter :: beta = 1.0 ! Coefficient of 'Y' (i.e. field_np1)

  ! External
  intrinsic :: eoshift, kind, any, min, max

  ! Input variables
  real(kind=time_precision), intent(in) :: & 
  dt ! Timestep    [s]

  character(len=2), intent(in) :: & 
  field_grid ! The grid of the field, either zt or zm

  real, dimension(gr%nnzp), intent(in) ::  & 
  field_n ! The field (e.g. rtm) at n, prior to n+1

  real, dimension(gr%nnzp), intent(out) ::  & 
  flux_pd,  & ! Budget of the change in the flux term due to the scheme
  field_pd ! Budget of the change in the mean term due to the scheme

  ! Output Variables

  real, intent(inout), dimension(gr%nnzp) :: & 
  field_np1,   & ! Field at n+1 (e.g. rtm in [kg/kg])
  flux_np1    ! Flux applied to field

  ! Local Variables
  integer :: & 
  kabove,   & ! # of vertical levels the flux higher point resides
  kbelow   ! # of vertical levels the flux lower point resides

!       real :: alpha ! Coefficient for matrix 'A'

  integer ::  & 
  k, kmhalf, kp1, kphalf ! Loop indices

  real, dimension(gr%nnzp) :: & 
  flux_plus, flux_minus, & ! [F_i+1/2]^+ [F_i+1/2]^- in Smolarkiewicz 
  fout,                  & ! (A4) F_i{}^OUT, or the sum flux_plus+flux_minus
  flux_lim,              & ! Correction applied to flux at n+1
  field_nonlim          ! Temporary variable for calculation

  real, dimension(gr%nnzp) ::  & 
  dz_over_dt ! Conversion factor  [m/s]


!-----------------------------------------------------------------------

  ! If all the values are positive or the values at the previous 
  ! timestep were negative, then just return
  if ( .not. any( field_np1 < 0. ) .or. any( field_n < 0. ) ) then
    flux_pd  = 0.
    field_pd = 0.
    return
  end if

  if ( field_grid == "zm" ) then
    kabove = 0
    kbelow = 1
  else if ( field_grid == "zt" ) then
    kabove = 1
    kbelow = 0
  else
    ! This is only necessary to avoid a compiler warning in g95
    kabove = -1
    kbelow = -1
    ! Joshua Fasching June 2008

    stop "Error in pos_def_adj"
  end if

  print *, "Correcting flux"

  do k = 1, gr%nnzp, 1

    ! Def. of F+ and F- from eqn 2 Smolarkowicz
    flux_plus(k)  =  max( zero_threshold, flux_np1(k) ) ! defined on flux levels
    flux_minus(k) = -min( zero_threshold, flux_np1(k) ) ! defined on flux levels

    if ( field_grid == "zm" ) then
      dz_over_dt(k) = real( ( 1./gr%dzm(k) ) / dt )

    else if ( field_grid == "zt" ) then
      dz_over_dt(k) = real( ( 1./gr%dzt(k) ) / dt )

    end if

  end do 

  do k = 1, gr%nnzp, 1
    ! If the scalar variable is on the kth t-level, then
    ! Smolarkowicz's k+1/2 flux level is the kth m-level in CLUBB.

    ! If the scalar variable is on the kth m-level, then
    ! Smolarkowicz's k+1/2 flux level is the k+1 t-level in CLUBB.

    kphalf = min( k+kabove, gr%nnzp ) ! k+1/2 flux level
    kmhalf = max( k-kbelow, 1 )       ! k-1/2 flux level

    ! Eqn A4 from Smolarkowicz
    ! We place a limiter of eps to prevent a divide by zero, and
    !   after this calculation fout is on the scalar level, and
    !   fout is the total outward flux for the scalar level k.

    fout(k) = max( flux_plus(kphalf) + flux_minus(kmhalf), eps )

  end do 


  do k = 1, gr%nnzp, 1
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! FIXME:
    ! We haven't tested this for negative values at the gr%nnzp level
    ! -dschanen 13 June 2008
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    kphalf = min( k+kabove, gr%nnzp ) ! k+1/2 flux level
    kp1    = min( k+1, gr%nnzp )      ! k+1 scalar level

    ! Eqn 10 from Smolarkowicz (1989)

    flux_lim(kphalf) & 
    = max( min( flux_np1(kphalf), & 
                ( flux_plus(kphalf)/fout(k) ) * field_n(k) & 
                  * dz_over_dt(k) & 
              ), & 
           -( ( flux_minus(kphalf)/fout(kp1) ) * field_n(kp1) & 
                * dz_over_dt(k) ) & 
         )
  end do

  ! Boundary conditions
  flux_lim(1) = flux_np1(1) 
  flux_lim(gr%nnzp) = flux_np1(gr%nnzp)

  flux_pd = real( ( flux_lim - flux_np1 ) / dt )

  field_nonlim = field_np1

  !alpha = -dt

  ! Apply change to field at n+1
  if ( field_grid == "zt" ) then

    field_np1 =  & 
        real( -dt * ddzm( flux_lim - flux_np1 ) + field_np1 )

  else if ( field_grid == "zm" ) then

    field_np1 =  & 
        real( -dt * ddzt( flux_lim - flux_np1 ) + field_np1 )

  end if

  ! Removed in favor of calling ddzt and ddzm.  It turns out this
  ! isn't quite generalized for use with either grid.
!       call band_mult( 'N', gr%nnzp, gr%nnzp, 1, 1, 1, 1,
!    .                  alpha, beta, flux_lhs, flux_pd, 
!    .                  field_np1 ) 

  ! Determine the total time tendency in field due to this calculation
  ! (for diagnostic purposes)
  field_pd = real( ( field_np1 - field_nonlim ) / dt )

  ! Replace the non-limited flux with the limited flux
  flux_np1 = flux_lim

  return
  end subroutine pos_definite_adj

  end module pos_definite_mod
