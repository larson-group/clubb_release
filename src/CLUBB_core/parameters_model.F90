!-------------------------------------------------------------------------------
! $Id$
!===============================================================================
module parameters_model

! Description:
! Contains model parameters that are determined at run time rather than
! compile time.
!
! References:
! None
!-------------------------------------------------------------------------------
  implicit none

  private ! Default scope

  ! Maximum allowable value for Lscale [m].
  ! Value depends on whether the model is run by itself or as part of a
  ! host model.
  real, public :: Lscale_max

!$omp threadprivate(Lscale_max)

  ! Maximum magnitude of PDF parameter 'mixt_frac'. 
  real, public :: mixt_frac_max_mag

!$omp threadprivate(mixt_frac_max_mag)

  ! Model parameters and constraints setup in the namelists
  real, public ::  & 
    T0,       & ! Reference temperature (usually 300)  [K]
    ts_nudge    ! Timescale of u/v nudging             [s]

#ifdef GFDL
 real, public ::  &   ! h1g, 2010-06-15
    cloud_frac_min    ! minimum cloud fraction for droplet #
#endif


!$omp threadprivate(T0, ts_nudge)
  integer, public :: & 
    sclr_dim,        & ! Number of passive scalars
    edsclr_dim,      & ! Number of eddy-diff. passive scalars
    hydromet_dim       ! Number of hydrometeor species

!$omp threadprivate(sclr_dim, edsclr_dim, hydromet_dim)

  real, dimension(:), allocatable, public :: & 
    sclr_tol ! Threshold(s) on the passive scalars  [units vary]

!$omp threadprivate(sclr_tol)

  public :: setup_parameters_model 

  contains

!-------------------------------------------------------------------------------
  subroutine setup_parameters_model &
             ( T0_in, ts_nudge_in, &
               hydromet_dim_in, & 
               sclr_dim_in, sclr_tol_in, edsclr_dim_in, &
               Lscale_max_in &

#ifdef GFDL
              , cloud_frac_min_in &    ! hlg, 2010-6-15
#endif
    
              )

! Description:
!   Sets parameters to their initial values
!
! References:
!   None
!-------------------------------------------------------------------------------
    use constants_clubb, only: Skw_max_mag, Skw_max_mag_sqd

    implicit none

    ! External
    intrinsic :: sqrt, allocated

    ! Input Variables
    real, intent(in) ::  & 
      T0_in,        & ! Ref. temperature             [K]
      ts_nudge_in,  & ! Timescale for u/v nudging    [s]
      Lscale_max_in   ! Largest value for Lscale     [m]

#ifdef GFDL
    real, intent(in) ::  cloud_frac_min_in  ! h1g, 2010-06-15
#endif


    integer, intent(in) :: & 
      hydromet_dim_in,  & ! Number of hydrometeor species
      sclr_dim_in,      & ! Number of passive scalars
      edsclr_dim_in       ! Number of eddy-diff. passive scalars

    real, intent(in), dimension(sclr_dim_in) :: & 
      sclr_tol_in     ! Threshold on passive scalars

    ! --- Begin Code --- !
     
    ! Formula from subroutine pdf_closure, where sigma_sqd_w = 0.4 and Skw =
    ! Skw_max_mag in this formula.  Note that this is constant, but can't appear
    ! with a Fortran parameter attribute, so we define it here. 
    mixt_frac_max_mag = 1.0 &
      - ( 0.5 * ( 1.0 - Skw_max_mag / sqrt( 4.0 * ( 1.0 - 0.4 )**3 + Skw_max_mag_sqd ) ) )

    Lscale_max = Lscale_max_in

    T0       = T0_in
    ts_nudge = ts_nudge_in

    hydromet_dim = hydromet_dim_in
    sclr_dim     = sclr_dim_in
    edsclr_dim   = edsclr_dim_in

    ! In a tuning run, this array has the potential to be allocated already
    if ( .not. allocated( sclr_tol ) ) then
      allocate( sclr_tol(1:sclr_dim) )
    else
      deallocate( sclr_tol )
      allocate( sclr_tol(1:sclr_dim) )
    end if

    sclr_tol(1:sclr_dim) = sclr_tol_in(1:sclr_dim)


#ifdef GFDL
     cloud_frac_min = cloud_frac_min_in  ! h1g, 2010-06-15
#endif


    return
  end subroutine setup_parameters_model
!-------------------------------------------------------------------------------

end module parameters_model
