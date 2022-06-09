!-------------------------------------------------------------------------
!$Id$
!===============================================================================
module pos_definite_module

  implicit none

  public :: pos_definite_adj

  private ! Default Scope

  contains
!-----------------------------------------------------------------------
  subroutine pos_definite_adj( nz, ngrdcol, gr, dt, field_grid, & 
                               field_np1, flux_np1, field_n, &
                               field_pd, flux_pd )
! Description:
!   Applies a  flux conservative positive definite scheme to a variable

!   There are two possible grids:
!   (1) flux on zm  field on zt
!   then
!   flux_zt(k) = ( flux_zm(k) + flux_zm(k-1) ) / 2

!         CLUBB grid                  Smolarkiewicz grid
!   m +-- flux  zm(k)  --+               flux        k + 1/2
!   t +-- field zt(k)  --+               field, fout k
!   m +-- flux  zm(k-1) --+              flux        k - 1/2
!   t +-- field zt(k-1) --+

!   (1) flux on zt field on zm
!   then
!   flux_zm(k) = ( flux_zt(k) + flux_zt(k+1) ) / 2

!         CLUBB grid                  Smolarkiewicz grid
!   m +-- field  (k+1)  --+
!   t +-- flux   (k+1)  --+               flux        k + 1/2
!   m +-- field  (k)    --+               field, fout k
!   t +-- flux   (k)    --+               flux        k - 1/2


! References:
!   ``A Positive Definite Advection Scheme Obtained by
!     Nonlinear Renormalization of the Advective Fluxes'' Smolarkiewicz (1989)
!     Monthly Weather Review, Vol. 117, pp. 2626--2632
!-----------------------------------------------------------------------

    use grid_class, only: & 
        grid, & ! Type
        ddzt, & ! Function
        ddzm    ! Function

    use constants_clubb, only :  & 
        eps, & ! Variable(s)
        zero_threshold, &
        zero

    use clubb_precision, only:  & 
        core_rknd ! Variable(s)

    use error_code, only: &
        clubb_at_least_debug_level   ! Procedure

    implicit none


    ! -------------------- Input variables --------------------
    integer, intent(in) :: &
      nz, &
      ngrdcol
    
    type (grid), target, intent(in) :: gr
        
    real( kind = core_rknd ), intent(in) :: & 
      dt ! Timestep    [s]

    character(len=2), intent(in) :: & 
      field_grid ! The grid of the field, either zt or zm

    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) ::  & 
      field_n ! The field (e.g. rtm) at n, prior to n+1

    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(out) ::  & 
      flux_pd,  & ! Budget of the change in the flux term due to the scheme
      field_pd    ! Budget of the change in the mean term due to the scheme

    ! Output Variables

    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,nz) :: & 
      field_np1,   & ! Field at n+1 (e.g. rtm in [kg/kg])
      flux_np1    ! Flux applied to field

    ! Local Variables
    integer :: & 
      kabove,   & ! # of vertical levels the flux higher point resides
      kbelow   ! # of vertical levels the flux lower point resides

    integer ::  & 
      i, k, kmhalf, kp1, kphalf ! Loop indices

    real( kind = core_rknd ), dimension(ngrdcol,nz) :: & 
      flux_plus, flux_minus, & ! [F_i+1/2]^+ [F_i+1/2]^- in Smolarkiewicz 
      fout,                  & ! (A4) F_i{}^OUT, or the sum flux_plus+flux_minus
      flux_lim,              & ! Correction applied to flux at n+1
      field_nonlim          ! Temporary variable for calculation

    real( kind = core_rknd ), dimension(ngrdcol,nz) ::  & 
      dz_over_dt ! Conversion factor  [m/s]


!-----------------------------------------------------------------------

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

      error stop "Error in pos_def_adj"
    end if

    if ( clubb_at_least_debug_level( 1 ) ) then
      print *, "Correcting flux"
    end if

    do k = 1, nz, 1
      do i = 1, ngrdcol

        ! Def. of F+ and F- from eqn 2 Smolarkowicz
        flux_plus(i,k)  =  max( zero_threshold, flux_np1(i,k) ) ! defined on flux levels
        flux_minus(i,k) = -min( zero_threshold, flux_np1(i,k) ) ! defined on flux levels

        if ( field_grid == "zm" ) then
          dz_over_dt(i,k) = ( 1._core_rknd/gr%invrs_dzm(i,k) ) / dt

        else if ( field_grid == "zt" ) then
          dz_over_dt(i,k) = ( 1._core_rknd/gr%invrs_dzt(i,k) ) / dt

        end if
        
      end do
    end do

    do k = 1, nz, 1
      do i = 1, ngrdcol
        ! If the scalar variable is on the kth t-level, then
        ! Smolarkowicz's k+1/2 flux level is the kth m-level in CLUBB.

        ! If the scalar variable is on the kth m-level, then
        ! Smolarkowicz's k+1/2 flux level is the k+1 t-level in CLUBB.

        kphalf = min( k+kabove, nz ) ! k+1/2 flux level
        kmhalf = max( k-kbelow, 1 )       ! k-1/2 flux level

        ! Eqn A4 from Smolarkowicz
        ! We place a limiter of eps to prevent a divide by zero, and
        !   after this calculation fout is on the scalar level, and
        !   fout is the total outward flux for the scalar level k.

        fout(i,k) = max( flux_plus(i,kphalf) + flux_minus(i,kmhalf), eps )
      end do
    end do


    do k = 1, nz, 1
      do i = 1, ngrdcol
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! FIXME:
        ! We haven't tested this for negative values at the nz level
        ! -dschanen 13 June 2008
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        kphalf = min( k+kabove, nz ) ! k+1/2 flux level
        kp1    = min( k+1, nz )      ! k+1 scalar level

        ! Eqn 10 from Smolarkowicz (1989)

        flux_lim(i,kphalf) & 
        = max( min( flux_np1(i,kphalf), & 
                    ( flux_plus(i,kphalf)/fout(i,k) ) * field_n(i,k) & 
                      * dz_over_dt(i,k) & 
                  ), & 
               -( ( flux_minus(i,kphalf)/fout(i,kp1) ) * field_n(i,kp1) & 
                    * dz_over_dt(i,k) ) & 
             )
      end do
    end do

    ! Boundary conditions
    flux_lim(:,1) = flux_np1(:,1)
    flux_lim(:,nz) = flux_np1(:,nz)

    do i = 1, ngrdcol
      ! Only set flux_pd for a column if there is a below zero value in that column
      if ( any( field_np1(i,:) < zero ) ) then
        flux_pd(i,:) = ( flux_lim(i,:) - flux_np1(i,:) ) / dt
      else
        flux_pd(i,:) = zero
      end if
    end do

    field_nonlim = field_np1

    ! Apply change to field at n+1
    if ( field_grid == "zt" ) then

      field_np1 = -dt * ddzm( nz, ngrdcol, gr, flux_lim - flux_np1 ) + field_np1

    else if ( field_grid == "zm" ) then

      field_np1 = -dt * ddzt( nz, ngrdcol, gr, flux_lim - flux_np1 ) + field_np1

    end if

    ! Determine the total time tendency in field due to this calculation
    ! (for diagnostic purposes)
    do i = 1, ngrdcol
      ! Only set flux_pd for a column if there is a below zero value in that column
      if ( any( field_np1(i,:) < zero ) ) then
        field_pd(i,:) = ( field_np1(i,:) - field_nonlim(i,:) ) / dt
      else
        field_pd(i,:) = zero
      end if
    end do

    ! Replace the non-limited flux with the limited flux
    flux_np1 = flux_lim

    return
  end subroutine pos_definite_adj

end module pos_definite_module
