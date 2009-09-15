! $Id$
module output_2D_samples_mod

  use stat_file_module, only : stat_file ! Type

  implicit none

  public :: open_2D_samples_file, output_2D_samples_file, close_2D_samples_file

  private ! Default scope

  type(stat_file), private :: sample_file

  contains
!-------------------------------------------------------------------------------
  subroutine open_2D_samples_file( nnzp, n_micro_calls, n_2D_variables, &
                                   fname_prefix, fdir, &
                                   time, dtwrite, zgrid, variable_names, &
                                   variable_descriptions, variable_units )
! Description:
!   Open a 2D sample file
! References:
!   None
!-------------------------------------------------------------------------------
    use output_netcdf, only: open_netcdf ! Procedure(s)

    use stats_precision, only: time_precision ! Constant(s)

    implicit none

    ! Parameter Constants
    integer, parameter :: &
      day   = 1, &
      month = 1, &
      year  = 1900

    ! Input Variables
    integer, intent(in) :: &
      nnzp,          & ! Number of vertical levels
      n_micro_calls, & ! Number of calls to the microphysics
      n_2D_variables   ! Number variables to output

    character(len=*), intent(in) :: &
      fdir,      & ! Output directory
      fname_prefix ! Prefix for the netCDF output

    character(len=*), intent(in), dimension(n_2D_variables) :: &
      variable_names,        & ! Names of the variables to be used in the 2D netCDF file
      variable_descriptions, & ! Description of the variables in the 2D file
      variable_units           ! Units on the variables  

    real(kind=time_precision), intent(in) :: &
      time,   & ! Start time                      [s]
      dtwrite   ! Interval for writing to disk    [s]

    real, intent(in), dimension(nnzp) :: &
      zgrid ! Vertical grid levels [m]

    ! Local Variables
    integer :: nlat, nlon ! Not actually latitudes and longitudes

    real, dimension(n_micro_calls) :: rlat

    real, dimension(1) :: rlon

    character(len=100) :: fname
    integer :: i

    ! ---- Begin Code ----
    fname = trim( fname_prefix )//"_LH_sample_points_2D"

    ! We need to set this like a latitude to trick GrADS and allow of viewing of
    ! the sample points with the GrADS application and sdfopen.
    nlat = n_micro_calls
    nlon = 1

    allocate( sample_file%rlat(n_micro_calls), sample_file%rlon(1) )
    allocate( sample_file%var( n_2D_variables ) )
    allocate( sample_file%z( nnzp ) )

    forall( i=1:n_micro_calls )
      rlat(i) = real( i ) ! Use made up arbitrary values for degrees north
    end forall

    rlon = 1.0 ! Also made up

    forall( i=1:n_2D_variables )
      sample_file%var(i)%name = trim( variable_names(i) )
      sample_file%var(i)%description = trim( variable_descriptions(i) )
      sample_file%var(i)%units = trim( variable_units(i) )
    end forall

    call open_netcdf( nlat, nlon, fdir, fname, 1, nnzp, zgrid, &
                      day, month, year, rlat, rlon, &
                      time, dtwrite, n_2D_variables, sample_file )

    return
  end subroutine open_2D_samples_file

!-------------------------------------------------------------------------------
  subroutine output_2D_samples_file( nnzp, n_micro_calls, d_variables, X_nl_all_levs, &
                                     LH_rt, LH_thl )
! Description:
!   Output a 2D snapshot of latin hypercube samples
! References:
!   None
!-------------------------------------------------------------------------------

    use output_netcdf, only: write_netcdf ! Procedure(s)

    use stats_precision, only: stat_rknd ! Constant(s)

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nnzp,          & ! Number of vertical levels
      n_micro_calls, & ! Number of calls to the microphysics
      d_variables      ! Number variates being sampled

    real(kind=stat_rknd), intent(in), dimension(nnzp,n_micro_calls,d_variables) :: &
      X_nl_all_levs ! Sample that is transformed ultimately to normal-lognormal

    real(kind=stat_rknd), intent(in), dimension(nnzp,n_micro_calls) :: &
      LH_rt, & ! Sample of total water mixing ratio             [kg/kg]
      LH_thl   ! Sample of liquid potential temperature         [K]

    integer :: i, j

    do j = 1, d_variables+2
      allocate( sample_file%var(j)%ptr(n_micro_calls,1,nnzp) )
    end do

    do i = 1, n_micro_calls
      do j = 1, d_variables
        sample_file%var(j)%ptr(i,1,1:nnzp) = X_nl_all_levs(1:nnzp,i,j)
      end do
    end do

    ! Append rt, thl at the end of the variables
    j = d_variables+1
    do i = 1, n_micro_calls
      sample_file%var(j)%ptr(i,1,1:nnzp) = LH_rt(1:nnzp,i)
    end do

    j = d_variables+2
    do i = 1, n_micro_calls
      sample_file%var(j)%ptr(i,1,1:nnzp) = LH_thl(1:nnzp,i)
    end do

    call write_netcdf( sample_file )

    do j = 1, d_variables+2
      deallocate( sample_file%var(j)%ptr )
    end do

    return
  end subroutine output_2D_samples_file

!-------------------------------------------------------------------------------
  subroutine close_2D_samples_file( )
! Description:
!   Close a 2D sample file
! References:
!   None
!-------------------------------------------------------------------------------
    use output_netcdf, only: close_netcdf ! Procedure

    implicit none

    call close_netcdf( sample_file )

    deallocate( sample_file%rlat, sample_file%rlon )
    deallocate( sample_file%var)
    deallocate( sample_file%z)

    return
  end subroutine close_2D_samples_file

end module output_2D_samples_mod
