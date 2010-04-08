! $Id$
module output_2D_samples_module

  use stat_file_module, only : stat_file ! Type

  implicit none

  public :: open_2D_samples_file, close_2D_samples_file, &
    output_2D_uniform_dist_file, output_2D_lognormal_dist_file

  private ! Default scope

  type(stat_file), public :: &
    lognormal_sample_file, &
    uniform_sample_file

  contains
!-------------------------------------------------------------------------------
  subroutine open_2D_samples_file( nnzp, n_micro_calls, n_2D_variables, &
                                   fname_prefix, fdir, &
                                   time, dtwrite, zgrid, variable_names, &
                                   variable_descriptions, variable_units, &
                                   sample_file )
! Description:
!   Open a 2D sample file
! References:
!   None
!-------------------------------------------------------------------------------
#ifdef NETCDF
    use output_netcdf, only: open_netcdf ! Procedure(s)
#endif

    use stats_precision, only: time_precision ! Constant(s)

    implicit none

    ! Parameter Constants
    integer, parameter :: &
      day   = 1, & ! Made up times for GrADS
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

    ! Input/Output Variables
    type(stat_file), intent(inout) :: &
      sample_file ! File that is being opened

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
    allocate( sample_file%var(n_2D_variables) )
    allocate( sample_file%z(nnzp) )

    forall( i=1:n_micro_calls )
      rlat(i) = real( i ) ! Use made up arbitrary values for degrees north
    end forall

    rlon = 1.0 ! Also made up

    forall( i=1:n_2D_variables )
      sample_file%var(i)%name = trim( variable_names(i) )
      sample_file%var(i)%description = trim( variable_descriptions(i) )
      sample_file%var(i)%units = trim( variable_units(i) )
    end forall

#ifdef NETCDF
    call open_netcdf( nlat, nlon, fdir, fname, 1, nnzp, zgrid, &
                      day, month, year, rlat, rlon, &
                      time, dtwrite, n_2D_variables, sample_file )
#else
    stop "This version of CLUBB was not compiled for netCDF output"
#endif

    return
  end subroutine open_2D_samples_file

!-------------------------------------------------------------------------------
  subroutine output_2D_lognormal_dist_file &
             ( nnzp, n_micro_calls, d_variables, X_nl_all_levs, &
               LH_rt, LH_thl )
! Description:
!   Output a 2D snapshot of latin hypercube samples
! References:
!   None
!-------------------------------------------------------------------------------
#ifdef NETCDF
    use output_netcdf, only: write_netcdf ! Procedure(s)
#endif

    use stats_precision, only: stat_rknd ! Constant(s)

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nnzp,          & ! Number of vertical levels
      n_micro_calls, & ! Number of calls to the microphysics
      d_variables      ! Number variates being sampled

    real(kind=stat_rknd), intent(in), dimension(nnzp,n_micro_calls,d_variables) :: &
      X_nl_all_levs ! Sample that is transformed ultimately to normal-lognormal

    real, intent(in), dimension(nnzp,n_micro_calls) :: &
      LH_rt, & ! Sample of total water mixing ratio             [kg/kg]
      LH_thl   ! Sample of liquid potential temperature         [K]

    integer :: sample, j

    ! ---- Begin Code ----

    do j = 1, d_variables+2
      allocate( lognormal_sample_file%var(j)%ptr(n_micro_calls,1,nnzp) )
    end do

    do sample = 1, n_micro_calls
      do j = 1, d_variables
        lognormal_sample_file%var(j)%ptr(sample,1,1:nnzp) = X_nl_all_levs(1:nnzp,sample,j)
      end do
    end do

    ! Append rt, thl at the end of the variables
    j = d_variables+1
    do sample = 1, n_micro_calls
      lognormal_sample_file%var(j)%ptr(sample,1,1:nnzp) = LH_rt(1:nnzp,sample)
    end do

    j = d_variables+2
    do sample = 1, n_micro_calls
      lognormal_sample_file%var(j)%ptr(sample,1,1:nnzp) = LH_thl(1:nnzp,sample)
    end do

#ifdef NETCDF
    call write_netcdf( lognormal_sample_file )
#else
    stop "This version of CLUBB was not compiled for netCDF output"
#endif

    do j = 1, d_variables+2
      deallocate( lognormal_sample_file%var(j)%ptr )
    end do

    return
  end subroutine output_2D_lognormal_dist_file

!-------------------------------------------------------------------------------
  subroutine output_2D_uniform_dist_file &
             ( nnzp, n_micro_calls, dp1, X_u_all_levs, X_mixt_comp_all_levs, &
               p_matrix_s_element )
! Description:
!   Output a 2D snapshot of latin hypercube uniform distribution, i.e. (0,1)
! References:
!   None
!-------------------------------------------------------------------------------
#ifdef NETCDF
    use output_netcdf, only: write_netcdf ! Procedure(s)
#endif

    use mt95, only: genrand_real ! Constant(s)

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nnzp,          & ! Number of vertical levels
      n_micro_calls, & ! Number of calls to the microphysics
      dp1              ! Number of variates being sampled + 1

    real(kind=genrand_real), intent(in), dimension(nnzp,n_micro_calls,dp1) :: &
      X_u_all_levs ! Uniformly distributed numbers between (0,1)

    integer, intent(in), dimension(nnzp,n_micro_calls) :: &
      X_mixt_comp_all_levs ! Either 1 or 2

    integer, intent(in), dimension(n_micro_calls) :: &
      p_matrix_s_element ! P matrix at the s_mellor column

    integer :: sample, j, k

    ! ---- Begin Code ----

    do j = 1, dp1+2
      allocate( uniform_sample_file%var(j)%ptr(n_micro_calls,1,nnzp) )
    end do

    do sample = 1, n_micro_calls
      do j = 1, dp1
        uniform_sample_file%var(j)%ptr(sample,1,1:nnzp) = X_u_all_levs(1:nnzp,sample,j)
      end do
      uniform_sample_file%var(dp1+1)%ptr(sample,1,1:nnzp) = &
        real( X_mixt_comp_all_levs(1:nnzp,sample) )
      do k = 1, nnzp 
        uniform_sample_file%var(dp1+2)%ptr(sample,1,k) = real( p_matrix_s_element(sample) )
      end do
    end do

#ifdef NETCDF
    call write_netcdf( uniform_sample_file )
#else
    stop "This version of CLUBB was not compiled for netCDF output"
#endif

    do j = 1, dp1+2
      deallocate( uniform_sample_file%var(j)%ptr )
    end do

    return
  end subroutine output_2D_uniform_dist_file

!-------------------------------------------------------------------------------
  subroutine close_2D_samples_file( sample_file )
! Description:
!   Close a 2D sample file
! References:
!   None
!-------------------------------------------------------------------------------
#ifdef NETCDF
    use output_netcdf, only: close_netcdf ! Procedure
#endif

    implicit none

    type(stat_file), intent(inout) :: &
      sample_file ! File that is being opened

    ! ---- Begin Code ----

#ifdef NETCDF
    call close_netcdf( sample_file )
#else
    stop "This version of CLUBB was not compiled for netCDF output"
#endif

    deallocate( sample_file%rlat, sample_file%rlon )
    deallocate( sample_file%var)
    deallocate( sample_file%z)

    return
  end subroutine close_2D_samples_file

end module output_2D_samples_module
