!-------------------------------------------------------------------------------
! $Id$
!===============================================================================
module output_2D_samples_module

  use stat_file_module, only : stat_file ! Type

  implicit none

  public :: open_2D_samples_file, close_2D_samples_file, &
    output_2D_uniform_dist_file, output_2D_lognormal_dist_file

  private ! Default scope

  type(stat_file), public :: &
    lognormal_sample_file, &
    uniform_sample_file

  !$omp threadprivate( lognormal_sample_file, uniform_sample_file )

  contains
!-------------------------------------------------------------------------------
  subroutine open_2D_samples_file( nzt, num_samples, n_2D_variables, &
                                   fname_prefix, fdir, &
                                   time, dtwrite, zgrid, variable_names, &
                                   variable_descriptions, variable_units, &
                                   nlon, nlat, lon_vals, lat_vals, &
                                   clubb_params, sclr_dim, sclr_tol, &
                                   l_uv_nudge, &
                                   l_tke_aniso, &
                                   l_standard_term_ta, &
                                   sample_file )
! Description:
!   Open a 2D sample file
! References:
!   None
!-------------------------------------------------------------------------------
#ifdef NETCDF
    use output_netcdf, only: &
        open_netcdf_for_writing, & ! Procedure
        first_write
#endif

    use clubb_precision, only: &
        time_precision, &
        core_rknd ! Constant(s)

    use parameter_indices, only: &
        nparams    ! Variable(s)

    implicit none

    ! Parameter Constants
    integer, parameter :: &
      day   = 1, & ! Made up times for GrADS
      month = 1, &
      year  = 1900

    ! Input Variables
    integer, intent(in) :: &
      nzt,          & ! Number of vertical levels
      num_samples, & ! Number of samples per variable
      n_2D_variables ! Number variables to output

    character(len=*), intent(in) :: &
      fdir,      & ! Output directory
      fname_prefix ! Prefix for the netCDF output

    character(len=*), intent(in), dimension(n_2D_variables) :: &
      variable_names,        & ! Names of the variables to be used in the 2D netCDF file
      variable_descriptions, & ! Description of the variables in the 2D file
      variable_units           ! Units on the variables

    real(kind=time_precision), intent(in) :: & 
      time      ! Start time                      [s] 
    
    real(kind=core_rknd), intent(in) :: &
      dtwrite   ! Interval for writing to disk    [s] 

    real( kind = core_rknd ), intent(in), dimension(nzt) :: & 
      zgrid ! Vertical grid levels [m]

    ! Input/Output Variables
    type(stat_file), intent(inout) :: &
      sample_file ! File that is being opened

    integer, intent(in) :: &
      nlon, & ! Number of points in the X direction [-]
      nlat    ! Number of points in the Y direction [-]

    real( kind = core_rknd ), dimension(nlon), intent(in) ::  &
      lon_vals  ! Longitude values [Degrees E]

    real( kind = core_rknd ), dimension(nlat), intent(in) ::  &
      lat_vals  ! Latitude values  [Degrees N]

    real( kind = core_rknd ), dimension(nparams), intent(in) :: &
      clubb_params    ! Array of CLUBB's tunable parameters    [units vary]

    integer, intent(in) :: &
      sclr_dim 

    real( kind = core_rknd ), dimension(sclr_dim), intent(in) :: &
      sclr_tol

    logical, intent(in) :: &
      l_uv_nudge,         & ! For wind speed nudging
      l_tke_aniso,        & ! For anisotropic turbulent kinetic energy, i.e. TKE = 1/2
                            ! (u'^2 + v'^2 + w'^2)
      l_standard_term_ta    ! Use the standard discretization for the turbulent advection terms.
                            ! Setting to .false. means that a_1 and a_3 are pulled outside of the
                            ! derivative in advance_wp2_wp3_module.F90 and in
                            ! advance_xp2_xpyp_module.F90.

    character(len=100) :: fname
    integer :: i

    ! ---- Begin Code ----

    fname = trim( fname_prefix )//"_lh_sample_points_2D"
    i = 1  ! This assignment prevents a g 95 compiler warning

    allocate( sample_file%samples_of_var(n_2D_variables) )
    allocate( sample_file%z(nzt) )

    do i=1, n_2D_variables
      sample_file%samples_of_var(i)%name = trim( variable_names(i) )
      sample_file%samples_of_var(i)%description = trim( variable_descriptions(i) )
      sample_file%samples_of_var(i)%units = trim( variable_units(i) )
    end do

#ifdef NETCDF
    call open_netcdf_for_writing( nlat, nlon, fdir, fname, 1, nzt, zgrid, &
                      day, month, year, lat_vals, lon_vals, &
                      time, dtwrite, n_2D_variables, sample_file, num_samples )

    ! Finalize the variable definitions
    call first_write( clubb_params, sclr_dim, sclr_tol, & ! intent(in)
                      l_uv_nudge, & ! intent(in)
                      l_tke_aniso, & ! intent(in)
                      l_standard_term_ta, & ! intent(in)
                      sample_file ) ! intent(inout)
#else
    error stop "This version of CLUBB was not compiled for netCDF output"
#endif

    return
  end subroutine open_2D_samples_file

!-------------------------------------------------------------------------------
  subroutine output_2D_lognormal_dist_file ( nzt, num_samples, pdf_dim, &
                                             X_nl_all_levs, &
                                             stats_metadata )
! Description:
!   Output a 2D snapshot of latin hypercube samples
! References:
!   None
!-------------------------------------------------------------------------------
#ifdef NETCDF
    use output_netcdf, only: &
      write_netcdf ! Procedure(s)
#endif

    use clubb_precision, only: &
      stat_rknd ! Constant

    use stats_variables, only: &
      stats_metadata_type

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nzt,          & ! Number of vertical levels
      num_samples, & ! Number of samples per variable
      pdf_dim   ! Number variates being sampled

    real(kind=stat_rknd), intent(in), dimension(num_samples,nzt,pdf_dim) :: &
      X_nl_all_levs ! Sample that is transformed ultimately to normal-lognormal
      
    type (stats_metadata_type), intent(in) :: &
      stats_metadata

    integer :: sample, j

    ! ---- Begin Code ----

    if ( .not. stats_metadata%l_stats_last ) return

    do j = 1, pdf_dim
      allocate( lognormal_sample_file%samples_of_var(j)%ptr(num_samples,1,1,nzt) )
    end do

    do sample = 1, num_samples
      do j = 1, pdf_dim
        lognormal_sample_file%samples_of_var(j)%ptr(sample,1,1,1:nzt) = &
                                                                      X_nl_all_levs(sample,1:nzt,j)
      end do
    end do

#ifdef NETCDF
    call write_netcdf( lognormal_sample_file )
#else
    error stop "This version of CLUBB was not compiled for netCDF output"
#endif

    do j = 1, pdf_dim
      deallocate( lognormal_sample_file%samples_of_var(j)%ptr )
    end do

    return
  end subroutine output_2D_lognormal_dist_file

!-------------------------------------------------------------------------------
  subroutine output_2D_uniform_dist_file( nzt, num_samples, dp2, &
                                          X_u_all_levs, &
                                          X_mixt_comp_all_levs, &
                                          lh_sample_point_weights, &
                                          stats_metadata )
! Description:
!   Output a 2D snapshot of latin hypercube uniform distribution, i.e. (0,1)
! References:
!   None
!-------------------------------------------------------------------------------
#ifdef NETCDF
    use output_netcdf, only: write_netcdf ! Procedure(s)
#endif

    use clubb_precision, only: &
      core_rknd, &          ! Precision(s)
      stat_rknd

    use stats_variables, only: &
      stats_metadata_type

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nzt,          & ! Number of vertical levels
      num_samples, & ! Number of samples per variable
      dp2            ! Number of variates being sampled + 2

    real(kind=core_rknd), intent(in), dimension(num_samples,nzt,dp2) :: &
      X_u_all_levs ! Uniformly distributed numbers between (0,1)

    integer, intent(in), dimension(num_samples,nzt) :: &
      X_mixt_comp_all_levs ! Either 1 or 2

    real( kind = core_rknd ), dimension(num_samples,nzt), intent(in) :: &
      lh_sample_point_weights ! Weight of each sample

    type (stats_metadata_type), intent(in) :: &
      stats_metadata

    integer :: sample, j, k

    ! ---- Begin Code ----

    if ( .not. stats_metadata%l_stats_last ) return

    do j = 1, dp2+2
      allocate( uniform_sample_file%samples_of_var(j)%ptr(num_samples,1,1,nzt) )
    end do

    do sample = 1, num_samples
      do j = 1, dp2
        uniform_sample_file%samples_of_var(j)%ptr(sample,1,1,1:nzt) = &
          real( X_u_all_levs(sample,1:nzt,j), kind = stat_rknd )
      end do
      uniform_sample_file%samples_of_var(dp2+1)%ptr(sample,1,1,1:nzt) = &
        real( X_mixt_comp_all_levs(sample,1:nzt), kind=stat_rknd )
      do k = 1, nzt 
        uniform_sample_file%samples_of_var(dp2+2)%ptr(sample,1,1,k) = &
          real( lh_sample_point_weights(sample,k), kind=stat_rknd )
      end do
    end do

#ifdef NETCDF
    call write_netcdf( uniform_sample_file )
#else
    error stop "This version of CLUBB was not compiled for netCDF output"
#endif

    do j = 1, dp2+2
      deallocate( uniform_sample_file%samples_of_var(j)%ptr )
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
    error stop "This version of CLUBB was not compiled for netCDF output"
#endif

    deallocate( sample_file%lat_vals, sample_file%lon_vals, sample_file%samp_idx )
    deallocate( sample_file%samples_of_var)
    deallocate( sample_file%z)

    return
  end subroutine close_2D_samples_file

end module output_2D_samples_module
