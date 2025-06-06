!-----------------------------------------------------------------------
! $Id$
!===============================================================================
module stats_type

  ! Description:
  !   Contains derived data type 'stats'.
  !   Used for storing output statistics to disk.
  !-----------------------------------------------------------------------

  use stat_file_module, only: & 
      stat_file ! Type

  use clubb_precision, only: & 
      stat_rknd,  & ! Variable(s)
      stat_nknd,  &
      core_rknd

  implicit none

  private ! Set Default Scope

  public :: stats

  ! Derived data types to store GrADS/netCDF statistics
  type stats

    ! Number of fields to sample
    integer ::  num_output_fields    ! Number of variables being output to disk (e.g.
                     ! cloud_frac, rain rate, etc.)

    integer :: &
      ii, & ! Horizontal extent of the variables (Usually 1 for the single-column model)
      jj, & ! Horizontal extent of the variables (Usually 1 for the single-column model)
      kk    ! Vertical extent of the variables (Usually gr%nzm from grid_class)

    ! Vertical levels
    real( kind = core_rknd ), allocatable, dimension(:) :: z ! altitude [m]

    ! Array to store sampled fields

    real(kind=stat_rknd), allocatable, dimension(:,:,:,:) :: accum_field_values
        ! The variable accum_field_values contains the cumulative sums
        ! of accum_num_samples sample values of each
        ! of the num_output_fields (e.g. the sum of the sampled rain rate values)

    integer(kind=stat_nknd), allocatable, dimension(:,:,:,:) :: accum_num_samples
        ! accum_num_samples is the number of samples for each of the num_output_fields fields
        ! and each of the kk vertical levels

    ! Tracks if a field is in the process of an update
    logical, allocatable, dimension(:,:,:,:) :: l_in_update

    ! Data for GrADS / netCDF output

    type (stat_file) ::  file

  end type stats

end module stats_type


