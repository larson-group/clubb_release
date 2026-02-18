!-------------------------------------------------------------------------------
! stats_netcdf.F90
!
! Single-file NetCDF stats output with name-based updates.
!
! Namelist format (entries stored as "name | grid | units | long_name"):
!   &clubb_stats_nl
!     entry(1) = "thlm | zt | K | Liquid water potential temperature (theta_l)"
!     entry(2) = "thvm | zt | K | Virtual potential temperature"
!     entry(3) = "rtm | zt | kg/kg | Total (vapor+liquid) water mixing ratio"
!     ...
!   /
!
!
! NetCDF layout:
!   dimensions (col, altitude, time): 
!     time: tracks time sample was taken
!     col: the column of the sample (ngrdcol)
!     altitude can be one of:
!       zt: thermodynamic levels
!       zm: momentum levels
!       rad_zt: radiation thermodynamic levels
!       rad_zm: radiation momentum levels
!       sfc: surface level
!   variables:
!       surface variables: (col, time)
!       field variables: (col, altitude, time)
!       optional clubb_params(col, nparams)
!
!
! Usage:
!
!   Call once:
!     - stats_init_api(): initialize stats context and optional netcdf output file
!
!     - stats_finalize_api(): deallocate stats resources and close netcdf file
!
!
!   Call each timestep
!     - stats_begin_timestep_api(): begin timestep and set sample/write flags
! 
!     - stats_end_timestep_api(): finalize timestep and write averaged samples
!               
!
!   Call for each variable each timestep
!     - stats_update(): accumulate sample values for a registered variable
!
!     - stats_begin_budget(): start budget accumulation for a registered variable
!
!     - stats_update_budget(): add incremental contribution to active budget term
!
!     - stats_finalize_budget(): close budget term and mark as complete
!
!
!   Call to query if a variable is on the output list 
!     - var_on_stats_list(): query whether a variable is registered for output
!
!
!   Call with updated grid data, only if using adaptive gridding
!     - stats_update_grid(): update adaptive-grid remapping inputs for this step
!
!
!   Call to output silhs subcolumn samples
!     - stats_lh_samples_init(): define optional SILHS sample output variables
!
!     - stats_lh_samples_write_lognormal(): write SILHS lognormal sample fields
!
!     - stats_lh_samples_write_uniform(): write SILHS uniform sample fields
!
!-------------------------------------------------------------------------------
module stats_netcdf

  use clubb_precision, only: &
      core_rknd, &
      stat_rknd, &
      time_precision
      
  use constants_clubb, only: &
      fstderr

  use netcdf, only: &
      NF90_NOERR, NF90_CLOBBER, NF90_UNLIMITED, NF90_DOUBLE, NF90_CHAR, &
      nf90_create, nf90_redef, nf90_def_dim, nf90_def_var, nf90_def_var_fill, &
      nf90_put_att, nf90_enddef, &
      nf90_put_var, nf90_close, nf90_strerror

  use remapping_module, only: &
      remap_vals_to_target

  use err_info_type_module, only: &
      err_info_type

  use error_code, only: &
      clubb_fatal_error

  use grid_class, only: &
      grid

  implicit none

  private

  public :: stats_type, &
            stats_init_api, stats_finalize_api, &
            stats_begin_timestep_api, stats_end_timestep_api, &
            stats_update, &
            stats_begin_budget, stats_update_budget, stats_finalize_budget, &
            stats_update_grid, &
            stats_lh_samples_init, stats_lh_samples_write_lognormal, &
            stats_lh_samples_write_uniform, &
            var_on_stats_list

  interface stats_update
    module procedure stats_update_scalar
    module procedure stats_update_1d
    module procedure stats_update_2d
  end interface

  interface stats_begin_budget
    module procedure stats_begin_budget_scalar
    module procedure stats_begin_budget_1d
    module procedure stats_begin_budget_2d
  end interface

  interface stats_update_budget
    module procedure stats_update_budget_scalar
    module procedure stats_update_budget_1d
    module procedure stats_update_budget_2d
  end interface

  interface stats_finalize_budget
    module procedure stats_finalize_budget_scalar
    module procedure stats_finalize_budget_1d
    module procedure stats_finalize_budget_2d
  end interface

  ! Various constants
  integer, parameter :: NAME_LEN = 64
  integer, parameter :: GRID_LEN = 16
  integer, parameter :: UNITS_LEN = 32
  integer, parameter :: DESC_LEN = 128
  integer, parameter :: REG_LINE_LEN = 512
  integer, parameter :: NML_REG_MAX_ENTRIES = 12000

  ! Enumerators for the grid type
  integer, parameter :: GRID_ZT     = 1
  integer, parameter :: GRID_ZM     = 2
  integer, parameter :: GRID_RAD_ZT = 3
  integer, parameter :: GRID_RAD_ZM = 4
  integer, parameter :: GRID_SFC    = 5
  integer, parameter :: GRID_LH_ZT  = 6
  integer, parameter :: GRID_LH_SFC = 7

  ! A stats variable type, one of these types is used for each tracked stats variable
  type stats_var
    character(len=NAME_LEN) :: name = ''      ! Name of variable
    character(len=GRID_LEN) :: grid = ''      ! Name of grid the variable lives on (zt,zm,..)
    character(len=UNITS_LEN) :: units = ''    ! Variable unit string
    character(len=DESC_LEN) :: long_name = '' ! Long name/description of variable
    integer :: grid_id = GRID_ZT              ! NetCDF ID of the grid this variable is on
    integer :: varid = -1                     ! NetCDF ID of variable itself
    integer :: nz = 0                         ! Number of vertical levels in sample
    integer :: out_nz = 0                     ! Number of vertical levels in netcdf file 
                                              ! usually out_nz=nz, unless using adaptive gridding
    logical :: l_budget = .false.             ! True if variable is in the middle of a budget
    logical :: l_in_budget = .false.          ! Set to true when first used in budget mode,
                                              ! only for warnings about mixing up call types

    real(kind=stat_rknd), allocatable, dimension(:,:) :: &
      buffer                                  ! Buffers store data for columns and levels in a
                                              ! sample dimension( ncol, nz )
    integer, allocatable, dimension(:,:) :: &
      nsamples                                ! Number of samples, tracked pointwise

  end type stats_var

  ! A module used exclusively for adaptive regridding
  type stats_grid_module
    logical :: is_initialized = .false.           ! True when grid module is ready
    integer :: grid_remap_method = 0              ! Selected remapping algorithm

    type(grid) :: gr_source                       ! Running grid definition
    type(grid) :: gr_output                       ! Output grid definition

    real(kind=core_rknd), allocatable, dimension(:,:) :: &
      rho_lin_spline_vals,  &                     ! (col, idx) precomputed density values
      rho_lin_spline_levels                       ! (col, idx) density levels for spline fit

    real(kind=core_rknd), allocatable, dimension(:) :: &
      p_sfc                                       ! (col) surface pressure for remapping

  end type stats_grid_module

  ! A module to output stats
  type stats_lh_2d
    logical :: is_initialized = .false.           ! True when SILHS output is ready
    integer :: sample_dimid = -1                  ! NetCDF dim id for sample index
    integer :: lh_zt_dimid = -1                   ! NetCDF dim id for lh_zt
    integer :: lh_zt_varid = -1                   ! NetCDF var id for lh_zt
    integer :: num_samples = 0                    ! Number of SILHS samples
    integer :: nzt = 0                            ! Number of lh_zt levels
    integer :: n_nl_vars = 0                      ! Number of lognormal sample vars
    integer :: n_u_vars = 0                       ! Number of uniform sample vars

    integer, allocatable, dimension(:) :: &
      nl_varids, &                                ! NetCDF var ids for lognormal fields
      u_varids                                    ! NetCDF var ids for uniform fields
      
    character(len=NAME_LEN), allocatable, dimension(:) :: &
      nl_names, &                                 ! Lognormal field names
      u_names                                     ! Uniform field names

  end type stats_lh_2d

  ! A module that tracks variable IDs in a cache, used to fast-match a string to an ID
  ! otherwise the string lookups are expensive
  type stats_varname_lookup_cache
    integer :: head = 1                           ! Next write slot in hit cache
    integer :: cache_len = 0                      ! Number of valid hit entries
    integer :: reject_head = 1                    ! Next write slot in miss cache
    integer :: reject_cache_len = 0               ! Number of valid miss entries

    integer, allocatable, dimension(:) :: &
      cache                                       ! Ring buffer of recent hits, stores matched IDs

    character(len=NAME_LEN), allocatable, dimension(:) :: &
      reject_cache                                ! Recent misses, stores rejected strings
      
  end type stats_varname_lookup_cache

  ! Overall stats type, containing other various types.
  type stats_type
    logical :: enabled = .false.                  ! True if this type has been setup
    logical :: l_sample = .false.                 ! True if we want to sample the current timestep
    logical :: l_last_sample = .false.            ! True if we want to output at the end of this sample
    logical :: l_netcdf_output = .true.           ! True if netcdf file I/O is enabled
    logical :: l_different_output_grid = .false.  ! True if we need to remap the variables before output
    logical :: l_output_rad_files = .false.       ! True if we want to output radiation grid variables
    integer :: time_index = 0                     ! Current time index
    integer :: nvars = 0                          ! Number of variables setup to be output
    integer :: ncol = 0                           ! Number of columns tracked
    integer :: stats_nsamp = 1                    ! Sampling interval in timesteps
    integer :: stats_nout = 1                     ! Output interval in timesteps
    integer :: samples_per_write = 1              ! Number of samples to take before writing to disk

    type(stats_var), allocatable, dimension(:) :: &
      vars                                        ! Array of variable types
      
    type(stats_grid_module) :: grid               ! Module used for regridding
    type(stats_lh_2d) :: lh_2d                    ! Special structure to output silhs data
    type(stats_varname_lookup_cache) :: lookup    ! Cache to speed up string searches

    ! Various ids. might be nicer in its own type, but we've exhausted out percent-symbol budget
    integer :: ncid = -1                          ! NetCDF file id
    integer :: time_dimid = -1                    ! NetCDF dim id for time
    integer :: time_varid = -1                    ! NetCDF var id for time
    integer :: col_dimid = -1                     ! NetCDF dim id for columns
    integer :: col_varid = -1                     ! NetCDF var id for columns
    integer :: zt_dimid = -1                      ! NetCDF dim id for zt
    integer :: zm_dimid = -1                      ! NetCDF dim id for zm
    integer :: rad_zt_dimid = -1                  ! NetCDF dim id for rad_zt
    integer :: rad_zm_dimid = -1                  ! NetCDF dim id for rad_zm
    integer :: sfc_dimid = -1                     ! NetCDF dim id for sfc
    integer :: zt_varid = -1                      ! NetCDF var id for zt
    integer :: zm_varid = -1                      ! NetCDF var id for zm
    integer :: rad_zt_varid = -1                  ! NetCDF var id for rad_zt
    integer :: rad_zm_varid = -1                  ! NetCDF var id for rad_zm
    integer :: param_strlen_dimid = -1            ! NetCDF dim id for parameter name length
    integer :: param_dimid = -1                   ! NetCDF dim id for parameter index
    integer :: param_varid = -1                   ! NetCDF var id for parameter index
    integer :: param_name_varid = -1              ! NetCDF var id for parameter names
    integer :: clubb_params_varid = -1            ! NetCDF var id for parameter values
  end type stats_type

contains

  !----------------------------------------------------------------------------
  !                   Initialization/Finalize
  !----------------------------------------------------------------------------

  subroutine stats_init_api( registry_path, output_path, ncol, stats_tsamp, stats_tout, dt_main, &
                         day_in, month_in, year_in, time_initial, zt, zm, stats, err_info, &
                         clubb_params, param_names, rad_zt, rad_zm, hydromet_list, &
                         sclr_dim, edsclr_dim, output_zt, output_zm, grid_remap_method )

    ! Description:
    !   Initialize the stats context for one output file. This routine validates
    !   user-provided options, reads and expands the registry, configures optional
    !   grid-remapping state, and opens/defines the NetCDF file and variables.
    !   On success, stats is fully configured and ready for timestep sampling/writes.
    !-------------------------------------------------------------------------------

    ! ------------------------ Inputs ------------------------
    integer, intent(in) :: &
      ncol, &               ! Number of horizontal columns [-]
      day_in, &             ! Start day of month [day]
      month_in, &           ! Start month of year [month]
      year_in               ! Start calendar year [year]

    character(len=*), intent(in) :: &
      registry_path, &      ! Path to stats registry namelist file [-]
      output_path           ! NetCDF output path (empty disables file output) [-]

    real(kind=core_rknd), intent(in) :: &
      stats_tsamp, &        ! Sampling period [s]
      stats_tout, &         ! Output write period [s]
      dt_main               ! Model timestep [s]

    real(kind=core_rknd), dimension(:), intent(in) :: &
      zt, &                 ! Thermodynamic grid heights [m]
      zm                    ! Momentum grid heights [m]

    real(kind=time_precision), intent(in) :: &
      time_initial          ! Start time relative to start date [s]

    ! ------------------------ InOuts ------------------------
    type(stats_type), intent(inout) :: &
      stats                 ! Stats runtime context to initialize/update [-]

    type(err_info_type), intent(inout) :: &
      err_info              ! Shared CLUBB error state container [-]

    ! ------------------------ Optional Ins ------------------------
    real(kind=core_rknd), dimension(:,:), intent(in), optional :: &
      clubb_params          ! Parameter values to write as metadata [model dependent]

    real(kind=core_rknd), dimension(:), intent(in), optional :: &
      rad_zt, &             ! Radiation thermodynamic heights [m]
      rad_zm, &             ! Radiation momentum heights [m]
      output_zt, &          ! Optional output thermodynamic heights for remapping [m]
      output_zm             ! Optional output momentum heights for remapping [m]

    character(len=*), dimension(:), intent(in), optional :: &
      param_names, &        ! Names matching clubb_params columns [-]
      hydromet_list         ! Hydrometeor species names for registry expansion [-]

    integer, intent(in), optional :: &
      sclr_dim, &           ! Number of scalar species for expansions [-]
      edsclr_dim, &         ! Number of extra diagnostic scalar species [-]
      grid_remap_method     ! Remapping-method selector enum [-]

    ! ------------------------ Locals ------------------------
    type(stats_var), allocatable :: &
      defs(:)

    character(len=64) :: &
      time_units

    integer :: &
      ndefs, grid_err, ierr

    ! ------------------------ Begin Code ------------------------

    if ( stats%enabled ) then
      write( fstderr,* ) "stats error: stats_init_api called on enabled stats variable"
      err_info%err_code = clubb_fatal_error
      return
    end if

    ierr = 0

    stats%l_different_output_grid = .false.

    stats%stats_nsamp = max( 1, nint( stats_tsamp / dt_main ) )
    stats%stats_nout  = max( 1, nint( stats_tout / dt_main ) )
    stats%samples_per_write = max( 1, stats%stats_nout / stats%stats_nsamp )
    
    ! We need stats_nout/stats_nsamp to be an integer
    if ( mod( stats%stats_nout, stats%stats_nsamp ) /= 0 ) then
      write( fstderr,* ) "stats init: stats_nout not a multiple of stats_nsamp"
      err_info%err_code = clubb_fatal_error
      return
    end if

    ! Either both rad_zt and rad_zm can be defined, or neither
    if ( present( rad_zt ) .neqv. present( rad_zm ) ) then
      write( fstderr,* ) "stats error: rad_zt and rad_zm must be used together, or not at all"
      err_info%err_code = clubb_fatal_error
      return
    else
      stats%l_output_rad_files = present( rad_zt ) .and. present( rad_zm )
    end if

    ! Either both clubb_params and param_names can be defined, or neither
    if ( present( clubb_params ) .neqv. present( param_names ) ) then
      write( fstderr,* ) "stats error: clubb_params and param_names must be used together"
      err_info%err_code = clubb_fatal_error
      return
    end if

    ! If we're initializing with grid remapping, we need 
    ! grid_remap_method, output_zt, and output_zm all defined together
    if ( ( present( grid_remap_method ) .neqv. present( output_zt ) ) & 
         .or. ( present( grid_remap_method ) .neqv. present( output_zm ) ) ) then
      write( fstderr,* ) "stats init: inconsistent initialization for grid remapping"
      err_info%err_code = clubb_fatal_error
      return
    end if

    ! Figure out which variables are defined for output.
    call stats_read_registry_namelist( registry_path, defs, ndefs, ierr )

    if ( ierr /= 0 ) then
      write( fstderr,* ) err_info%err_header_global
      write( fstderr,* ) "stats init: error reading registry"
      err_info%err_code = clubb_fatal_error
      return
    end if

    ! Some variable definitions needs to be expanded into multiple variables,
    ! e.g defining "sclrm" in the registry expands to "sclr1m", "sclr2m", ...
    call stats_expand_registry( sclr_dim, edsclr_dim, &
                                defs, ndefs, hydromet_list )

    ! Get a nice string for the starting date
    call format_date( day_in, month_in, year_in, time_initial, time_units )

    ! Empty output path means collect stats only, without netcdf I/O.
    stats%l_netcdf_output = ( len_trim( output_path ) > 0 )

    ! If using grid remapping, we need to initialize the stats%grid info
    if ( stats%l_netcdf_output .and. &
         present( grid_remap_method ) .and. present( output_zt ) .and. present( output_zm ) ) then
      
      stats%l_different_output_grid = .true.

      call stats_grid_init( ncol, zt, zm, output_zt, output_zm, &
                            grid_remap_method, stats%grid, grid_err )

      if ( grid_err /= 0 ) then
        write( fstderr,* ) err_info%err_header_global
        write( fstderr,* ) "stats init: failed to initialize output grid remap"
        err_info%err_code = clubb_fatal_error
        return
      end if

    end if

    ! Initialize runtime data with output_zm/zt if using adaptive gridding.
    if ( stats%l_different_output_grid .and. stats%l_netcdf_output ) then
      call stats_type_initialize( ncol, defs, ndefs, output_zt, output_zm,  & ! Ins
                                  stats, ierr,                              & ! Inouts
                                  rad_zt, rad_zm )                            ! Optional Ins
    else
      call stats_type_initialize( ncol, defs, ndefs, zt, zm,  & ! Ins
                                  stats, ierr,                & ! Inouts
                                  rad_zt, rad_zm )              ! Optional Ins
    end if

    if ( ierr /= 0 ) then
      write( fstderr,* ) err_info%err_header_global
      write( fstderr,* ) "stats init: failed to initialize runtime state"
      err_info%err_code = clubb_fatal_error
      return
    end if

    if ( stats%l_netcdf_output ) then

      if ( stats%l_different_output_grid ) then

        ! Outputting with a different grid, using "output_zt" and "output_zm" for the netcdf
        ! file size definitions
        call stats_open_netcdf( output_path, time_units, &
                                output_zt, output_zm, stats, ierr, &
                                clubb_params=clubb_params, param_names=param_names, &
                                rad_zt=rad_zt, rad_zm=rad_zm )
      else
        
        ! Outputting with normal grid
        call stats_open_netcdf( output_path, time_units, &
                                zt, zm, stats, ierr, &
                                clubb_params=clubb_params, param_names=param_names, &
                                rad_zt=rad_zt, rad_zm=rad_zm )
      end if

      if ( ierr /= 0 ) then
        write( fstderr,* ) err_info%err_header_global
        write( fstderr,* ) "stats init: failed to open netcdf output"
        err_info%err_code = clubb_fatal_error
        return
      end if

    end if

    stats%enabled = .true.

  end subroutine stats_init_api
  
  subroutine stats_finalize_api( stats, err_info )

    ! Description:
    !   Finalize the stats context by releasing internal allocations and closing
    !   the active NetCDF file handle. This routine is safe to call when stats
    !   are disabled and leaves stats in a non-active state.
    !-------------------------------------------------------------------------------

    ! ------------------------ InOuts ------------------------
    type(stats_type), intent(inout) :: &
      stats       ! Stats runtime context to finalize [-]

    type(err_info_type), intent(inout) :: &
      err_info    ! Shared CLUBB error state container [-]

    ! ------------------------ Locals ------------------------
    integer :: ret_code, i

    ! ------------------------ Begin Code ------------------------
    ! Release variable storage even if stats was only partially initialized.
    if ( allocated( stats%vars ) ) then
      do i = 1, size( stats%vars )
        if ( allocated( stats%vars(i)%buffer ) ) deallocate( stats%vars(i)%buffer )
        if ( allocated( stats%vars(i)%nsamples ) ) deallocate( stats%vars(i)%nsamples )
      end do
      deallocate( stats%vars )
    end if

    if ( stats%lh_2d%is_initialized ) then
      if ( allocated( stats%lh_2d%nl_varids ) ) deallocate( stats%lh_2d%nl_varids )
      if ( allocated( stats%lh_2d%u_varids ) ) deallocate( stats%lh_2d%u_varids )
      if ( allocated( stats%lh_2d%nl_names ) ) deallocate( stats%lh_2d%nl_names )
      if ( allocated( stats%lh_2d%u_names ) ) deallocate( stats%lh_2d%u_names )
      stats%lh_2d%is_initialized = .false.
    end if

    if ( stats%grid%is_initialized ) then
      if ( allocated( stats%grid%rho_lin_spline_vals ) ) deallocate( stats%grid%rho_lin_spline_vals )
      if ( allocated( stats%grid%rho_lin_spline_levels ) ) deallocate( stats%grid%rho_lin_spline_levels )
      if ( allocated( stats%grid%p_sfc ) ) deallocate( stats%grid%p_sfc )
      if ( allocated( stats%grid%gr_source%zm ) ) deallocate( stats%grid%gr_source%zm )
      if ( allocated( stats%grid%gr_source%zt ) ) deallocate( stats%grid%gr_source%zt )
      if ( allocated( stats%grid%gr_output%zm ) ) deallocate( stats%grid%gr_output%zm )
      if ( allocated( stats%grid%gr_output%zt ) ) deallocate( stats%grid%gr_output%zt )
      stats%grid%is_initialized = .false.
    end if

    if ( allocated( stats%lookup%cache ) ) deallocate( stats%lookup%cache )
    if ( allocated( stats%lookup%reject_cache ) ) deallocate( stats%lookup%reject_cache )

    stats%enabled = .false.
    stats%l_sample = .false.
    stats%l_last_sample = .false.
    stats%l_netcdf_output = .true.
    stats%nvars = 0
    stats%ncol = 0
    stats%time_index = 0

    if ( stats%ncid >= 0 ) then
      ret_code = nf90_close( stats%ncid )
      if ( ret_code /= NF90_NOERR ) then
        write( fstderr,* ) err_info%err_header_global
        write( fstderr,* ) "stats finalize: netcdf close failed"
        err_info%err_code = clubb_fatal_error
        return
      end if
    end if

    stats%ncid = -1

  end subroutine stats_finalize_api

  subroutine stats_type_initialize( ncol, defs, ndefs, zt, zm, stats, ierr, rad_zt, rad_zm )

    ! Description:
    !   Initialize the mutable members of stats_type from namelist
    !   definitions and active grid sizes. Allocates per-variable runtime
    !   storage, resets name-lookup caches, and assigns nz/out_nz for each
    !   registered variable based on grid metadata.
    !-------------------------------------------------------------------------------

    ! ------------------------ Inputs ------------------------

    integer, intent(in) :: &
      ncol, &         ! Number of horizontal columns [-]
      ndefs           ! Number of registry variable definitions [-]

    type(stats_var), dimension(:), intent(in) :: &
      defs            ! Registry variable definitions to copy into runtime state [-]

    real(kind=core_rknd), dimension(:), intent(in) :: &
      zt, &           ! Thermodynamic grid heights used for output sizing [m]
      zm              ! Momentum grid heights used for output sizing [m]

    ! ------------------------ InOuts ------------------------
    type(stats_type), intent(inout) :: &
      stats           ! Stats runtime context to allocate/reset [-]

    ! ------------------------ Outputs ------------------------
    integer, intent(out) :: &
      ierr            ! Nonzero on initialization failure [-]

    ! ------------------------ Optional Ins ------------------------
    real(kind=core_rknd), dimension(:), intent(in), optional :: &
      rad_zt, &       ! Radiation thermodynamic heights [m]
      rad_zm          ! Radiation momentum heights [m]

    ! ------------------------ Locals ------------------------
    integer :: i, nzt, nzm, nrad_zt, nrad_zm

    ! ------------------------ Begin Code ------------------------
    ierr = 0

    ! Refresh global runtime sizes for this stats context.
    stats%nvars = ndefs
    stats%ncol = ncol

    ! Drop stale per-variable buffers before rebuilding variable runtime state.
    if ( allocated( stats%vars ) ) then
      do i = 1, size( stats%vars )
        if ( allocated( stats%vars(i)%buffer ) ) deallocate( stats%vars(i)%buffer )
        if ( allocated( stats%vars(i)%nsamples ) ) deallocate( stats%vars(i)%nsamples )
      end do
      deallocate( stats%vars )
    end if

    allocate( stats%vars(stats%nvars) )

    ! Rebuild positive/negative lookup caches used by stats_get_id fast paths.
    if ( allocated( stats%lookup%cache ) ) deallocate( stats%lookup%cache )
    allocate( stats%lookup%cache(max( 16,stats%nvars )) )
    stats%lookup%cache = 0

    if ( allocated( stats%lookup%reject_cache ) ) deallocate( stats%lookup%reject_cache )
    allocate( stats%lookup%reject_cache(max( 16,stats%nvars )) )
    stats%lookup%reject_cache = ''

    ! Collect active vertical sizes for each supported grid family.
    nzt = size( zt )
    nzm = size( zm )

    if ( present( rad_zt ) .and. present( rad_zm ) ) then
      nrad_zt = size( rad_zt )
      nrad_zm = size( rad_zm )
    else
      nrad_zt = 0
      nrad_zm = 0
    end if

    ! The number of vertical levels depends on the grid type, and whether or not 
    ! remapping is active
    do i = 1, ndefs

      stats%vars(i) = defs(i)

      select case ( stats%vars(i)%grid_id )

      case ( GRID_ZT )

        stats%vars(i)%out_nz = nzt

        if ( stats%l_different_output_grid ) then
          stats%vars(i)%nz = stats%grid%gr_source%nzt
        else
          stats%vars(i)%nz = stats%vars(i)%out_nz
        end if

      case ( GRID_ZM )

        stats%vars(i)%out_nz = nzm

        if ( stats%l_different_output_grid ) then
          stats%vars(i)%nz = stats%grid%gr_source%nzm
        else
          stats%vars(i)%nz = stats%vars(i)%out_nz
        end if

      case ( GRID_RAD_ZT )

        stats%vars(i)%out_nz  = nrad_zt
        stats%vars(i)%nz      = nrad_zt

      case ( GRID_RAD_ZM )

        stats%vars(i)%out_nz  = nrad_zm
        stats%vars(i)%nz      = nrad_zm

      case ( GRID_SFC, GRID_LH_SFC )
        
        stats%vars(i)%out_nz  = 1
        stats%vars(i)%nz      = 1

      case ( GRID_LH_ZT )

        ! lh_zt uses the source grid length even when adaptive remapping is active.
        if ( stats%l_different_output_grid ) then
          stats%vars(i)%out_nz  = stats%grid%gr_source%nzt
          stats%vars(i)%nz      = stats%grid%gr_source%nzt
        else
          stats%vars(i)%out_nz  = nzt
          stats%vars(i)%nz      = nzt
        end if
        
      case default
        ierr = -1
        write( fstderr,* ) "stats error in grid definition for: ", trim( stats%vars(i)%name )
        return
      end select

      ! Allocate accumulation arrays and reset sample state.
      allocate( stats%vars(i)%buffer(stats%ncol,stats%vars(i)%nz) )
      allocate( stats%vars(i)%nsamples(stats%ncol,stats%vars(i)%nz) )
      stats%vars(i)%buffer = 0.0_stat_rknd
      stats%vars(i)%nsamples = 0

    end do

  end subroutine stats_type_initialize
  
  subroutine stats_open_netcdf( output_path, time_units, &
                               zt, zm, stats, ierr, clubb_params, param_names, rad_zt, rad_zm )
    
    ! Description:
    !   Create/overwrite the NetCDF output file, define dimensions/variables
    !   from the already-prepared runtime state, and write coordinate and
    !   parameter metadata arrays.
    !-------------------------------------------------------------------------------

    ! ------------------------ Inputs ------------------------
    character(len=*), intent(in) :: &
      output_path, &            ! NetCDF output file path [-]
      time_units                ! UDUNITS time-units string for time variable [-]

    real(kind=core_rknd), dimension(:), intent(in) :: &
      zt, &                     ! Thermodynamic grid heights for NetCDF coordinate [m]
      zm                        ! Momentum grid heights for NetCDF coordinate [m]
    ! ------------------------ InOuts ------------------------
    type(stats_type), intent(inout) :: &
      stats                     ! Stats runtime context containing var definitions and IDs [-]

    ! ------------------------ Outputs ------------------------
    integer, intent(out) :: &
      ierr                      ! NetCDF/library error code (0 on success) [-]

    ! ------------------------ Optional Ins ------------------------
    character(len=*), dimension(:), intent(in), optional :: &
      param_names               ! Names matching columns of clubb_params [-]

    real(kind=core_rknd), dimension(:,:), intent(in), optional :: &
      clubb_params              ! Optional parameter matrix to write [model dependent]
      
    real(kind=core_rknd), dimension(:), intent(in), optional :: &
      rad_zt, &                 ! Radiation thermodynamic heights [m]
      rad_zm                    ! Radiation momentum heights [m]
    
    ! ------------------------ Locals ------------------------
    integer :: i, ret_code, nzt, nzm, nrad_zt, nrad_zm, nparams, lh_nzt, nparam_strlen

    ! ------------------------ Begin Code ------------------------

    ierr = 0
    if ( .not. allocated( stats%vars ) ) then
      ierr = -1
      write( fstderr,* ) "stats netcdf error: runtime state not prepared"
      return
    end if

    ! Define size of zt (nzt) and zm (nzm)
    nzt = size( zt )
    nzm = size( zm )

    ! Define size of rad_zt (nrad_zt) and rad_zm (nrad_zm)
    if ( present( rad_zt ) .and. present( rad_zm ) ) then
      nrad_zt = size( rad_zt )
      nrad_zm = size( rad_zm )
    else
      nrad_zt = 0
      nrad_zm = 0
    end if

    ! Define size of clubb_params (nparams)
    if ( present( clubb_params ) .and. present( param_names ) ) then
      nparams = size( clubb_params, 2 )
    else
      nparams = 0
    end if

    ! Define size of lh_zt (lh_nzt)
    if ( stats%l_different_output_grid ) then
      ! lh_zt variables are not setup to work be output with adaptive gridding
      ! so we use gr_source%nzt
      lh_nzt = stats%grid%gr_source%nzt
    else
      lh_nzt = nzt
    end if

    ! Create the netcdf file, NF90_CLOBBER means anything existing is overwritten
    ret_code = nf90_create( trim( output_path ), NF90_CLOBBER, stats%ncid )
    if ( ret_code /= NF90_NOERR ) then
      write( fstderr,* ) "stats netcdf error at nf90_create: ", trim( nf90_strerror( ret_code ) )
      ierr = ret_code
      return
    end if

    ! Define time dimension
    ret_code = nf90_def_dim( stats%ncid, "time", NF90_UNLIMITED, stats%time_dimid )
    ret_code = nf90_def_var( stats%ncid, "time", NF90_DOUBLE, (/ stats%time_dimid /), stats%time_varid )
    ret_code = nf90_put_att( stats%ncid, stats%time_varid, "units", trim( time_units ) )

    ! Define columns dimensions
    ret_code = nf90_def_dim( stats%ncid, "col", stats%ncol, stats%col_dimid )
    ret_code = nf90_def_var( stats%ncid, "col", NF90_DOUBLE, (/ stats%col_dimid /), stats%col_varid )
    ret_code = nf90_put_att( stats%ncid, stats%col_varid, "units", "index" )

    ! Define surface dimensions
    ret_code = nf90_def_dim( stats%ncid, "sfc", 1, stats%sfc_dimid )

    ! Define zm (momentum) dimension
    ret_code = nf90_def_dim( stats%ncid, "zm", nzm, stats%zm_dimid )
    ret_code = nf90_def_var( stats%ncid, "zm", NF90_DOUBLE, (/ stats%zm_dimid /), stats%zm_varid )
    ret_code = nf90_put_att( stats%ncid, stats%zm_varid, "units", "m" )

    ! Define zt (thermodynamic) dimension
    ret_code = nf90_def_dim( stats%ncid, "zt", nzt, stats%zt_dimid )
    ret_code = nf90_def_var( stats%ncid, "zt", NF90_DOUBLE, (/ stats%zt_dimid /), stats%zt_varid )
    ret_code = nf90_put_att( stats%ncid, stats%zt_varid, "units", "m" )

    ! Define lh_zt dimension (only if we actually have a lh_zt variable)
    if ( any( stats%vars(:)%grid_id == GRID_LH_ZT ) ) then

      ret_code = nf90_def_dim( stats%ncid, "lh_zt", lh_nzt, stats%lh_2d%lh_zt_dimid )
      ret_code = nf90_def_var( stats%ncid, "lh_zt", NF90_DOUBLE, (/ stats%lh_2d%lh_zt_dimid /), stats%lh_2d%lh_zt_varid )
      ret_code = nf90_put_att( stats%ncid, stats%lh_2d%lh_zt_varid, "units", "m" )
    end if

    ! Define rad_zt and rad_zm dimensions
    if ( nrad_zt > 0 .and. nrad_zm > 0 ) then
      ret_code = nf90_def_dim( stats%ncid, "rad_zt", nrad_zt, stats%rad_zt_dimid )
      ret_code = nf90_def_var( stats%ncid, "rad_zt", NF90_DOUBLE, (/ stats%rad_zt_dimid /), stats%rad_zt_varid )
      ret_code = nf90_put_att( stats%ncid, stats%rad_zt_varid, "units", "m" )

      ret_code = nf90_def_dim( stats%ncid, "rad_zm", nrad_zm, stats%rad_zm_dimid )
      ret_code = nf90_def_var( stats%ncid, "rad_zm", NF90_DOUBLE, (/ stats%rad_zm_dimid /), stats%rad_zm_varid )
      ret_code = nf90_put_att( stats%ncid, stats%rad_zm_varid, "units", "m" )
    end if

    ! Define the entries to store clubb_params
    if ( nparams > 0 ) then
      nparam_strlen = len( param_names(1) )
      ret_code = nf90_def_dim( stats%ncid, "param_strlen", nparam_strlen, stats%param_strlen_dimid )

      ret_code = nf90_def_var( stats%ncid, "param_name", NF90_CHAR, &
                          (/ stats%param_strlen_dimid, stats%param_dimid /), stats%param_name_varid )

      ret_code = nf90_def_dim( stats%ncid, "param", nparams, stats%param_dimid )
      ret_code = nf90_def_var( stats%ncid, "param", NF90_DOUBLE, (/ stats%param_dimid /), stats%param_varid )
      ret_code = nf90_put_att( stats%ncid, stats%param_varid, "units", "index" )

      ret_code = nf90_def_var( stats%ncid, "clubb_params", NF90_DOUBLE, &
                          (/ stats%col_dimid, stats%param_dimid /), stats%clubb_params_varid )

      ret_code = nf90_put_att( stats%ncid, stats%clubb_params_varid, "long_name", "clubb_params" )
    end if

    ! Setup individual variables
    do i = 1, stats%nvars

      ! Setup variables depending on grid type
      select case ( stats%vars(i)%grid_id )
      case ( GRID_ZT ) ! --- zt

        ret_code = nf90_def_var( stats%ncid, trim( stats%vars(i)%name ), NF90_DOUBLE, &
                            (/ stats%col_dimid, stats%zt_dimid, stats%time_dimid /), &
                            stats%vars(i)%varid )

      case ( GRID_ZM ) ! --- zm

        ret_code = nf90_def_var( stats%ncid, trim( stats%vars(i)%name ), NF90_DOUBLE, &
                            (/ stats%col_dimid, stats%zm_dimid, stats%time_dimid /), &
                            stats%vars(i)%varid )

      case ( GRID_RAD_ZT ) ! --- rad_zt

        ret_code = nf90_def_var( stats%ncid, trim( stats%vars(i)%name ), NF90_DOUBLE, &
                            (/ stats%col_dimid, stats%rad_zt_dimid, stats%time_dimid /), &
                            stats%vars(i)%varid )

      case ( GRID_RAD_ZM ) ! --- rad_zm

        ret_code = nf90_def_var( stats%ncid, trim( stats%vars(i)%name ), NF90_DOUBLE, &
                            (/ stats%col_dimid, stats%rad_zm_dimid, stats%time_dimid /), &
                            stats%vars(i)%varid )

      case ( GRID_SFC, GRID_LH_SFC ) ! --- sfc or lh_sfc

        ret_code = nf90_def_var( stats%ncid, trim( stats%vars(i)%name ), NF90_DOUBLE, &
                            (/ stats%col_dimid, stats%time_dimid /), stats%vars(i)%varid )
      
      case ( GRID_LH_ZT ) ! --- lh_zt

        ret_code = nf90_def_var( stats%ncid, trim( stats%vars(i)%name ), NF90_DOUBLE, &
                            (/ stats%col_dimid, stats%lh_2d%lh_zt_dimid, stats%time_dimid /), &
                            stats%vars(i)%varid )
      
      case default

        ! We shouldn't get here
        ierr = -1
        write( fstderr,* ) "stats error in grid definition for: ", trim( stats%vars(i)%name )
        return

      end select

      ! Add name and units to variable definitions
      ret_code = nf90_def_var_fill( stats%ncid, stats%vars(i)%varid, 0, 0.0_core_rknd )
      ret_code = nf90_put_att( stats%ncid, stats%vars(i)%varid, "long_name", trim( stats%vars(i)%long_name ) )
      ret_code = nf90_put_att( stats%ncid, stats%vars(i)%varid, "units", trim( stats%vars(i)%units ) )
    end do

    ! Leave define mode before writing
    ret_code = nf90_enddef( stats%ncid )

    ! Write our zm, zt, and column dimension values 
    ret_code = nf90_put_var( stats%ncid, stats%col_varid, [( real( i, kind=core_rknd ), i=1, stats%ncol )] )
    ret_code = nf90_put_var( stats%ncid, stats%zm_varid, zm )
    ret_code = nf90_put_var( stats%ncid, stats%zt_varid, zt )

    ! 2D silhs variables are not setup to work be output with adaptive gridding
    if ( stats%l_different_output_grid ) then
      ret_code = nf90_put_var(stats%ncid,stats%lh_2d%lh_zt_varid,stats%grid%gr_source%zt(1,:))
    else
      ret_code = nf90_put_var( stats%ncid, stats%lh_2d%lh_zt_varid, zt )
    end if

    ! Only write rad_zm and rad_zt if we're using them
    if ( nrad_zt > 0 .and. nrad_zm > 0 ) then
      ret_code = nf90_put_var( stats%ncid, stats%rad_zt_varid, rad_zt )
      ret_code = nf90_put_var( stats%ncid, stats%rad_zm_varid, rad_zm )
    end if

    ! Write clubb_params to file, this requires the param dimensions (just 1,2,3,...), and
    ! there parameter names, and the parameter values
    if ( nparams > 0 ) then
      ret_code = nf90_put_var( stats%ncid, stats%param_varid, [( real( i, kind=core_rknd ), i=1, nparams )] )
      ret_code = nf90_put_var(stats%ncid,stats%param_name_varid,param_names(1:nparams))
      ret_code = nf90_put_var( stats%ncid, stats%clubb_params_varid, clubb_params )
    end if

  end subroutine stats_open_netcdf

  !----------------------------------------------------------------------------
  !                       Begin/End timestep
  !----------------------------------------------------------------------------

  subroutine stats_begin_timestep_api( itime, stats )

    ! Description:
    !   Update per-timestep sampling/write flags from model time index. This
    !   routine sets stats%l_sample and stats%l_last_sample according to the
    !   configured sampling/output cadence and resets lookup-cache cursors so
    !   ordered name lookups start at the expected position each timestep.
    !-------------------------------------------------------------------------------

    ! ------------------------ Inputs ------------------------
    integer, intent(in) :: &
      itime               ! Current model timestep index [-]

    ! ---------------------- Input/Output ---------------------
    type(stats_type), intent(inout) :: &
      stats               ! Stats runtime context to update sample flags/caches [-]

    integer :: &
      i

    ! ------------------------ Begin Code ------------------------

    if ( .not. stats%enabled ) then
      stats%l_sample = .false.
      stats%l_last_sample = .false.
      return
    end if

    ! Sample only on stats_nsamp intervals.
    stats%l_sample = ( mod( itime, stats%stats_nsamp ) == 0 )

    ! Write only on stats_nout intervals.
    stats%l_last_sample = ( mod( itime, stats%stats_nout ) == 0 )

    ! Start each output window with clean accumulation buffers.
    if ( stats%stats_nout > 0 ) then
      if ( mod( itime - 1, stats%stats_nout ) == 0 ) then
        do i = 1, stats%nvars
          stats%vars(i)%buffer = 0.0_stat_rknd
          stats%vars(i)%nsamples = 0
        end do
      end if
    end if

    ! Replay cached expected lookup order from the beginning each timestep.
    stats%lookup%head = 1
    stats%lookup%reject_head = 1

  end subroutine stats_begin_timestep_api

  subroutine stats_end_timestep_api( time_value, stats, err_info )

    ! Description:
    !   End-of-timestep driver for stats output. If this is a write boundary,
    !   it performs consistency cleanup for in-progress budget flags and writes
    !   the current averaged record via stats_avg_and_write.
    !-------------------------------------------------------------------------------

    ! ------------------------ Inputs ------------------------
    real(kind=core_rknd), intent(in) :: &
      time_value            ! Model time written to the NetCDF time variable [s]

    ! ---------------------- Input/Output ---------------------
    type(stats_type), intent(inout) :: &
      stats                 ! Stats runtime context for end-of-step processing [-]

    type(err_info_type), intent(inout) :: &
      err_info              ! Shared CLUBB error state container [-]

    ! ------------------------- Locals ------------------------
    integer :: ierr, i

    ! ------------------------ Begin Code ------------------------

    ierr = 0
    if ( .not. stats%enabled ) return
    if ( .not. stats%l_last_sample ) return

    ! Lightweight consistency check: no variable should remain mid-budget
    ! after timestep logic has completed.
    do i = 1, stats%nvars
      if ( .not. stats%vars(i)%l_in_budget ) cycle
      write( fstderr,* ) "stats warning: budget left mid-update at end_timestep for ", &
                       trim( stats%vars(i)%name )
      stats%vars(i)%l_in_budget = .false.
    end do

    ! Divide the buffer storage by the number of samples, 
    ! and write to the netcdf file if stats%l_netcdf_output=.true.
    call stats_avg_and_write( time_value, stats, ierr )
    
    ! Reset sample and write flags
    stats%l_sample = .false.
    stats%l_last_sample = .false.

  end subroutine stats_end_timestep_api

  !----------------------------------------------------------------------------
  !                         Sampling/Writing subroutines
  !----------------------------------------------------------------------------

  subroutine stats_update_scalar( name, values, stats, icol, level )

    ! Description:
    !   Add sampled scalar values to the accumulation buffer for a named
    !   variable. Optional icol/level indexing determines whether this updates
    !   a specific surface column point or a specific profile point.
    !-------------------------------------------------------------------------------

    ! ------------------------ Inputs ------------------------
    character(len=*), intent(in) :: &
      name              ! Registered variable name to update [-]

    real(kind=core_rknd), intent(in) :: &
      values            ! Sample value in variable units [var-dependent]

    ! ---------------------- Input/Output ---------------------
    type(stats_type), intent(inout) :: &
      stats             ! Stats runtime context and accumulation buffers [-]

    ! ------------------------ Optional Ins ------------------------
    integer, intent(in), optional :: &
      icol, &           ! Column index for single-point update [-]
      level             ! Vertical index for profile single-point update [-]

    ! ------------------------- Locals ------------------------
    integer :: id

    ! ------------------------ Begin Code ------------------------

    if ( .not. stats%enabled ) return
    id = stats_get_id( stats, name )
    if ( id <= 0 ) return
    if ( .not. stats%l_sample ) return

    if ( stats%vars(id)%l_budget ) then
      write( fstderr,* ) "stats error: called stat_update for variable" // &
                       " currently in budget: ", trim( stats%vars(id)%name )
      return
    end if

    !$acc update host( values ) if_present

    if ( present( level ) ) then
      ! Updating a specific level and column in 2D field
      stats%vars(id)%buffer(icol,level) = stats%vars(id)%buffer(icol,level) + values
      stats%vars(id)%nsamples(icol,level) = stats%vars(id)%nsamples(icol,level) + 1
    else if ( present( icol ) ) then
      ! Updating a specific column in a surface field
      stats%vars(id)%buffer(icol,1) = stats%vars(id)%buffer(icol,1) + values
      stats%vars(id)%nsamples(icol,1) = stats%vars(id)%nsamples(icol,1) + 1
    end if

    if ( any( stats%vars(id)%nsamples > stats%samples_per_write ) ) then
      write( fstderr,* ) "stats oversampling warning for ", trim( stats%vars(id)%name ), &
                       ": max_count=", maxval( stats%vars(id)%nsamples ), &
                       " expected <=", stats%samples_per_write
    end if

  end subroutine stats_update_scalar

  subroutine stats_update_1d( name, values, stats, icol, level )

    ! Description:
    !   Add sampled 1D values to the accumulation buffer for a named variable.
    !   For profile variables with icol present this updates one column profile;
    !   without icol this updates all surface columns.
    !-------------------------------------------------------------------------------

    ! ------------------------ Inputs ------------------------
    character(len=*), intent(in) :: &
      name                ! Registered variable name to update [-]

    real(kind=core_rknd), intent(in), dimension(:) :: &
      values              ! 1D payload in variable units [var-dependent]

    ! ---------------------- Input/Output ---------------------
    type(stats_type), intent(inout) :: &
      stats               ! Stats runtime context and accumulation buffers [-]

    ! ------------------------ Optional Ins ------------------------
    integer, intent(in), optional :: &
      icol, &             ! Column index when values represents one-column profile [-]
      level               ! Level index when values represents all columns at one level [-]

    ! ------------------------- Locals ------------------------
    integer :: id

    ! ------------------------ Begin Code ------------------------

    if ( .not. stats%enabled ) return
    id = stats_get_id( stats, name )
    if ( id <= 0 ) return
    if ( .not. stats%l_sample ) return

    if ( stats%vars(id)%l_budget ) then
      write( fstderr,* ) "stats error: called stat_update for variable" // &
                       " currently in budget: ", trim( stats%vars(id)%name )
      return
    end if

    !$acc update host( values ) if_present

    if ( present( icol ) ) then
      ! Updating a specific column in a profile field
      stats%vars(id)%buffer(icol,:) = stats%vars(id)%buffer(icol,:) + values(:)
      stats%vars(id)%nsamples(icol,:) = stats%vars(id)%nsamples(icol,:) + 1
    else
      ! Updating all columns in a surface field
      stats%vars(id)%buffer(:,1) = stats%vars(id)%buffer(:,1) + values(:)
      stats%vars(id)%nsamples(:,1) = stats%vars(id)%nsamples(:,1) + 1
    end if

    if ( any( stats%vars(id)%nsamples > stats%samples_per_write ) ) then
      write( fstderr,* ) "stats oversampling warning for ", trim( stats%vars(id)%name ), &
                       ": max_count=", maxval( stats%vars(id)%nsamples ), &
                       " expected <=", stats%samples_per_write
    end if

  end subroutine stats_update_1d

  subroutine stats_update_2d( name, values, stats, icol, level )

    ! Description:
    !   Add sampled full-field 2D values to the accumulation buffer for a named
    !   profile variable. This path updates all columns and all levels at once.
    !-------------------------------------------------------------------------------

    ! ------------------------ Inputs ------------------------
    character(len=*), intent(in) :: &
      name            ! Registered variable name to update [-]

    real(kind=core_rknd), intent(in), dimension(:,:) :: &
      values          ! 2D payload (col,level) in variable units [var-dependent]

    ! ---------------------- Input/Output ---------------------
    type(stats_type), intent(inout) :: &
      stats           ! Stats runtime context and accumulation buffers [-]

    ! ------------------------ Optional Ins ------------------------
    integer, intent(in), optional :: &
      icol, &         ! Optional single-column selector (unused in 2D path) [-]
      level           ! Optional single-level selector (unused in 2D path) [-]

    ! ------------------------- Locals ------------------------
    integer :: id

    ! ------------------------ Begin Code ------------------------

    if ( .not. stats%enabled ) return
    id = stats_get_id( stats, name )
    if ( id <= 0 ) return
    if ( .not. stats%l_sample ) return

    if ( stats%vars(id)%l_budget ) then
      write( fstderr,* ) "stats error: called stat_update for variable" // &
                       " currently in budget: ", trim( stats%vars(id)%name )
      return
    end if
    
    if ( size( values, 1 ) /= size( stats%vars(id)%buffer, 1 ) &
         .or. size( values, 2 ) /= size( stats%vars(id)%buffer, 2 ) ) then
      write( fstderr,* ) "stats shape mismatch for ", trim( stats%vars(id)%name ), &
                       ": input nz =", size( values, 1 ), &
                       " expected nz =", size( stats%vars(id)%buffer, 1 ), &
                       ": input ngrdcol =", size( values, 2 ), &
                       " expected ngrdcol ", size( stats%vars(id)%buffer, 2 )
      return
    end if
    
    !$acc update host( values ) if_present

    ! Updating all columns and levels in a profile field
    stats%vars(id)%buffer(:,:) = stats%vars(id)%buffer(:,:) + values(:,:)
    stats%vars(id)%nsamples(:,:) = stats%vars(id)%nsamples(:,:) + 1

    if ( any( stats%vars(id)%nsamples > stats%samples_per_write ) ) then
      write( fstderr,* ) "stats oversampling warning for ", trim( stats%vars(id)%name ), &
                       ": max_count=", maxval( stats%vars(id)%nsamples ), &
                       " expected <=", stats%samples_per_write
    end if

  end subroutine stats_update_2d

  subroutine stats_begin_budget_scalar( name, values, stats, icol )

    ! Description:
    !   Begin budget accumulation for scalar payloads by subtracting the current
    !   state from the buffer. This is used for per-column point updates.
    !-------------------------------------------------------------------------------

    ! ------------------------ Inputs ------------------------
    character(len=*), intent(in) :: &
      name            ! Registered budget variable name [-]

    real(kind=core_rknd), intent(in) :: &
      values          ! Current state value in variable units [var-dependent]

    ! ---------------------- Input/Output ---------------------
    type(stats_type), intent(inout) :: &
      stats           ! Stats runtime context and budget buffer state [-]

    ! ------------------------ Optional Ins ------------------------
    integer, intent(in), optional :: &
      icol            ! Optional column index for single-point budget start [-]

    ! ------------------------- Locals ------------------------
    integer :: id

    ! ------------------------ Begin Code ------------------------

    if ( .not. stats%enabled ) return
    id = stats_get_id( stats, name )
    if ( id <= 0 ) return
    if ( .not. stats%l_sample ) return
    
    if ( stats%vars(id)%l_in_budget ) then
      write( fstderr,* ) "stats budget begin twice for ", trim( stats%vars(id)%name )
      return
    end if

    !$acc update host( values ) if_present

    ! Subtract the current state from the buffer to begin the budget accumulation.
    stats%vars(id)%buffer(icol,1) = stats%vars(id)%buffer(icol,1) - values
    stats%vars(id)%l_budget = .true.
    stats%vars(id)%l_in_budget = .true.

  end subroutine stats_begin_budget_scalar

  subroutine stats_begin_budget_1d( name, values, stats, icol )

    ! Description:
    !   Begin budget accumulation for 1D payloads by subtracting the current
    !   vector state from the buffer.
    !-------------------------------------------------------------------------------

    ! ------------------------ Inputs ------------------------
    character(len=*), intent(in) :: &
      name                  ! Registered budget variable name [-]

    real(kind=core_rknd), intent(in), dimension(:) :: &
      values                ! 1D state payload in variable units [var-dependent]

    ! ---------------------- Input/Output ---------------------
    type(stats_type), intent(inout) :: &
      stats                 ! Stats runtime context and budget buffer state [-]

    ! ------------------------ Optional Ins ------------------------
    integer, intent(in), optional :: &
      icol                  ! Optional column index for one-column profile budget start [-]

    ! ------------------------- Locals ------------------------
    integer ::  id

    ! ------------------------ Begin Code ------------------------

    if ( .not. stats%enabled ) return
    id = stats_get_id( stats, name )
    if ( id <= 0 ) return
    if ( .not. stats%l_sample ) return
    if ( stats%vars(id)%l_in_budget ) then
      write( fstderr,* ) "stats budget begin twice for ", trim( stats%vars(id)%name )
      return
    end if

    !$acc update host( values ) if_present

    if ( present( icol ) ) then
      ! Updating a specific column in a profile field.
      ! Subtract the current state from the buffer to begin the budget accumulation.
      stats%vars(id)%buffer(icol,:) = stats%vars(id)%buffer(icol,:) - values(:)
    else
      ! Updating all columns in a surface field
      ! Subtract the current state from the buffer to begin the budget accumulation.
      stats%vars(id)%buffer(:,1) = stats%vars(id)%buffer(:,1) - values(:)
    end if

    stats%vars(id)%l_budget = .true.
    stats%vars(id)%l_in_budget = .true.

  end subroutine stats_begin_budget_1d

  subroutine stats_begin_budget_2d( name, values, stats, icol )

    ! Description:
    !   Begin budget accumulation for full-field 2D payloads by subtracting the
    !   current state from the full buffer.
    !-------------------------------------------------------------------------------

    ! ------------------------ Inputs ------------------------
    character(len=*), intent(in) :: &
      name                  ! Registered budget variable name [-]

    real(kind=core_rknd), intent(in), dimension(:,:) :: &
      values                ! 2D state payload (col,level) in variable units [var-dependent]

    ! ---------------------- Input/Output ---------------------
    type(stats_type), intent(inout) :: &
      stats                 ! Stats runtime context and budget buffer state [-]

    ! ------------------------ Optional Ins ------------------------
    integer, intent(in), optional :: &
      icol                  ! Optional column selector (unused in 2D path) [-]

    ! ------------------------- Locals ------------------------
    integer :: id

    ! ------------------------ Begin Code ------------------------ 

    if ( .not. stats%enabled ) return
    id = stats_get_id( stats, name )
    if ( id <= 0 ) return
    if ( .not. stats%l_sample ) return
    if ( stats%vars(id)%l_in_budget ) then
      write( fstderr,* ) "stats budget begin twice for ", trim( stats%vars(id)%name )
      return
    end if

    !$acc update host( values ) if_present

    ! Subtract the current state from the buffer to begin the budget accumulation.
    stats%vars(id)%buffer(:,:) = stats%vars(id)%buffer(:,:) - values(:,:)
    stats%vars(id)%l_budget = .true.
    stats%vars(id)%l_in_budget = .true.

  end subroutine stats_begin_budget_2d

  subroutine stats_update_budget_scalar( name, values, stats, icol, level )

    ! Description:
    !   Add intermediate budget contributions from scalar payloads during an
    !   active budget window.
    !-------------------------------------------------------------------------------

    ! ------------------------ Inputs ------------------------
    character(len=*), intent(in) :: &
      name                ! Registered budget variable name [-]

    real(kind=core_rknd), intent(in) :: &
      values              ! Incremental budget contribution in variable units [var-dependent]

    ! ---------------------- Input/Output ---------------------
    type(stats_type), intent(inout) :: &
      stats               ! Stats runtime context and budget buffer state [-]

    ! ------------------------ Optional Ins ------------------------
    integer, intent(in), optional :: &
      icol, &             ! Optional column index for scalar update [-]
      level               ! Optional level index for profile-point update [-]

    ! ------------------------- Locals ------------------------
    integer :: id

    ! ------------------------ Begin Code ------------------------

    if ( .not. stats%enabled ) return
    if ( .not. stats%l_sample ) return
    id = stats_get_id( stats, name )
    if ( id <= 0 ) return

    if ( .not. stats%vars(id)%l_in_budget ) then
      write( fstderr,* ) "stats budget update without begin for ", trim( stats%vars(id)%name )
      return
    end if

    !$acc update host( values ) if_present

    if ( present( level ) ) then
      ! Updating specific level and column in 2D field
      stats%vars(id)%buffer(icol,level) = stats%vars(id)%buffer(icol,level) + values
    else if ( present( icol ) ) then
      ! Updating specific column in surface field
      stats%vars(id)%buffer(icol,1) = stats%vars(id)%buffer(icol,1) + values
    end if

  end subroutine stats_update_budget_scalar

  subroutine stats_update_budget_1d( name, values, stats, icol, level )

    ! Description:
    !   Add intermediate budget contributions from 1D payloads during an active
    !   budget window.
    !-------------------------------------------------------------------------------

    ! ------------------------ Inputs ------------------------
    character(len=*), intent(in) :: &
      name                ! Registered budget variable name [-]

    real(kind=core_rknd), intent(in), dimension(:) :: &
      values              ! 1D incremental budget contribution [var-dependent]

    ! ---------------------- Input/Output ---------------------
    type(stats_type), intent(inout) :: &
      stats               ! Stats runtime context and budget buffer state [-]

    ! ------------------------ Optional Ins ------------------------
    integer, intent(in), optional :: &
      icol, &             ! Optional column index for one-column profile update [-]
      level               ! Optional level index for all-column single-level update [-]

    ! ------------------------- Locals ------------------------
    integer :: id

    ! ------------------------ Begin Code ------------------------

    if ( .not. stats%enabled ) return
    if ( .not. stats%l_sample ) return
    id = stats_get_id( stats, name )
    if ( id <= 0 ) return

    if ( .not. stats%vars(id)%l_in_budget ) then
      write( fstderr,* ) "stats budget update without begin for ", trim( stats%vars(id)%name )
      return
    end if

    !$acc update host( values ) if_present

    if ( present( icol ) ) then
      ! Updating specific column in 2D field 
      stats%vars(id)%buffer(icol,:) = stats%vars(id)%buffer(icol,:) + values(:)
    else
      ! Updating surface variable
      stats%vars(id)%buffer(:,1) = stats%vars(id)%buffer(:,1) + values(:)
    end if

  end subroutine stats_update_budget_1d

  subroutine stats_update_budget_2d( name, values, stats, icol, level )

    ! Description:
    !   Add intermediate budget contributions from full-field 2D payloads during
    !   an active budget window.
    !-------------------------------------------------------------------------------

    ! ------------------------ Inputs ------------------------
    character(len=*), intent(in) :: &
      name                ! Registered budget variable name [-]

    real(kind=core_rknd), intent(in), dimension(:,:) :: &
      values              ! 2D incremental budget contribution [var-dependent]

    ! ---------------------- Input/Output ---------------------
    type(stats_type), intent(inout) :: &
      stats               ! Stats runtime context and budget buffer state [-]

    ! ------------------------ Optional Ins ------------------------
    integer, intent(in), optional :: &
      icol, &             ! Optional column selector (unused in 2D path) [-]
      level               ! Optional level selector (unused in 2D path) [-]

    ! ------------------------- Locals ------------------------
    integer :: &
      id

    ! ------------------------ Begin Code ------------------------

    if ( .not. stats%enabled ) return
    if ( .not. stats%l_sample ) return
    id = stats_get_id( stats, name )
    if ( id <= 0 ) return

    !$acc update host( values ) if_present

    if ( .not. stats%vars(id)%l_in_budget ) then
      write( fstderr,* ) "stats budget update without begin for ", trim( stats%vars(id)%name )
      return
    end if

    stats%vars(id)%buffer(:,:) = stats%vars(id)%buffer(:,:) + values(:,:)

  end subroutine stats_update_budget_2d

  subroutine stats_finalize_budget_scalar( name, values, stats, icol, l_count_sample )

    ! Description:
    !   Close an active budget window using scalar payloads, and optionally
    !   increment per-element sample counters for averaging.
    !-------------------------------------------------------------------------------

    ! ------------------------ Inputs ------------------------
    character(len=*), intent(in) :: &
      name                ! Registered budget variable name [-]

    real(kind=core_rknd), intent(in) :: &
      values              ! Closing state value in variable units [var-dependent]

    ! ---------------------- Input/Output ---------------------
    type(stats_type), intent(inout) :: &
      stats               ! Stats runtime context and budget/sample state [-]

    ! ------------------------ Optional Ins ------------------------
    integer, intent(in), optional :: &
      icol                ! Optional column index for single-point finalize [-]

    logical, intent(in), optional :: &
      l_count_sample      ! If true, increment sample counters for averaging [-]

    ! ------------------------- Locals ------------------------
    integer :: &
      id

    logical :: &
      do_count_sample

    ! ------------------------ Begin Code ------------------------

    if ( .not. stats%enabled ) return
    if ( .not. stats%l_sample ) return
    id = stats_get_id( stats, name )
    if ( id <= 0 ) return

    ! Don't sample if l_count_sample is false, useful for modifying without counting a sample
    do_count_sample = .true.
    if ( present( l_count_sample ) ) do_count_sample = l_count_sample

    if ( .not. stats%vars(id)%l_in_budget ) then
      write( fstderr,* ) "stats budget finalize without begin for ", trim( stats%vars(id)%name )
      return
    end if

    !$acc update host( values ) if_present

    stats%vars(id)%buffer(icol,1) = stats%vars(id)%buffer(icol,1) + values
    if ( do_count_sample ) stats%vars(id)%nsamples(icol,1) = stats%vars(id)%nsamples(icol,1) + 1

    stats%vars(id)%l_budget = .true.
    stats%vars(id)%l_in_budget = .false.

    if ( do_count_sample ) then
      if ( any( stats%vars(id)%nsamples > stats%samples_per_write ) ) then
        write( fstderr,* ) "stats oversampling warning for ", trim( stats%vars(id)%name ), &
                        ": max_count=", maxval( stats%vars(id)%nsamples ), &
                        " expected <=", stats%samples_per_write
      end if
    end if

  end subroutine stats_finalize_budget_scalar

  subroutine stats_finalize_budget_1d( name, values, stats, icol, l_count_sample )

    ! Description:
    !   Close an active budget window using 1D payloads, and optionally
    !   increment per-element sample counters for averaging.
    !-------------------------------------------------------------------------------

    ! ------------------------ Inputs ------------------------
    character(len=*), intent(in) :: &
      name                ! Registered budget variable name [-]

    real(kind=core_rknd), intent(in), dimension(:) :: &
      values              ! 1D closing state payload in variable units [var-dependent]

    ! ---------------------- Input/Output ---------------------
    type(stats_type), intent(inout) :: &
      stats               ! Stats runtime context and budget/sample state [-]

    ! ------------------------ Optional Ins ------------------------
    integer, intent(in), optional :: &
      icol                ! Optional column index for one-column profile finalize [-]

    logical, intent(in), optional :: &
      l_count_sample      ! If true, increment sample counters for averaging [-]

    ! ------------------------- Locals ------------------------
    integer :: &
      id

    logical :: &
      do_count_sample

    ! ------------------------ Begin Code ------------------------

    if ( .not. stats%enabled ) return
    if ( .not. stats%l_sample ) return
    id = stats_get_id( stats, name )
    if ( id <= 0 ) return

    ! Don't sample if l_count_sample is false, useful for modifying without counting a sample
    do_count_sample = .true.
    if ( present( l_count_sample ) ) do_count_sample = l_count_sample

    if ( .not. stats%vars(id)%l_in_budget ) then
      write( fstderr,* ) "stats budget finalize without begin for ", trim( stats%vars(id)%name )
      return
    end if

    !$acc update host( values ) if_present

    if ( present( icol ) ) then
      ! Updating a specific column in a profile field
      stats%vars(id)%buffer(icol,:) = stats%vars(id)%buffer(icol,:) + values(:)
      if ( do_count_sample ) stats%vars(id)%nsamples(icol,:) = stats%vars(id)%nsamples(icol,:) + 1
    else
      ! Updating all columns in a surface field
      stats%vars(id)%buffer(:,1) = stats%vars(id)%buffer(:,1) + values(:)
      if ( do_count_sample ) stats%vars(id)%nsamples(:,1) = stats%vars(id)%nsamples(:,1) + 1
    end if

    stats%vars(id)%l_budget = .true.
    stats%vars(id)%l_in_budget = .false.

    if ( do_count_sample ) then
      if ( any( stats%vars(id)%nsamples > stats%samples_per_write ) ) then
        write( fstderr,* ) "stats oversampling warning for ", trim( stats%vars(id)%name ), &
                        ": max_count=", maxval( stats%vars(id)%nsamples ), &
                        " expected <=", stats%samples_per_write
      end if
    end if

  end subroutine stats_finalize_budget_1d

  subroutine stats_finalize_budget_2d( name, values, stats, icol, l_count_sample )

    ! Description:
    !   Close an active budget window using full-field 2D payloads, and
    !   optionally increment per-element sample counters for averaging.
    !-------------------------------------------------------------------------------

    ! ------------------------ Inputs ------------------------
    character(len=*), intent(in) :: &
      name                ! Registered budget variable name [-]

    real(kind=core_rknd), intent(in), dimension(:,:) :: &
      values              ! 2D closing state payload in variable units [var-dependent]

    ! ---------------------- Input/Output ---------------------
    type(stats_type), intent(inout) :: &
      stats               ! Stats runtime context and budget/sample state [-]

    ! ------------------------ Optional Ins ------------------------
    integer, intent(in), optional :: &
      icol                ! Optional column selector (unused in 2D path) [-]

    logical, intent(in), optional :: &
      l_count_sample      ! If true, increment sample counters for averaging [-]

    ! ------------------------- Locals ------------------------
    integer :: &
      id

    logical :: &
      do_count_sample

    ! ------------------------ Begin Code ------------------------

    if ( .not. stats%enabled ) return
    if ( .not. stats%l_sample ) return
    id = stats_get_id( stats, name )
    if ( id <= 0 ) return

    ! Don't sample if l_count_sample is false, useful for modifying without counting a sample
    do_count_sample = .true.
    if ( present( l_count_sample ) ) do_count_sample = l_count_sample

    if ( .not. stats%vars(id)%l_in_budget ) then
      write( fstderr,* ) "stats budget finalize without begin for ", trim( stats%vars(id)%name )
      return
    end if

    !$acc update host( values ) if_present

    stats%vars(id)%buffer(:,:) = stats%vars(id)%buffer(:,:) + values(:,:)
    if ( do_count_sample ) stats%vars(id)%nsamples(:,:) = stats%vars(id)%nsamples(:,:) + 1

    stats%vars(id)%l_budget = .true.
    stats%vars(id)%l_in_budget = .false.

    if ( do_count_sample ) then
      if ( any( stats%vars(id)%nsamples > stats%samples_per_write ) ) then
        write( fstderr,* ) "stats oversampling warning for ", trim( stats%vars(id)%name ), &
                        ": max_count=", maxval( stats%vars(id)%nsamples ), &
                        " expected <=", stats%samples_per_write
      end if
    end if

  end subroutine stats_finalize_budget_2d

  subroutine stats_avg_and_write( time_value, stats, ierr )

    ! Description:
    !   Write one output record to NetCDF from buffered sampled values. The
    !   routine averages by sample counts, applies optional vertical remapping
    !   for zt/zm profile fields, and writes all registered variables.
    !-------------------------------------------------------------------------------

    ! ------------------------ In ------------------------

    real(kind=core_rknd), intent(in) :: &
      time_value          ! Model time to write for this record [s]

    ! ------------------------ InOut ------------------------
    type(stats_type), intent(inout) :: &
      stats               ! Stats runtime context and variable buffers [-]

    ! ------------------------ Out ------------------------
    integer, intent(out) :: &
      ierr                ! Nonzero on write/remap/NetCDF failure [-]

    ! ------------------------ Locals ------------------------
    integer :: ret_code, i, iv, grid_id

    real(kind=core_rknd), dimension(1) :: &
      time_buf

    ! ------------------------ Begin Code ------------------------

    ierr = 0
    if ( .not. stats%enabled ) return

    ! Increment time record and write time coordinate.
    stats%time_index = stats%time_index + 1

    ret_code = NF90_NOERR
    if ( stats%l_netcdf_output ) then
      ! Write time coordinate
      time_buf(1) = time_value
      ret_code = nf90_put_var( stats%ncid, stats%time_varid, time_buf, &
                           start=(/ stats%time_index /), count=(/ 1 /) )
    end if

    ! Loop through all variables
    do i = 1, stats%nvars

      ! Skip if variable had no samples
      if ( maxval( stats%vars(i)%nsamples ) <= 0 ) cycle

      ! Average over the sampled slots collected since the last write.
      ! Scalars are stored as (col,time); profiles as (col,altitude,time).
      grid_id = stats%vars(i)%grid_id

      ! Special case for sfc and lh_sfc variables because they have one less dimension
      if ( grid_id == GRID_SFC .or. grid_id == GRID_LH_SFC ) then

        ! Average buffer over sample count
        where ( stats%vars(i)%nsamples(:,1) > 0 )
          stats%vars(i)%buffer(:,1) = stats%vars(i)%buffer(:,1) / &
                                       real(stats%vars(i)%nsamples(:,1),kind=stat_rknd)
        end where

        if ( stats%l_netcdf_output ) then
          ! netcdf write
          ret_code = nf90_put_var( stats%ncid, stats%vars(i)%varid, stats%vars(i)%buffer(:,1), &
                              start=(/ 1, stats%time_index /), count=(/ stats%ncol, 1 /) )
        end if
      else

        ! Average buffer over sample count
        where ( stats%vars(i)%nsamples(:,:) > 0 )
          stats%vars(i)%buffer(:,:) = stats%vars(i)%buffer(:,:) / &
                                       real(stats%vars(i)%nsamples(:,:),kind=stat_rknd)
        end where

        ! Regridding requires transforming output using "remap_vals_to_target"
        if ( stats%l_different_output_grid &
             .and. ( grid_id == GRID_ZM .or. grid_id == GRID_ZT ) ) then

          iv = 1

          if ( stats%l_netcdf_output ) then
            ! netcdf write
            ret_code = nf90_put_var( &
                stats%ncid, stats%vars(i)%varid, &
                remap_vals_to_target( stats%ncol, &
                                      stats%grid%gr_source, stats%grid%gr_output, stats%vars(i)%nz, &
                                      real( stats%vars(i)%buffer, kind=core_rknd ), stats%vars(i)%out_nz, &
                                      size( stats%grid%rho_lin_spline_vals, 2 ), stats%grid%rho_lin_spline_vals, &
                                      stats%grid%rho_lin_spline_levels, iv, stats%grid%p_sfc, &
                                      stats%grid%grid_remap_method, ( grid_id == GRID_ZT ) ), &
                start=(/ 1, 1, stats%time_index /), &
                count=(/ stats%ncol, stats%vars(i)%out_nz, 1 /) )
          end if
        else

          if ( stats%l_netcdf_output ) then
            ! netcdf write
            ret_code = nf90_put_var( stats%ncid, stats%vars(i)%varid, &
                                stats%vars(i)%buffer, &
                                start=(/ 1, 1, stats%time_index /), &
                                count=(/ stats%ncol, stats%vars(i)%out_nz, 1 /) )
          end if

        end if
      end if

      if ( stats%l_netcdf_output .and. ret_code /= NF90_NOERR ) then
        write( fstderr,* ) "stats write error for ", trim( stats%vars(i)%name ), &
                         " grid=", trim( stats%vars(i)%grid ), &
                         " nz=", stats%vars(i)%nz, " out_nz=", stats%vars(i)%out_nz, &
                         " ncol=", stats%ncol, &
                         " time_index=", stats%time_index
        write( fstderr,* ) "stats write error: ", trim( nf90_strerror( ret_code ) )
        ierr = ret_code
        return
      end if

    end do

  end subroutine stats_avg_and_write

  !----------------------------------------------------------------------------
  !                         Grid Remapping Helpers
  !----------------------------------------------------------------------------

  subroutine stats_grid_init( ncol, zt_src, zm_src, zt_tgt, zm_tgt, &
                              grid_remap_method, grid_ctx, ierr )

    ! Description:
    !   Initialize source/target grid objects used by remapping. This routine
    !   validates staggered-grid consistency (nzt = nzm-1), allocates source
    !   and target altitude arrays for all columns, and stores remap-method
    !   metadata in the grid context.
    !-------------------------------------------------------------------------------

    ! ------------------------ Inputs ------------------------
    integer, intent(in) :: &
      ncol, &             ! Number of horizontal columns [-]
      grid_remap_method   ! Remapping-method selector enum [-]

    real(kind=core_rknd), dimension(:), intent(in) :: &
      zt_src, &           ! Source thermodynamic heights [m]
      zm_src, &           ! Source momentum heights [m]
      zt_tgt, &           ! Target/output thermodynamic heights [m]
      zm_tgt              ! Target/output momentum heights [m]

    ! ---------------------- Input/Output ---------------------
    type(stats_grid_module), intent(inout) :: &
      grid_ctx            ! Adaptive-grid remapping context to initialize [-]

    ! ------------------------ Outputs ------------------------
    integer, intent(out) :: &
      ierr                ! Nonzero when grid compatibility checks fail [-]

    ! ------------------------- Locals ------------------------
    integer :: &
      nzm_src, &
      nzt_src, &
      nzm_tgt, &
      nzt_tgt

    ierr = 0
    nzm_src = size( zm_src )
    nzt_src = size( zt_src )
    nzm_tgt = size( zm_tgt )
    nzt_tgt = size( zt_tgt )

    grid_ctx%grid_remap_method = grid_remap_method

    grid_ctx%gr_source%nzm = nzm_src
    grid_ctx%gr_source%nzt = nzt_src
    allocate( grid_ctx%gr_source%zm(ncol,nzm_src) )
    allocate( grid_ctx%gr_source%zt(ncol,nzt_src) )
    grid_ctx%gr_source%zm = spread( zm_src, dim=1, ncopies=ncol )
    grid_ctx%gr_source%zt = spread( zt_src, dim=1, ncopies=ncol )

    grid_ctx%gr_output%nzm = nzm_tgt
    grid_ctx%gr_output%nzt = nzt_tgt
    allocate( grid_ctx%gr_output%zm(ncol,nzm_tgt) )
    allocate( grid_ctx%gr_output%zt(ncol,nzt_tgt) )
    grid_ctx%gr_output%zm = spread( zm_tgt, dim=1, ncopies=ncol )
    grid_ctx%gr_output%zt = spread( zt_tgt, dim=1, ncopies=ncol )

    grid_ctx%is_initialized = .true.

  end subroutine stats_grid_init

  subroutine stats_update_grid( zt_src, zm_src, rho_vals, rho_levels, p_sfc, stats )

    ! Description:
    !   Refresh per-column remapping inputs for adaptive-grid runs. Updates the
    !   source grid altitudes, density-spline arrays, and surface pressure values
    !   that are consumed by profile remapping at write time.
    !-------------------------------------------------------------------------------

    ! ------------------------ Inputs ------------------------
    real(kind=core_rknd), dimension(:,:), intent(in) :: &
      zt_src, &     ! Per-column source thermodynamic heights [m]
      zm_src, &     ! Per-column source momentum heights [m]
      rho_vals, &   ! Density values for remap spline [kg m-3]
      rho_levels    ! Vertical coordinate for rho spline values [m]

    real(kind=core_rknd), dimension(:), intent(in) :: &
      p_sfc               ! Surface pressure by column [Pa]

    ! ---------------------- Input/Output ---------------------
    type(stats_type), intent(inout) :: &
      stats               ! Stats runtime context carrying adaptive-grid state [-]

    ! ------------------------- Locals ------------------------
    integer :: &
      nzt, &
      nzm, &
      nrho

    if ( .not. stats%enabled ) return
    if ( .not. stats%l_different_output_grid ) return
    if ( .not. stats%grid%is_initialized ) return

    if ( size( zt_src, 1 ) /= stats%ncol ) return
    if ( size( zm_src, 1 ) /= stats%ncol ) return
    if ( size( rho_vals, 1 ) /= stats%ncol ) return
    if ( size( rho_levels, 1 ) /= stats%ncol ) return
    if ( size( p_sfc ) /= stats%ncol ) return

    nzt = size( zt_src, 2 )
    nzm = size( zm_src, 2 )
    nrho = size( rho_vals, 2 )

    if ( nzt /= nzm - 1 ) return
    if ( size( rho_levels, 2 ) /= nrho ) return

    if ( stats%grid%gr_source%nzt /= nzt .or. stats%grid%gr_source%nzm /= nzm ) then
      if ( allocated( stats%grid%gr_source%zt ) ) deallocate( stats%grid%gr_source%zt )
      if ( allocated( stats%grid%gr_source%zm ) ) deallocate( stats%grid%gr_source%zm )
      stats%grid%gr_source%nzt = nzt
      stats%grid%gr_source%nzm = nzm
      allocate( stats%grid%gr_source%zt(stats%ncol,nzt) )
      allocate( stats%grid%gr_source%zm(stats%ncol,nzm) )
    end if

    if ( .not. allocated( stats%grid%rho_lin_spline_vals ) .or. &
        size( stats%grid%rho_lin_spline_vals, 2 ) /= nrho ) then
      if ( allocated( stats%grid%rho_lin_spline_vals ) ) deallocate( stats%grid%rho_lin_spline_vals )
      if ( allocated( stats%grid%rho_lin_spline_levels ) ) deallocate( stats%grid%rho_lin_spline_levels )
      allocate( stats%grid%rho_lin_spline_vals(stats%ncol,nrho) )
      allocate( stats%grid%rho_lin_spline_levels(stats%ncol,nrho) )
    end if

    if ( .not. allocated( stats%grid%p_sfc ) ) then
      allocate( stats%grid%p_sfc(stats%ncol) )
    end if

    stats%grid%gr_source%zt = zt_src
    stats%grid%gr_source%zm = zm_src
    stats%grid%rho_lin_spline_vals = rho_vals
    stats%grid%rho_lin_spline_levels = rho_levels
    stats%grid%p_sfc = p_sfc

  end subroutine stats_update_grid

  !----------------------------------------------------------------------------
  !                         SILHS samples handlers
  !----------------------------------------------------------------------------

  subroutine stats_lh_samples_init( num_samples, nzt, nl_var_names, u_var_names, zt_vals, &
                                    stats, err_info )

    ! Description:
    !   Define SILHS sample dimensions/variables in the main stats file.
    !   Creates lh_sample_number and lh_zt coordinates, registers requested
    !   lognormal/uniform variable names, and writes the lh_zt coordinate
    !   values used for subsequent SILHS sample output.
    !-------------------------------------------------------------------------------

    ! ------------------------ Inputs ------------------------
    integer, intent(in) :: &
      num_samples, &      ! Number of SILHS samples per column [-]
      nzt                 ! Number of SILHS vertical levels [-]

    character(len=*), dimension(:), intent(in) :: &
      nl_var_names, &     ! Names of lognormal SILHS outputs [-]
      u_var_names         ! Names of uniform SILHS outputs [-]

    real(kind=core_rknd), dimension(nzt), intent(in) :: &
      zt_vals             ! SILHS vertical coordinate values [m]

    ! ------------------------ InOuts ------------------------
    type(stats_type), intent(inout) :: &
      stats               ! Stats runtime context and NetCDF IDs [-]

    type(err_info_type), intent(inout) :: &
      err_info            ! Shared CLUBB error state container [-]

    ! ------------------------ Locals ------------------------
    integer :: ret_code, i, nnl, nu

    character(len=NAME_LEN) :: &
      var_name

    ! ------------------------ Begin Code ------------------------

    if ( .not. stats%enabled ) return
    if ( .not. stats%l_netcdf_output ) return
    if ( num_samples <= 0 .or. nzt <= 0 ) return
    if ( stats%zt_dimid < 0 ) return

    nnl = size( nl_var_names )
    nu = size( u_var_names )

    if ( stats%lh_2d%is_initialized ) then
      if ( stats%lh_2d%num_samples /= num_samples .or. stats%lh_2d%nzt /= nzt ) then
        write( fstderr,* ) err_info%err_header_global
        write( fstderr,* ) "stats lh init mismatch with existing definition"
        err_info%err_code = clubb_fatal_error
      end if
      return
    end if

    ret_code = nf90_redef( stats%ncid )

    ret_code = nf90_def_dim( stats%ncid, "lh_sample_number", num_samples, stats%lh_2d%sample_dimid )

    if ( stats%lh_2d%lh_zt_dimid < 0 ) then
      ret_code = nf90_def_dim( stats%ncid, "lh_zt", nzt, stats%lh_2d%lh_zt_dimid )
    end if

    if ( stats%lh_2d%lh_zt_varid < 0 ) then
      ret_code = nf90_def_var( stats%ncid, "lh_zt", NF90_DOUBLE, &
                          (/ stats%lh_2d%lh_zt_dimid /), stats%lh_2d%lh_zt_varid )
      ret_code = nf90_put_att( stats%ncid, stats%lh_2d%lh_zt_varid, "units", "m" )
    end if

    stats%lh_2d%num_samples = num_samples
    stats%lh_2d%nzt = nzt
    stats%lh_2d%n_nl_vars = nnl
    stats%lh_2d%n_u_vars = nu

    if ( nnl > 0 ) then
      allocate( stats%lh_2d%nl_varids(nnl) )
      allocate( stats%lh_2d%nl_names(nnl) )
      do i = 1, nnl

        var_name = "lh_nl_" // trim( adjustl( nl_var_names(i) ) )
        stats%lh_2d%nl_names(i) = var_name

        ret_code = nf90_def_var( stats%ncid, trim( var_name ), NF90_DOUBLE, &
                            (/ stats%col_dimid, stats%lh_2d%sample_dimid, &
                               stats%lh_2d%lh_zt_dimid, stats%time_dimid /), &
                            stats%lh_2d%nl_varids(i) )

        ret_code = nf90_def_var_fill( stats%ncid, stats%lh_2d%nl_varids(i), 0, 0.0_core_rknd )

      end do
    end if

    if ( nu > 0 ) then
      allocate( stats%lh_2d%u_varids(nu) )
      allocate( stats%lh_2d%u_names(nu) )

      do i = 1, nu

        var_name = "lh_u_" // trim( adjustl( u_var_names(i) ) )
        stats%lh_2d%u_names(i) = var_name

        ret_code = nf90_def_var( stats%ncid, trim( var_name ), NF90_DOUBLE, &
                            (/ stats%col_dimid, stats%lh_2d%sample_dimid, &
                               stats%lh_2d%lh_zt_dimid, stats%time_dimid /), &
                            stats%lh_2d%u_varids(i) )

        ret_code = nf90_def_var_fill( stats%ncid, stats%lh_2d%u_varids(i), 0, 0.0_core_rknd )
      end do

    end if

    ret_code = nf90_enddef( stats%ncid )

    ! Write
    ret_code = nf90_put_var( stats%ncid, stats%lh_2d%lh_zt_varid, zt_vals )

    stats%lh_2d%is_initialized = .true.

  end subroutine stats_lh_samples_init

  subroutine stats_lh_samples_write_lognormal( samples, stats, err_info )

    ! Description:
    !   Write one timestep of SILHS lognormal sample fields to NetCDF. Input
    !   samples are provided as (col, sample, zt, var) and are written directly
    !   to the stored variable layout (col, sample, zt, time) for each variable.
    !-------------------------------------------------------------------------------

    ! ------------------------ Inputs ------------------------
    real(kind=core_rknd), dimension(:,:,:,:), intent(in) :: &
      samples             ! Lognormal samples (col,sample,zt,var) [var-dependent]

    ! ------------------------ InOuts ------------------------
    type(stats_type), intent(inout) :: &
      stats               ! Stats runtime context and NetCDF IDs [-]

    type(err_info_type), intent(inout) :: &
      err_info            ! Shared CLUBB error state container [-]

    ! ------------------------ Locals ------------------------
    integer :: t, v, ret_code

    ! ------------------------ Begin Code ------------------------

    if ( .not. stats%enabled ) return
    if ( .not. stats%l_netcdf_output ) return
    if ( .not. stats%lh_2d%is_initialized ) return
    if ( .not. stats%l_last_sample ) return
    if ( stats%lh_2d%n_nl_vars <= 0 ) return

    if ( size( samples, 1 ) /= stats%ncol ) return
    if ( size( samples, 2 ) /= stats%lh_2d%num_samples ) return
    if ( size( samples, 3 ) /= stats%lh_2d%nzt ) return
    if ( size( samples, 4 ) /= stats%lh_2d%n_nl_vars ) return

    t = stats%time_index + 1
    do v = 1, stats%lh_2d%n_nl_vars

      ret_code = nf90_put_var( stats%ncid, stats%lh_2d%nl_varids(v), samples(:,:,:,v:v), &
                           start = (/ 1, 1, 1, t /), &
                           count = (/ stats%ncol, stats%lh_2d%num_samples, stats%lh_2d%nzt, 1 /) )

      if ( ret_code /= NF90_NOERR ) then
        write( fstderr,* ) err_info%err_header_global
        write( fstderr,* ) "stats lh lognormal write failed for ", trim( stats%lh_2d%nl_names(v) )
        err_info%err_code = clubb_fatal_error
        return
      end if

    end do

  end subroutine stats_lh_samples_write_lognormal

  subroutine stats_lh_samples_write_uniform( uniform_vals, mixture_comp, sample_weights, &
                                            stats, err_info )

    ! Description:
    !   Write one timestep of SILHS uniform sample outputs, including per-var
    !   uniform variates plus mixture-component and sample-weight fields. Input
    !   arrays are written to the NetCDF layout (col, sample, zt, time).
    !-------------------------------------------------------------------------------

    ! ------------------------ Inputs ------------------------
    real(kind=core_rknd), dimension(:,:,:,:), intent(in) :: &
      uniform_vals        ! Uniform variates (col,sample,zt,var) [0-1]

    integer, dimension(:,:,:), intent(in) :: &
      mixture_comp        ! Mixture-component index (col,sample,zt) [-]

    real(kind=core_rknd), dimension(:,:,:), intent(in) :: &
      sample_weights      ! Sample weight (col,sample,zt) [0-1]

    ! ------------------------ InOuts ------------------------
    type(stats_type), intent(inout) :: &
      stats               ! Stats runtime context and NetCDF IDs [-]

    type(err_info_type), intent(inout) :: &
      err_info            ! Shared CLUBB error state container [-]

    ! ------------------------ Locals ------------------------
    integer :: t, v, ret_code, dp2

    ! ------------------------ Begin Code ------------------------
      
    if ( .not. stats%enabled ) return
    if ( .not. stats%l_netcdf_output ) return
    if ( .not. stats%lh_2d%is_initialized ) return
    if ( .not. stats%l_last_sample ) return
    if ( stats%lh_2d%n_u_vars <= 0 ) return

    ! Sanity check for sizes
    if ( size( uniform_vals, 1 ) /= stats%ncol ) return
    if ( size( uniform_vals, 2 ) /= stats%lh_2d%num_samples ) return
    if ( size( uniform_vals, 3 ) /= stats%lh_2d%nzt ) return
    if ( size( mixture_comp, 1 ) /= stats%ncol ) return
    if ( size( mixture_comp, 2 ) /= stats%lh_2d%num_samples ) return
    if ( size( mixture_comp, 3 ) /= stats%lh_2d%nzt ) return
    if ( size( sample_weights, 1 ) /= stats%ncol ) return
    if ( size( sample_weights, 2 ) /= stats%lh_2d%num_samples ) return
    if ( size( sample_weights, 3 ) /= stats%lh_2d%nzt ) return

    dp2 = size( uniform_vals, 4 )
    if ( stats%lh_2d%n_u_vars /= dp2 + 2 ) return

    t = stats%time_index + 1
    do v = 1, dp2

      ! Write
      ret_code = nf90_put_var( stats%ncid, stats%lh_2d%u_varids(v), uniform_vals(:,:,:,v:v), &
                           start = (/ 1, 1, 1, t /), &
                           count = (/ stats%ncol, stats%lh_2d%num_samples, stats%lh_2d%nzt, 1 /) )
      
      if ( ret_code /= NF90_NOERR ) then
        write( fstderr,* ) err_info%err_header_global
        write( fstderr,* ) "stats lh uniform write failed for ", trim( stats%lh_2d%u_names(v) )
        err_info%err_code = clubb_fatal_error
        return
      end if

    end do

    ! Write
    ret_code = nf90_put_var( stats%ncid, stats%lh_2d%u_varids(dp2+1), &
                         spread( real( mixture_comp, kind=core_rknd ), dim = 4, ncopies = 1 ), &
                         start = (/ 1, 1, 1, t /), &
                         count = (/ stats%ncol, stats%lh_2d%num_samples, stats%lh_2d%nzt, 1 /) )

    if ( ret_code /= NF90_NOERR ) then
      write( fstderr,* ) err_info%err_header_global
      write( fstderr,* ) "stats lh uniform write failed for ", trim( stats%lh_2d%u_names(dp2+1) )
      err_info%err_code = clubb_fatal_error
      return
    end if

    ! Write
    ret_code = nf90_put_var( stats%ncid, stats%lh_2d%u_varids(dp2+2), &
                         spread( sample_weights, dim = 4, ncopies = 1 ), &
                         start = (/ 1, 1, 1, t /), &
                         count = (/ stats%ncol, stats%lh_2d%num_samples, stats%lh_2d%nzt, 1 /) )

    if ( ret_code /= NF90_NOERR ) then
      write( fstderr,* ) err_info%err_header_global
      write( fstderr,* ) "stats lh uniform write failed for ", trim( stats%lh_2d%u_names(dp2+2) )
      err_info%err_code = clubb_fatal_error
      return
    end if

  end subroutine stats_lh_samples_write_uniform

  !----------------------------------------------------------------------------
  !                         Helpers
  !----------------------------------------------------------------------------

  logical function var_on_stats_list(stats, name)

    ! Description:
    !   Return true when a variable name exists in the active registry. This is
    !   a read-only membership query that does not mutate lookup-cache state.
    !-------------------------------------------------------------------------------

    ! ------------------------ Inputs ------------------------
    type(stats_type), intent(in) :: &
      stats               ! Stats runtime context to query [-]
    character(len=*), intent(in) :: &
      name                ! Candidate variable name [-]

    ! ------------------------ Locals ------------------------
    integer :: i

    ! ------------------------ Begin Code ------------------------

    if ( .not. stats%enabled .or. len_trim( name ) == 0 ) then
      var_on_stats_list = .false.
      return
    end if
    
    do i = 1, stats%nvars
      if ( trim( stats%vars(i)%name ) == trim( adjustl( name ) ) ) then
        var_on_stats_list = .true.
        return
      end if
    end do

    var_on_stats_list = .false.

  end function var_on_stats_list
  
  integer function stats_get_id(stats, name) result(id)

    ! Description:
    !   Resolve a variable name to its internal index in stats%vars. The lookup
    !   first uses ordered positive/reject caches for fast-path hits and then
    !   falls back to linear search. Misses are recorded in the reject cache.
    !-------------------------------------------------------------------------------

    ! ---------------------- Input/Output ---------------------
    type(stats_type), intent(inout) :: &
      stats               ! Stats runtime context and lookup caches [-]

    ! ------------------------ Inputs ------------------------
    character(len=*), intent(in) :: &
      name                ! Variable name to resolve to internal ID [-]

    ! ------------------------- Locals ------------------------
    character(len=NAME_LEN) :: &
      needle

    integer :: &
      i, &
      cached_id

    ! ------------------------ Begin Code ------------------------

    id = 0
    if ( .not. stats%enabled .or. len_trim( name ) == 0 ) return

    needle = trim( adjustl( name ) )

    ! Fast path: expected lookup at current cache head.
    if ( stats%lookup%head <= stats%lookup%cache_len ) then
      cached_id = stats%lookup%cache(stats%lookup%head)
      if ( trim( stats%vars(cached_id)%name ) == needle ) then
        id = cached_id
        stats%lookup%head = stats%lookup%head + 1
        return
      end if
    end if

    !print *, "stats_get_id cache miss for ", trim(needle)

    ! Fast reject path: expected missing name at current reject-cache head.
    if ( stats%lookup%reject_head <= stats%lookup%reject_cache_len ) then
      if ( trim( stats%lookup%reject_cache(stats%lookup%reject_head) ) == needle ) then
        stats%lookup%reject_head = stats%lookup%reject_head + 1
        return
      end if
    end if

    !print *, "stats_get_id reject cache miss for ", trim(needle)

    ! Fallback reject path: repeated miss not in expected order.
    do i = 1, stats%lookup%reject_cache_len
      if ( trim( stats%lookup%reject_cache(i) ) == needle ) then
        return
      end if
    end do

    !print *, "stats_get_id reject cache full miss for ", trim(needle)

    ! Fallback: linear search.
    do i = 1, stats%nvars
      if ( trim( stats%vars(i)%name ) == needle ) then
        id = i
        if ( stats%lookup%head > stats%lookup%cache_len ) then
          call stats_lookup_append( stats, id )
          stats%lookup%head = stats%lookup%head + 1
        else if ( stats%lookup%cache(stats%lookup%head) == id ) then
          stats%lookup%head = stats%lookup%head + 1
        end if
        return
      end if
    end do

    !print *, "stats_get_id linear search miss for ", trim(needle), " -- adding to reject cache"

    call stats_lookup_reject_append( stats, needle )

  end function stats_get_id

  subroutine stats_lookup_append( stats, id )

    ! Description:
    !   Append a successful variable-id lookup to the ordered positive cache,
    !   growing the cache allocation when needed.
    !-------------------------------------------------------------------------------

    ! ------------------------ Inputs ------------------------
    integer, intent(in) :: &
      id                  ! Resolved variable ID to append to hit cache [-]

    ! ------------------------ InOuts ------------------------
    type(stats_type), intent(inout) :: &
      stats               ! Stats runtime context and hit cache [-]

    ! ------------------------ Locals ------------------------
    integer, allocatable, dimension(:) :: tmp
    integer :: new_cap

    ! ------------------------ Begin Code ------------------------

    ! Positive cache stores the expected lookup order within a timestep.
    if ( .not. allocated( stats%lookup%cache ) ) then
      ! Case 1: first use this run; allocate an initial cache block.
      new_cap = max( 16, stats%nvars )
      allocate( stats%lookup%cache(new_cap) )
      stats%lookup%cache = 0
      stats%lookup%cache_len = 0
    else if ( stats%lookup%cache_len >= size( stats%lookup%cache ) ) then
      ! Case 2: cache is full; grow capacity geometrically to keep appends cheap.
      new_cap = max( 16, 2 * size( stats%lookup%cache ) )
      allocate( tmp(new_cap) )
      tmp = 0
      if ( stats%lookup%cache_len > 0 ) tmp(1:stats%lookup%cache_len) = stats%lookup%cache(1:stats%lookup%cache_len)
      call move_alloc( tmp, stats%lookup%cache )
    end if

    ! Append newly confirmed id at the tail to preserve observed order.
    stats%lookup%cache_len = stats%lookup%cache_len + 1
    stats%lookup%cache(stats%lookup%cache_len) = id
  end subroutine stats_lookup_append

  subroutine stats_lookup_reject_append( stats, name )

    ! Description:
    !   Append a failed variable-name lookup to the ordered reject cache,
    !   growing the cache allocation when needed.
    !-------------------------------------------------------------------------------

    ! ------------------------ Inputs ------------------------
    character(len=*), intent(in) :: &
      name                ! Missing variable name to append to miss cache [-]

    ! ------------------------ InOuts ------------------------
    type(stats_type), intent(inout) :: &
      stats               ! Stats runtime context and miss cache [-]

    ! ------------------------ Locals ------------------------
    character(len=NAME_LEN), allocatable, dimension(:) :: tmp
    integer :: new_cap

    ! ------------------------ Begin Code ------------------------

    if ( len_trim( name ) == 0 ) return

    ! Reject cache tracks repeated misses to avoid redundant linear scans.
    if ( .not. allocated( stats%lookup%reject_cache ) ) then
      ! Case 1: first reject entry; allocate the initial reject-cache block.
      new_cap = max( 16, stats%nvars )
      allocate( stats%lookup%reject_cache(new_cap) )
      stats%lookup%reject_cache = ''
      stats%lookup%reject_cache_len = 0
    else if ( stats%lookup%reject_cache_len >= size( stats%lookup%reject_cache ) ) then
      ! Case 2: reject cache is full; grow geometrically and preserve entries.
      new_cap = max( 16, 2 * size( stats%lookup%reject_cache ) )
      allocate( tmp(new_cap) )
      tmp = ''
      if ( stats%lookup%reject_cache_len > 0 ) then
        tmp(1:stats%lookup%reject_cache_len) = stats%lookup%reject_cache(1:stats%lookup%reject_cache_len)
      end if
      call move_alloc( tmp, stats%lookup%reject_cache )
    end if

    ! Store normalized names so reject comparisons are stable.
    stats%lookup%reject_cache_len = stats%lookup%reject_cache_len + 1
    stats%lookup%reject_cache(stats%lookup%reject_cache_len) = trim( adjustl( name ) )

  end subroutine stats_lookup_reject_append

  subroutine format_date &
             ( day_in, month_in, year_in, time_in, &
               date )

    ! Description:
    !   Convert model start date/time inputs into a UDUNITS-compatible
    !   reference-time string used for the NetCDF time variable units.
    !
    ! Notes:
    !   Adapted from the original GrADS version written by Chris Golaz.
    !   Uses Fortran `internal' files to write the string output.
    !-------------------------------------------------------------------------------

    use calendar, only: &
        compute_current_date_api ! Procedure(s)

    use clubb_precision, only: &
        time_precision ! Variable(s)

    implicit none

    ! External
    intrinsic :: floor, int, mod, nint

    ! Input Variables
    integer, intent(in) ::  & 
      day_in,           & ! Day of Month at Model Start   [dd]
      month_in,         & ! Month of Year at Model Start  [mm]
      year_in             ! Year at Model Start         [yyyy]

    real(kind=time_precision), intent(in) :: time_in ! Start time [s]

    ! Output Variables
    character(len=35), intent(out) :: date ! UDUNITS time-units string ("seconds since ...") [-]

    integer::  & 
      iday, imonth, iyear  ! Integer for day, month and year.

    real(kind=time_precision) :: st_time ! Start time [s]

    ! ------------------------ Begin Code ------------------------

    call compute_current_date_api( day_in, month_in,  & ! intent(in)
                                   year_in, &  ! intent(in)
                                   time_in, & ! intent(in)
                                   iday, imonth, & ! intent(out)
                                   iyear, &  ! intent(out)
                                   st_time ) ! intent(out)

    date = "seconds since YYYY-MM-DD HH:MM:00.0"
    write(date(15:18),'(i4.4)') iyear
    write(date(20:21),'(i2.2)') imonth
    write(date(23:24),'(i2.2)') iday
    write(date(26:27),'(i2.2)') floor( st_time / 3600._time_precision )
    write(date(29:30),'(i2.2)') int( mod( nint( st_time ),3600 ) / 60 )
    
    write(date(32:33),'(i2.2)') nint( ( ( real( mod( nint( st_time ),3600 ),kind=time_precision ) / &
                    60._time_precision ) - ( real( int( mod( nint( st_time ),3600 ) / 60 ), & 
                                              kind=time_precision ) ) )*60._time_precision )

    return
  end subroutine format_date

  !----------------------------------------------------------------------------
  !                         Namelist parsing
  !----------------------------------------------------------------------------

  subroutine stats_read_registry_namelist( path, defs, ndefs, ierr )

    ! Description:
    !   Read registry entries from `&clubb_stats_nl` with `entry(i)=...`,
    !   then filter empty/comment entries, allocate the registry definition
    !   array, and parse each active line into a stats_var definition structure.
    !-------------------------------------------------------------------------------

    ! ------------------------ Inputs ------------------------
    character(len=*), intent(in) :: &
      path                ! Registry namelist file path [-]

    ! ------------------------ Outputs ------------------------
    type(stats_var), allocatable, dimension(:), intent(out) :: &
      defs                ! Parsed registry definitions [-]
    integer, intent(out) :: &
      ndefs, &            ! Number of parsed active definitions [-]
      ierr                ! Nonzero on file/read/parse failure [-]

    ! ------------------------ Locals ------------------------
    integer :: iunit, ios, i, count
    logical :: l_parse_ok
    character(len=REG_LINE_LEN) , dimension(NML_REG_MAX_ENTRIES) :: entry
    namelist /clubb_stats_nl/ entry

    ! ------------------------ Begin Code ------------------------

    ierr = 0
    ndefs = 0
    entry(:) = ''

    open( newunit=iunit, file=trim( path ), status='old', action='read', iostat=ios )
    if ( ios /= 0 ) then
      ierr = ios
      return
    end if

    read( iunit, nml=clubb_stats_nl, iostat=ios )

    close( iunit )
    if ( ios /= 0 ) then
      ierr = ios
      return
    end if

    count = 0
    do i = 1, NML_REG_MAX_ENTRIES
      if ( len_trim( entry(i) ) == 0 ) cycle
      if ( is_line_comment( entry(i) ) ) cycle
      count = count + 1
    end do

    if ( count <= 0 ) return

    allocate( defs(count) )
    ndefs = 0
    do i = 1, NML_REG_MAX_ENTRIES

      if ( len_trim( entry(i) ) == 0 ) cycle
      if ( is_line_comment( entry(i) ) ) cycle
      ndefs = ndefs + 1

      call parse_registry_line( entry(i), defs(ndefs), l_parse_ok )

      if ( .not. l_parse_ok ) then
        ierr = -1
        write( fstderr,* ) "Invalid stats registry entry (expected: name | grid | units | long_name)"
        write( fstderr,* ) trim( entry(i) )
        return
      end if

    end do

  end subroutine stats_read_registry_namelist

  subroutine stats_expand_registry( sclr_dim, edsclr_dim, defs, ndefs, hydromet_list )

    ! Description:
    !   Expand shorthand/template registry names (e.g., hydrometeor/scalar
    !   families) into explicit variable definitions. Deduplicates names while
    !   preserving order and returns a compact definitions array.
    !-------------------------------------------------------------------------------

    ! ------------------------ Inputs ------------------------
    integer, intent(in) :: &
      sclr_dim, &         ! Number of scalar tracers for template expansion [-]
      edsclr_dim          ! Number of ED scalar tracers for expansion [-]

    ! ------------------------ InOuts ------------------------
    type(stats_var), allocatable, dimension(:), intent(inout) :: &
      defs                ! Definition array to expand in-place [-]
    integer, intent(inout) :: &
      ndefs               ! Input/output count of active definitions [-]

    ! ------------------------ Optional Ins ------------------------
    character(len=*), dimension(:), intent(in), optional :: &
      hydromet_list       ! Hydrometeor names list [-]

    ! ------------------------ Locals ------------------------
    type(stats_var), allocatable, dimension(:) :: out_defs
    type(stats_var) :: base_def
    character(len=16) :: idx_str
    integer :: max_defs
    integer :: i, j, j2, k, nhm

    ! ------------------------ Begin Code ------------------------

    nhm = 0
    if ( present( hydromet_list ) ) nhm = size( hydromet_list )

    ! Over-allocate once, then trim at the end.
    max_defs = max( 1, ndefs + 4096 )
    allocate( out_defs(max_defs) )

    k = 0
    do i = 1, ndefs

      base_def = defs(i)

      select case ( trim( base_def%name ) )
      case ( "hm_i" )
        do j = 1, nhm
          call add_expanded_def(base_def,trim(hydromet_list(j)( 1:2 ))//"_1",out_defs,k)
          call add_expanded_def(base_def,trim(hydromet_list(j)( 1:2 ))//"_2",out_defs,k)
        end do
      case ( "mu_hm_i" )
        do j = 1, nhm
          call add_expanded_def(base_def,"mu_"//trim(hydromet_list(j)( 1:2 ))//"_1",out_defs,k)
          call add_expanded_def(base_def,"mu_"//trim(hydromet_list(j)( 1:2 ))//"_2",out_defs,k)
        end do
      case ( "mu_Ncn_i" )
        call add_expanded_def( base_def, "mu_Ncn_1",out_defs, k )
        call add_expanded_def( base_def, "mu_Ncn_2",out_defs, k )
      case ( "mu_hm_i_n" )
        do j = 1, nhm
          call add_expanded_def(base_def,"mu_"//trim(hydromet_list(j)( 1:2 ))//"_1_n",out_defs,k)
          call add_expanded_def(base_def,"mu_"//trim(hydromet_list(j)( 1:2 ))//"_2_n",out_defs,k)
        end do
      case ( "mu_Ncn_i_n" )
        call add_expanded_def( base_def, "mu_Ncn_1_n",out_defs, k )
        call add_expanded_def( base_def, "mu_Ncn_2_n",out_defs, k )
      case ( "sigma_hm_i" )
        do j = 1, nhm
          call add_expanded_def(base_def,"sigma_"//trim(hydromet_list(j)( 1:2 ))//"_1",out_defs,k)
          call add_expanded_def(base_def,"sigma_"//trim(hydromet_list(j)( 1:2 ))//"_2",out_defs,k)
        end do
      case ( "sigma_Ncn_i" )
        call add_expanded_def( base_def, "sigma_Ncn_1",out_defs, k )
        call add_expanded_def( base_def, "sigma_Ncn_2",out_defs, k )
      case ( "sigma_hm_i_n" )
        do j = 1, nhm
          call add_expanded_def(base_def,"sigma_"//trim(hydromet_list(j)( 1:2 ))//"_1_n",out_defs,k)
          call add_expanded_def(base_def,"sigma_"//trim(hydromet_list(j)( 1:2 ))//"_2_n",out_defs,k)
        end do
      case ( "sigma_Ncn_i_n" )
        call add_expanded_def( base_def, "sigma_Ncn_1_n",out_defs, k )
        call add_expanded_def( base_def, "sigma_Ncn_2_n",out_defs, k )
      case ( "corr_w_hm_i" )
        do j = 1, nhm
          call add_expanded_def(base_def,"corr_w_"//trim(hydromet_list(j)( 1:2 ))//"_1",out_defs,k)
          call add_expanded_def(base_def,"corr_w_"//trim(hydromet_list(j)( 1:2 ))//"_2",out_defs,k)
        end do
      case ( "corr_w_Ncn_i" )
        call add_expanded_def( base_def, "corr_w_Ncn_1",out_defs, k )
        call add_expanded_def( base_def, "corr_w_Ncn_2",out_defs, k )
      case ( "corr_chi_hm_i" )
        do j = 1, nhm
          call add_expanded_def(base_def,"corr_chi_"//trim(hydromet_list(j)( 1:2 ))//"_1",out_defs,k)
          call add_expanded_def(base_def,"corr_chi_"//trim(hydromet_list(j)( 1:2 ))//"_2",out_defs,k)
        end do
      case ( "corr_chi_Ncn_i" )
        call add_expanded_def( base_def, "corr_chi_Ncn_1",out_defs, k )
        call add_expanded_def( base_def, "corr_chi_Ncn_2",out_defs, k )
      case ( "corr_eta_hm_i" )
        do j = 1, nhm
          call add_expanded_def(base_def,"corr_eta_"//trim(hydromet_list(j)( 1:2 ))//"_1",out_defs,k)
          call add_expanded_def(base_def,"corr_eta_"//trim(hydromet_list(j)( 1:2 ))//"_2",out_defs,k)
        end do
      case ( "corr_eta_Ncn_i" )
        call add_expanded_def( base_def, "corr_eta_Ncn_1",out_defs, k )
        call add_expanded_def( base_def, "corr_eta_Ncn_2",out_defs, k )
      case ( "corr_Ncn_hm_i" )
        do j = 1, nhm
          call add_expanded_def(base_def,"corr_Ncn_"//trim(hydromet_list(j)( 1:2 ))//"_1",out_defs,k)
          call add_expanded_def(base_def,"corr_Ncn_"//trim(hydromet_list(j)( 1:2 ))//"_2",out_defs,k)
        end do
      case ( "corr_hmx_hmy_i" )
        do j = 1, nhm
          do j2 = j + 1, nhm
            call add_expanded_def( base_def, &
                                   "corr_"//trim( hydromet_list(j)( 1:2 ) )//"_"// &
                                   trim( hydromet_list(j2)( 1:min( 2,len_trim( hydromet_list(j2) ) ) ) )//"_1", &
                                   out_defs, k )
            call add_expanded_def( base_def, &
                                   "corr_"//trim( hydromet_list(j)( 1:2 ) )//"_"// &
                                   trim( hydromet_list(j2)( 1:min( 2,len_trim( hydromet_list(j2) ) ) ) )//"_2", &
                                   out_defs, k )
          end do
        end do
      case ( "corr_w_hm_i_n" )
        do j = 1, nhm
          call add_expanded_def(base_def,"corr_w_"//trim(hydromet_list(j)( 1:2 ))//"_1_n",out_defs,k)
          call add_expanded_def(base_def,"corr_w_"//trim(hydromet_list(j)( 1:2 ))//"_2_n",out_defs,k)
        end do
      case ( "corr_w_Ncn_i_n" )
        call add_expanded_def( base_def, "corr_w_Ncn_1_n",out_defs, k )
        call add_expanded_def( base_def, "corr_w_Ncn_2_n",out_defs, k )
      case ( "corr_chi_hm_i_n" )
        do j = 1, nhm
          call add_expanded_def(base_def,"corr_chi_"//trim(hydromet_list(j)( 1:2 ))//"_1_n",out_defs,k)
          call add_expanded_def(base_def,"corr_chi_"//trim(hydromet_list(j)( 1:2 ))//"_2_n",out_defs,k)
        end do
      case ( "corr_chi_Ncn_i_n" )
        call add_expanded_def( base_def, "corr_chi_Ncn_1_n",out_defs, k )
        call add_expanded_def( base_def, "corr_chi_Ncn_2_n",out_defs, k )
      case ( "corr_eta_hm_i_n" )
        do j = 1, nhm
          call add_expanded_def(base_def,"corr_eta_"//trim(hydromet_list(j)( 1:2 ))//"_1_n",out_defs,k)
          call add_expanded_def(base_def,"corr_eta_"//trim(hydromet_list(j)( 1:2 ))//"_2_n",out_defs,k)
        end do
      case ( "corr_eta_Ncn_i_n" )
        call add_expanded_def( base_def, "corr_eta_Ncn_1_n",out_defs, k )
        call add_expanded_def( base_def, "corr_eta_Ncn_2_n",out_defs, k )
      case ( "corr_Ncn_hm_i_n" )
        do j = 1, nhm
          call add_expanded_def(base_def,"corr_Ncn_"//trim(hydromet_list(j)( 1:2 ))//"_1_n",out_defs,k)
          call add_expanded_def(base_def,"corr_Ncn_"//trim(hydromet_list(j)( 1:2 ))//"_2_n",out_defs,k)
        end do
      case ( "corr_hmx_hmy_i_n" )
        do j = 1, nhm
          do j2 = j + 1, nhm
            call add_expanded_def( base_def, &
                                   "corr_"//trim( hydromet_list(j)( 1:2 ) )//"_"// &
                                   trim( hydromet_list(j2)( 1:min( 2,len_trim( hydromet_list(j2) ) ) ) )//"_1_n", &
                                   out_defs, k )
            call add_expanded_def( base_def, &
                                   "corr_"//trim( hydromet_list(j)( 1:2 ) )//"_"// &
                                   trim( hydromet_list(j2)( 1:min( 2,len_trim( hydromet_list(j2) ) ) ) )//"_2_n", &
                                   out_defs, k )
          end do
        end do
      case ( "wp2hmp" )
        do j = 1, nhm
          call add_expanded_def(base_def,"wp2"//trim(hydromet_list(j)( 1:2 ))//"p",out_defs,k)
        end do
      case ( "K_hm" )
        do j = 1, nhm
          call add_expanded_def(base_def,"K_hm_"//trim(hydromet_list(j)( 1:2 )),out_defs,k)
        end do
      case ( "hydrometp2" )
        do j = 1, nhm
          call add_expanded_def(base_def,trim(hydromet_list(j)( 1:2 ))//"p2",out_defs,k)
        end do
      case ( "wphydrometp" )
        do j = 1, nhm
          call add_expanded_def(base_def,"wp"//trim(hydromet_list(j)( 1:2 ))//"p",out_defs,k)
        end do
      case ( "hmp2_zt" )
        do j = 1, nhm
          call add_expanded_def(base_def,trim(hydromet_list(j)( 1:2 ))//"p2_zt",out_defs,k)
        end do
      case ( "rtphmp" )
        do j = 1, nhm
          call add_expanded_def(base_def,"rtp"//trim(hydromet_list(j)( 1:2 ))//"p",out_defs,k)
        end do
      case ( "thlphmp" )
        do j = 1, nhm
          call add_expanded_def(base_def,"thlp"//trim(hydromet_list(j)( 1:2 ))//"p",out_defs,k)
        end do
      case ( "hmxphmyp" )
        do j = 1, nhm
          do j2 = j + 1, nhm
            call add_expanded_def( base_def, &
                                   trim( hydromet_list(j)( 1:2 ) )//"p"// &
                                   trim( hydromet_list(j2)( 1:min( 2,len_trim( hydromet_list(j2) ) ) ) )//"p", &
                                   out_defs, k )
          end do
        end do
      case ( "sclrm" )
        do j = 1, sclr_dim
          write( idx_str,'(I0)' ) j
          call add_expanded_def( base_def, "sclr"//trim( idx_str )//"m",out_defs, k )
        end do
      case ( "sclrm_f" )
        do j = 1, sclr_dim
          write( idx_str,'(I0)' ) j
          call add_expanded_def( base_def, "sclr"//trim( idx_str )//"m_f",out_defs, k )
        end do
      case ( "edsclrm" )
        do j = 1, edsclr_dim
          write( idx_str,'(I0)' ) j
          call add_expanded_def( base_def, "edsclr"//trim( idx_str )//"m",out_defs, k )
        end do
      case ( "edsclrm_f" )
        do j = 1, edsclr_dim
          write( idx_str,'(I0)' ) j
          call add_expanded_def( base_def, "edsclr"//trim( idx_str )//"m_f",out_defs, k )
        end do
      case ( "sclrprtp" )
        do j = 1, sclr_dim
          write( idx_str,'(I0)' ) j
          call add_expanded_def( base_def, "sclr"//trim( idx_str )//"prtp",out_defs, k )
        end do
      case ( "sclrp2" )
        do j = 1, sclr_dim
          write( idx_str,'(I0)' ) j
          call add_expanded_def( base_def, "sclr"//trim( idx_str )//"p2",out_defs, k )
        end do
      case ( "sclrpthvp" )
        do j = 1, sclr_dim
          write( idx_str,'(I0)' ) j
          call add_expanded_def( base_def, "sclr"//trim( idx_str )//"pthvp",out_defs, k )
        end do
      case ( "sclrpthlp" )
        do j = 1, sclr_dim
          write( idx_str,'(I0)' ) j
          call add_expanded_def( base_def, "sclr"//trim( idx_str )//"pthlp",out_defs, k )
        end do
      case ( "sclrprcp" )
        do j = 1, sclr_dim
          write( idx_str,'(I0)' ) j
          call add_expanded_def( base_def, "sclr"//trim( idx_str )//"prcp",out_defs, k )
        end do
      case ( "wpsclrp" )
        do j = 1, sclr_dim
          write( idx_str,'(I0)' ) j
          call add_expanded_def( base_def, "wpsclr"//trim( idx_str )//"p",out_defs, k )
        end do
      case ( "wpsclrp2" )
        do j = 1, sclr_dim
          write( idx_str,'(I0)' ) j
          call add_expanded_def( base_def, "wpsclr"//trim( idx_str )//"p2",out_defs, k )
        end do
      case ( "wp2sclrp" )
        do j = 1, sclr_dim
          write( idx_str,'(I0)' ) j
          call add_expanded_def( base_def, "wp2sclr"//trim( idx_str )//"p",out_defs, k )
        end do
      case ( "wpsclrprtp" )
        do j = 1, sclr_dim
          write( idx_str,'(I0)' ) j
          call add_expanded_def( base_def, "wpsclr"//trim( idx_str )//"prtp",out_defs, k )
        end do
      case ( "wpsclrpthlp" )
        do j = 1, sclr_dim
          write( idx_str,'(I0)' ) j
          call add_expanded_def( base_def, "wpsclr"//trim( idx_str )//"pthlp",out_defs, k )
        end do
      case ( "silhs_variance_category" )
        do j = 1, 8
          write( idx_str,'(I0)' ) j
          call add_expanded_def( base_def, "silhs_var_cat_"//trim( idx_str ),out_defs, k )
        end do
      case ( "lh_samp_frac_category" )
        do j = 1, 8
          write( idx_str,'(I0)' ) j
          call add_expanded_def( base_def, "lh_samp_frac_"//trim( idx_str ),out_defs, k )
        end do
      case ( "wpedsclrp" )
        do j = 1, edsclr_dim
          write( idx_str,'(I0)' ) j
          call add_expanded_def( base_def, "wpedsclr"//trim( idx_str )//"p",out_defs, k )
        end do
      case default
        call add_expanded_def( base_def, base_def%name,out_defs, k )
      end select
    end do

    deallocate( defs )
    allocate( defs( max( 1, k ) ) )
    if ( k > 0 ) defs(1:k) = out_defs(1:k)
    ndefs = k
    deallocate( out_defs )

  end subroutine stats_expand_registry

  subroutine add_expanded_def( base_def, name, out_defs, nout )

    ! Description:
    !   Append one expanded variable definition to an output definition array
    !   unless a definition with the same name is already present.
    !-------------------------------------------------------------------------------

    ! ------------------------ Inputs ------------------------
    type(stats_var), intent(in) :: &
      base_def            ! Template/base definition to copy from [-]
    character(len=*), intent(in) :: &
      name                ! Expanded variable name to add [-]

    ! ------------------------ InOuts ------------------------
    type(stats_var), dimension(:), intent(inout) :: &
      out_defs            ! Destination expanded definitions buffer [-]
    integer, intent(inout) :: &
      nout                ! Current number of used entries in out_defs [-]

    ! ------------------------ Locals ------------------------
    type(stats_var) :: def
    integer :: i
    character(len=NAME_LEN) :: needle
    logical, save :: l_warned_overflow = .false.

    needle = trim( name )
    if ( len_trim( needle ) == 0 ) return

    do i = 1, nout
      if ( trim( out_defs(i)%name ) == needle ) return
    end do

    if ( nout >= size( out_defs ) ) then
      if ( .not. l_warned_overflow ) then
        write( fstderr,* ) "stats registry expansion overflow; some template vars may be skipped"
        l_warned_overflow = .true.
      end if
      return
    end if

    def = base_def
    def%name = needle
    def%varid = -1
    def%nz = 0

    nout = nout + 1
    out_defs(nout) = def
    
  end subroutine add_expanded_def

  subroutine parse_registry_line( line, def, l_valid )

    ! Description:
    !   Parse one registry entry string of the form
    !   "name | grid | units | long_name" into a stats_var definition and
    !   map the textual grid token to an internal grid-id enumerator.
    !-------------------------------------------------------------------------------

    ! ------------------------ Inputs ------------------------
    character(len=*), intent(in) :: &
      line                ! Raw registry entry line to parse [-]

    ! ------------------------ Outputs ------------------------
    type(stats_var), intent(out) :: &
      def                 ! Parsed registry definition result [-]
    logical, intent(out) :: &
      l_valid             ! True if line parsed into a valid definition [-]

    ! ------------------------ Locals ------------------------
    character(len=256) , dimension(4) :: fields
    integer :: nfields

    ! ------------------------ Begin Code ------------------------

    def%name = ''
    def%grid = ''
    def%grid_id = GRID_ZT
    def%units = ''
    def%long_name = ''
    def%varid = -1
    def%nz = 0
    def%out_nz = 0
    def%l_budget = .false.
    def%l_in_budget = .false.
    fields = ''
    l_valid = .false.
    call split_registry_fields( line, fields, nfields )

    if ( nfields /= 3 .and. nfields /= 4 ) return

    def%name = trim( adjustl( fields(1) ) )
    def%grid = trim( adjustl( fields(2) ) )
    if ( nfields == 4 ) then
      def%units = trim( adjustl( fields(3) ) )
      def%long_name = trim( adjustl( fields(4) ) )
    else
      ! Legacy 3-field compatibility: name | grid | long_name
      def%long_name = trim( adjustl( fields(3) ) )
    end if
    if ( len_trim( def%name ) == 0 .or. len_trim( def%grid ) == 0 ) return

    select case ( trim( def%grid ) )
    case ( "zt", "nzt" )
      def%grid_id  = GRID_ZT
    case ( "lh_zt" )
      def%grid_id  = GRID_LH_ZT
    case ( "zm", "lh_zm", "nzm" )
      def%grid_id  = GRID_ZM
    case ( "rad_zt" )
      def%grid_id  = GRID_RAD_ZT
    case ( "rad_zm" )
      def%grid_id  = GRID_RAD_ZM
    case ( "sfc" )
      def%grid_id  = GRID_SFC
    case ( "lh_sfc" )
      def%grid_id  = GRID_LH_SFC
    case default
      def%grid_id  = GRID_ZT
    end select

    l_valid = .true.

  end subroutine parse_registry_line

  subroutine split_registry_fields( line, fields, nfields )

    ! Description:
    !   Split a pipe-delimited registry line into fields and trim whitespace
    !   around each field.
    !-------------------------------------------------------------------------------

    ! ------------------------ Inputs ------------------------
    character(len=*), intent(in) :: &
      line                ! Raw registry line to split [-]

    ! ------------------------ InOuts ------------------------
    character(len=*), dimension(:), intent(inout) :: &
      fields              ! Output field buffer (trimmed tokens) [-]

    ! ------------------------ Outputs ------------------------
    integer, intent(out) :: &
      nfields             ! Number of parsed fields from line [-]

    ! ------------------------ Locals ------------------------
    integer :: i, lenline, istart
    character :: ch

    ! ------------------------ Begin Code ------------------------

    nfields = 0
    fields = ''
    lenline = len_trim( line )
    if ( lenline <= 0 ) return

    istart = 1
    do i = 1, lenline
      ch = line(i:i)

      if ( ch == '|' ) then
        nfields = nfields + 1
        if ( nfields <= size( fields ) ) then
          if ( i > istart ) then
            fields(nfields) = trim(adjustl(line( istart:i-1 )))
          else
            fields(nfields) = ''
          end if
        end if
        istart = i + 1
      end if
    end do

    nfields = nfields + 1
    if ( nfields <= size( fields ) ) then
      if ( lenline >= istart ) then
        fields(nfields) = trim(adjustl(line( istart:lenline )))
      else
        fields(nfields) = ''
      end if
    end if

  end subroutine split_registry_fields

  logical function is_line_comment(line)

    ! Description:
    !   Return true for blank lines and lines beginning with comment markers
    !   (`#` or `!`) used in registry namelist entries.
    !-------------------------------------------------------------------------------

    ! ------------------------ Inputs ------------------------
    character(len=*), intent(in) :: &
      line                ! Candidate registry line to classify [-]

    ! ------------------------ Locals ------------------------
    character(len=1) :: c

    ! ------------------------ Begin Code ------------------------
    if ( len_trim( line ) == 0 ) then
      is_line_comment = .true.
      return
    end if
    c = line(1:1)
    is_line_comment = ( c == '#' .or. c == '!' )
  end function is_line_comment

end module stats_netcdf
