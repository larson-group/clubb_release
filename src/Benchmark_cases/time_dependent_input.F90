!$Id$
module time_dependent_input
!
!  Description:
!    This module is responsible for managing the reading in and
!    storage of time dependent information for a case.
!
!  References:
!    None
!--------------------------------------------------------------------------------------------------

  use constants_clubb, only: eps, one

  use clubb_precision, only: &
    core_rknd ! Variable(s)

  implicit none

  public :: finalize_t_dependent_input, time_select, &
            apply_time_dependent_forcings, &
            apply_time_dependent_forcings_from_dycore, &
            initialize_t_dependent_input

  private :: initialize_t_dependent_forcings, &
             finalize_t_dependent_forcings,   & 
             initialize_t_dependent_sfc,  &
             finalize_t_dependent_sfc,    &
             read_to_grid, &
             apply_time_dependent_forcings_from_array                      

  integer, parameter :: nforcings = 10 ! Number of columns in the input file

  ! Module variables used to describe 
  real( kind = core_rknd ), public, allocatable, dimension(:) :: &
    time_sfc_given, &                                  ! the surface over time.
    latent_ht_given, &
    sens_ht_given, &
    thlm_sfc_given, &
    rtm_sfc_given,  &
    CO2_sfc_given,  &
    upwp_sfc_given, &
    vpwp_sfc_given, &
    T_sfc_given, &
    wpthlp_sfc_given, &
    wpqtp_sfc_given

!$omp threadprivate( time_sfc_given, latent_ht_given, sens_ht_given, thlm_sfc_given, &
!$omp   rtm_sfc_given, CO2_sfc_given,  upwp_sfc_given, vpwp_sfc_given, &
!$omp   T_sfc_given, wpthlp_sfc_given, wpqtp_sfc_given )

  real( kind = core_rknd ), private, allocatable :: &
    t_dependent_forcing_data(:,:,:,:), & ! (Column, level, variable, time) [units vary]
    t_dependent_forcing_data_f_grid(:,:,:), & ! Source (level, variable, time) [units vary]
    time_forcing_given(:) ! Shared forcing times [s]

  character(len=30), private :: forcing_names(nforcings) ! Source field names and units

!$omp threadprivate( t_dependent_forcing_data, t_dependent_forcing_data_f_grid, &
!$omp                time_forcing_given, forcing_names )

  logical, public :: l_t_dependent ! Flag used to determine when
  !                                  time dependent information is read in.
  !                                  It is suggested that the flag be checked
  !                                  before using any of the variables stored
  !                                  in the module.
!$omp threadprivate( l_t_dependent )

  logical, public :: l_input_xpwp_sfc ! Flag used to determine whether or not to read 
                                      ! in the surface momentum fluxes, upwp_sfc and vpwp_sfc.
!$omp threadprivate( l_input_xpwp_sfc )

  logical, public :: l_ignore_forcings ! Flag used to determine if the forcings
                                       ! should be ignored for this case.
!$omp threadprivate( l_ignore_forcings )

  logical, public :: l_sfc_already_initialized = .false. ! Flag used to determine if the sfc
                                                         ! variable were already initialized
!$omp threadprivate( l_sfc_already_initialized )

  private

  contains

  !================================================================================================
  subroutine initialize_t_dependent_input( ngrdcol, iunit, & ! In
                                           forcing_file, sfc_file, iop_file, & ! In
                                           l_readiopdata, rho, & ! In
                                           nz, & ! In
                                           grid, p_in_Pa, & ! In
                                           l_add_dycore_grid, & ! In
                                           grid_adapt_in_time_method ) ! In
    !
    !  Description: 
    !    This subroutine reads in time dependent information about a
    !    case that is stored inside the module.
    !
    !  References:
    !    None
    !---------------------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! Input Variable(s)
    integer, intent(in) :: ngrdcol ! Number of model columns [-]
    integer, intent(in) :: iunit ! File I/O

    character(len=*), intent(in) :: &
      forcing_file, & ! Time-dependent forcing profiles
      sfc_file, & ! Time-dependent surface data
      iop_file ! Optional IOP NetCDF dataset

    integer, intent(in) :: nz ! Size of the model grid

    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: grid ! Model grid

    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: p_in_Pa ! Pressure[Pa]

    logical, intent(in) :: &
      l_add_dycore_grid, &     ! Flag to set remapping from dycore on or off
      l_readiopdata            ! Flag whether or not we are using a netcdf file for input

    real(kind=core_rknd), dimension(:,:), intent(in) :: rho

    integer, intent(in) :: &
      grid_adapt_in_time_method   ! flag that says whether the grid should be adapted over
                                  ! time and if so, how the grid density function should
                                  ! be formed

    ! ----------------- Begin Code --------------------

    if ( .not. l_ignore_forcings ) then
      call initialize_t_dependent_forcings &
                   ( ngrdcol, iunit, forcing_file, nz, grid, & ! In
                     p_in_Pa, & ! In
                     l_add_dycore_grid, & ! In
                     grid_adapt_in_time_method ) ! In
    end if

    if ( .not. l_sfc_already_initialized ) then
      call initialize_t_dependent_sfc &
                    ( iunit, l_readiopdata, rho, sfc_file, iop_file ) ! Intent(in)
    end if

  end subroutine initialize_t_dependent_input

  !================================================================================================
  subroutine finalize_t_dependent_input()
    !
    ! Description: 
    !   This subroutine frees memory stored after initilizing the
    !   time dependent data of this module.
    !
    ! References:
    !   None
    !-----------------------------------------------------------------------------

    implicit none

    ! ----------------- Begin Code --------------------

    if ( .not. l_ignore_forcings ) then
      call finalize_t_dependent_forcings()
    end if

    call finalize_t_dependent_sfc()

  end subroutine finalize_t_dependent_input

  !================================================================================================
  subroutine initialize_t_dependent_sfc( iunit, l_readiopdata, rho, input_file, iop_file )
    !
    !  Description: This subroutine reads in a file that details time dependent
    !  input values that vary in one dimension.
    !-----------------------------------------------------------------------------

    use input_reader, only: &
      read_one_dim_file, one_dim_read_var, & ! Procedure(s)
      fill_blanks_one_dim_vars, read_x_profile, &
      get_target_index, deallocate_one_dim_vars, &
      count_columns

    use readiopdata_module, only: readiopdata_sfc

    use sfc_flux, only: &
      convert_sens_ht_to_km_s, &
      convert_latent_ht_to_m_s

    use input_names, only: &
      time_name,     &
      thetal_name,   &
      rt_name,       &
      latent_ht_name,       &
      sens_ht_name,       &
      CO2_umol_name, &
      upwp_sfc_name, &
      vpwp_sfc_name, &
      T_sfc_name,    &
      wpthlp_sfc_name, &
      wpqtp_sfc_name

    implicit none

    ! Input Variable(s)
    integer, intent(in) :: iunit ! File I/O unit

    character(len=*), intent(in) :: input_file ! Path to surface input file
    character(len=*), intent(in) :: iop_file ! Optional IOP NetCDF dataset

    logical, intent(in) :: l_readiopdata
    real(kind=core_rknd), dimension(:,:), intent(in) :: rho

    ! Local Variable(s)

    integer :: i

    type(one_dim_read_var), allocatable, dimension(:) :: &
      retVars ! retVars stores the name of a variable (e.g. pressure),
              ! the name of the dimension the variable varies along (e.g. time), and
              ! the time-dependent values of the variable to be input into CLUBB.

    integer ::  &
      dim_size, & ! Number of time-dependent values of a variable to be input into CLUBB
      nforcings       ! Number of variables with time-dependent input data


    ! ----------------- Begin Code --------------------
    if ( l_readiopdata ) then
      call readiopdata_sfc( iop_file, & ! In
                            T_sfc_given, shsfc=sens_ht_given, & ! Optional Out
                            lhsfc=latent_ht_given,            & ! Intent(out)
                            timesfc=time_sfc_given )            ! Intent(out)
      allocate( wpthlp_sfc_given(size(time_sfc_given)), wpqtp_sfc_given(size(time_sfc_given)) )
      do i=1,size(time_sfc_given)
        wpthlp_sfc_given(i) = convert_sens_ht_to_km_s( sens_ht_given(i), rho(1,1) )   ! Intent(in)
        wpqtp_sfc_given(i) = convert_latent_ht_to_m_s( latent_ht_given(i), rho(1,1) ) ! Intent(in)
      end do
    else
      nforcings = count_columns( iunit, input_file )

      allocate( retVars(1:nforcings) )

      ! Read the surface.in file and store the necessary input information in retVars
      call read_one_dim_file( iunit, nforcings, input_file, retVars )

      ! Fill blank values stored as -999.9 using linear interpolation
      call fill_blanks_one_dim_vars( nforcings, retVars )

      ! dim_size is the number of values input for a particular variable
      dim_size = size( retVars(1)%values )

      ! Store the data read from the file in each [variable]_sfc_given

      if( get_target_index(nforcings, time_name, retVars) > 0 ) then
        allocate( time_sfc_given(1:dim_size) )
        time_sfc_given = read_x_profile( nforcings, dim_size, time_name, retVars, &
                                      input_file )
      end if

      if( get_target_index(nforcings, latent_ht_name, retVars) > 0 ) then
        allocate( latent_ht_given(1:dim_size) )
        latent_ht_given = read_x_profile( nforcings, dim_size, latent_ht_name, retVars, &
                                input_file )
      end if

      if( get_target_index(nforcings, sens_ht_name, retVars) > 0 ) then
        allocate( sens_ht_given(1:dim_size) )
        sens_ht_given = read_x_profile( nforcings, dim_size, sens_ht_name, retVars, &
                                input_file )
      end if

      if( get_target_index(nforcings, thetal_name, retVars) > 0 ) then
        allocate( thlm_sfc_given(1:dim_size) )
        thlm_sfc_given = read_x_profile( nforcings, dim_size, thetal_name, retVars, &
                                      input_file )
      end if

      if( get_target_index(nforcings, rt_name, retVars) > 0 ) then
        allocate( rtm_sfc_given(1:dim_size) )
        rtm_sfc_given = read_x_profile( nforcings, dim_size, rt_name, retVars, &
                                      input_file )
      end if

      ! As of July 2010, this is only in cobra
      if( get_target_index(nforcings, CO2_umol_name, retVars) > 0 ) then
        allocate( CO2_sfc_given(1:dim_size) )
        CO2_sfc_given = read_x_profile( nforcings, dim_size, CO2_umol_name, retVars, &
                                        input_file )
      end if

      ! As of July 2010, this is only in gabls3_night
      if( get_target_index(nforcings, upwp_sfc_name, retVars) > 0 ) then
        allocate( upwp_sfc_given(1:dim_size) )
        upwp_sfc_given = read_x_profile( nforcings, dim_size, upwp_sfc_name, retVars, &
                                        input_file )
      end if

      ! As of July 2010, this is only in gabls3_night
      if( get_target_index(nforcings, vpwp_sfc_name, retVars) > 0 ) then
        allocate( vpwp_sfc_given(1:dim_size) )
        vpwp_sfc_given = read_x_profile( nforcings, dim_size, vpwp_sfc_name, retVars, &
                                        input_file )
      end if

      ! As of July 2010, this is only in astex_a209
      if( get_target_index(nforcings, T_sfc_name, retVars) > 0 ) then
        allocate( T_sfc_given(1:dim_size) )
        T_sfc_given = read_x_profile( nforcings, dim_size, T_sfc_name, retVars, &
                                        input_file )
      end if


      if( get_target_index(nforcings, wpthlp_sfc_name, retVars) > 0 ) then
        allocate( wpthlp_sfc_given(1:dim_size) )
        wpthlp_sfc_given = read_x_profile( nforcings, dim_size, wpthlp_sfc_name, &
                                        retVars, input_file )
      end if

      if( get_target_index(nforcings, wpqtp_sfc_name, retVars) > 0 ) then
        allocate( wpqtp_sfc_given(1:dim_size) )
        wpqtp_sfc_given = read_x_profile( nforcings, dim_size, wpqtp_sfc_name, &
                                        retVars, input_file )
      end if

      ! Deallocate memory
      call deallocate_one_dim_vars( nforcings, retVars )
    end if

    l_sfc_already_initialized = .true.

    return
  end subroutine initialize_t_dependent_sfc

  !================================================================================================
  subroutine initialize_t_dependent_forcings( ngrdcol, iunit, input_file, nz, grid, & ! In
                                              p_in_Pa, & ! In
                                              l_add_dycore_grid, & ! In
                                              grid_adapt_in_time_method ) ! In
    !
    !  Description: This subroutine reads in a file that details time dependent
    !  input values that vary in two dimensions.
    !
    !-------------------------------------------------------------------------------------

    use input_reader, only: read_two_dim_file, fill_blanks_two_dim_vars, &
      two_dim_read_var, one_dim_read_var

    use input_names, only: &
      z_name, &
      pressure_name

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    use model_flags, only: &
      no_grid_adaptation

    implicit none

    ! External
    intrinsic :: spread

    ! Input Variable(s)
    integer, intent(in) :: ngrdcol ! Number of model columns [-]
    integer, intent(in) :: iunit ! File I/O

    character(len=*), intent(in) :: input_file ! Path to the input file

    integer, intent(in) :: nz  ! Size of Model Grid [-]

    real( kind = core_rknd ), intent(in), dimension(ngrdcol,nz) :: grid ! Altitudes of Grid [m]

    real( kind = core_rknd ), intent(in), dimension(ngrdcol,nz) :: p_in_Pa ! Pressure [Pa]

    logical, intent(in) :: &
      l_add_dycore_grid     ! Flag to set remapping from dycore on or off

    integer, intent(in) :: &
      grid_adapt_in_time_method   ! flag that says whether the grid should be adapted over
                                  ! time and if so, how the grid density function should
                                  ! be formed

    ! Local Variables

    integer :: i, n, source_index ! Column and forcing indices [-]
    integer :: n_f_grid_z, n_f_grid_t ! Source level and time counts [-]

    type(two_dim_read_var) :: forcing_input(nforcings) ! Temporary source-file records [units vary]
    type(one_dim_read_var) :: forcing_times ! Temporary source time coordinate [s]

    ! ----------------- Begin Code --------------------

    if ( .not. allocated(t_dependent_forcing_data_f_grid) ) then

      call read_two_dim_file( iunit, nforcings, input_file, & ! In
                              forcing_input, forcing_times ) ! Out

      ! Fill source-file gaps before packing the temporary records into arrays.
      call fill_blanks_two_dim_vars( nforcings, forcing_times, forcing_input )

      forcing_names = forcing_input%name
      time_forcing_given = forcing_times%values
      allocate( t_dependent_forcing_data_f_grid(size(forcing_input(1)%values,1), &
                                                nforcings,size(time_forcing_given)) )
      do n = 1, nforcings

        t_dependent_forcing_data_f_grid(:,n,:) = forcing_input(n)%values

      end do

    end if

    n_f_grid_z = size(t_dependent_forcing_data_f_grid,1)
    n_f_grid_t = size(time_forcing_given)

    if ( allocated(t_dependent_forcing_data) ) then

      if ( any(shape(t_dependent_forcing_data) /= [ngrdcol,nz,nforcings,n_f_grid_t]) ) then
        !$acc exit data delete( t_dependent_forcing_data )
        deallocate( t_dependent_forcing_data )
      end if

    end if

    if ( .not. allocated(t_dependent_forcing_data) ) then
      allocate( t_dependent_forcing_data(ngrdcol,nz,nforcings,n_f_grid_t) )
      !$acc enter data create( t_dependent_forcing_data )
    end if

    ! Preserve each column's interpolation coordinates on the model grid.
    do i = 1, ngrdcol

      select case ( forcing_names(1) )
      case ( z_name )
        t_dependent_forcing_data(i,:,1,:) = spread(grid(i,:),2,n_f_grid_t)
      case ( pressure_name )
        t_dependent_forcing_data(i,:,1,:) = spread(-p_in_Pa(i,:),2,n_f_grid_t)
      case default
        error stop "Incompatible grid type in first element of t_dependent_forcings."
      end select

    end do

    if ( forcing_names(1) == pressure_name ) then
      ! Retained raw pressure coordinates must stay negative across grid updates.
      t_dependent_forcing_data_f_grid(:,1,:) = -abs(t_dependent_forcing_data_f_grid(:,1,:))
    end if

    do n = 2, nforcings

      source_index = findloc( forcing_names, forcing_names(n), dim=1 )
      do i = 1, ngrdcol

        t_dependent_forcing_data(i,:,n,:) = read_to_grid( n_f_grid_z, n_f_grid_t, nz, & ! In
          t_dependent_forcing_data(i,:,1,1), t_dependent_forcing_data_f_grid(:,1,:), & ! In
          t_dependent_forcing_data_f_grid(:,source_index,:) ) ! In

      end do

    end do

    ! Refresh the device table after initialization or a model-grid change.
    !$acc update device( t_dependent_forcing_data )

    ! Grid adaptation needs the source profiles again when the model grid changes.
    if ( l_add_dycore_grid .or. grid_adapt_in_time_method == no_grid_adaptation ) then
      deallocate( t_dependent_forcing_data_f_grid )
    end if

    return
  end subroutine initialize_t_dependent_forcings


  !================================================================================================
  subroutine finalize_t_dependent_forcings()
    ! Description:
    ! Release the source forcing table, model-grid profiles and shared times.

    implicit none

    !--------------------------- Begin Code ---------------------------
    if ( allocated(t_dependent_forcing_data_f_grid) ) then
      deallocate( t_dependent_forcing_data_f_grid )
    end if

    if ( allocated(t_dependent_forcing_data) ) then

      !$acc exit data delete( t_dependent_forcing_data )
      deallocate( t_dependent_forcing_data )

    end if
    if ( allocated(time_forcing_given) ) deallocate( time_forcing_given )

  end subroutine finalize_t_dependent_forcings
  
  !================================================================================================
  subroutine finalize_t_dependent_sfc( )
    !
    !  Description: Clears memory initialized in initialize_t_dependent_surface.
    !  This should be called at the end of the model.
    !
    !------------------------------------------------------------------------------------

    implicit none

    ! ----------------- Begin Code --------------------

    if ( allocated( time_sfc_given ) ) deallocate( time_sfc_given )
    if ( allocated( latent_ht_given ) )       deallocate( latent_ht_given )
    if ( allocated( sens_ht_given ) )       deallocate( sens_ht_given )
    if ( allocated( thlm_sfc_given ) ) deallocate( thlm_sfc_given )
    if ( allocated( rtm_sfc_given ) )  deallocate( rtm_sfc_given )
    if ( allocated( CO2_sfc_given ) )  deallocate( CO2_sfc_given )
    if ( allocated( upwp_sfc_given ) ) deallocate( upwp_sfc_given )
    if ( allocated( vpwp_sfc_given ) ) deallocate( vpwp_sfc_given )
    if ( allocated( T_sfc_given ) )    deallocate( T_sfc_given )
    if ( allocated( wpthlp_sfc_given ) ) deallocate( wpthlp_sfc_given )
    if ( allocated( wpqtp_sfc_given ) )  deallocate( wpqtp_sfc_given )

  end subroutine finalize_t_dependent_sfc

  !================================================================================================
  function read_to_grid( dim_size, other_dim_size, nz, grid, grid_src, var_src ) result(var)
    ! Description:
    ! Interpolate a source field onto one column's model grid at every input time.
    ! Independent time samples form the batch passed to vertical interpolation.

    use interpolation, only: zlinterp

    implicit none

    !--------------------------- Input Variables ---------------------------
    integer, intent(in) :: dim_size ! Source level count [-]
    integer, intent(in) :: other_dim_size ! Time sample count [-]
    integer, intent(in) :: nz ! Model level count [-]
    real( kind = core_rknd ), intent(in) :: grid(nz) ! Model heights [m] or negative pressures [Pa]
    real( kind = core_rknd ), intent(in) :: &
      grid_src(dim_size,other_dim_size), & ! Source heights [m] or negative pressures [Pa]
      var_src(dim_size,other_dim_size) ! Source field [units vary]

    !--------------------------- Output Variables ---------------------------
    real( kind = core_rknd ) :: var(nz,other_dim_size) ! Interpolated field [units vary]

    !--------------------------- Begin Code ---------------------------
    var = transpose( zlinterp( other_dim_size, nz, dim_size, & ! In
                               spread(grid,1,other_dim_size), & ! In
                               transpose(grid_src), transpose(var_src) ) ) ! In

  end function read_to_grid

  !================================================================================================
  subroutine apply_time_dependent_forcings_from_array( &
              ngrdcol, nzm, nzt, &
              sclr_dim, edsclr_dim, sclr_idx, &
              l_do_expldiff_rtm_thlm, &
              gr, rtm, rho, exner,  &
              forcings_array, &
              thlm_f, rtm_f, um_ref, vm_ref, um_f, vm_f, &
              wm_zt, wm_zm,  ug, vg, &
              sclrm_forcing, edsclrm_forcing )
    !
    !  Description: This subroutine is a helper subroutine, that takes the forcings
    !               in a 3D array (one set of forcing profiles per column) and
    !               initializes the forcing variables.
    !
    !---------------------------------------------------------------------------------

    use constants_clubb, only: &
      grav, & ! Variable(s)
      sec_per_hr, &
      pascal_per_mb, &
      
      zero, &
      fstderr

    use input_names, only: &
      temperature_f_name, &  ! Variable(s)
      rt_f_name,&
      sp_humidity_f_name, &
      thetal_f_name, &
      theta_f_name, &
      wm_name, &
      omega_name, &
      um_ref_name, &
      vm_ref_name, &
      um_f_name, &
      vm_f_name, &
      ug_name,&
      vg_name, &
      omega_mb_hr_name

    use grid_class, only : &
      grid,     & ! Type
      zt2zm_api   ! Procedure(s)

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    use array_index, only: &
      sclr_idx_type

    implicit none

    !--------------------- Input Variables ---------------------
    integer, intent(in) :: &
      ngrdcol, &
      nzm, &
      nzt, &
      sclr_dim, & 
      edsclr_dim

    type (sclr_idx_type), intent(in) :: &
      sclr_idx

    logical, intent(in) :: &
      l_do_expldiff_rtm_thlm ! thlm/rtm eddy scalars are copies diffused in
                             ! advance_windm_edsclrm

    type (grid), intent(in) :: &
      gr

    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
      exner,   & ! Exner Function                             [-]
      rho,     & ! Air Density                                [kg/m^3]
      rtm        ! Total Water Mixing Ratio                   [kg/kg]

    real( kind = core_rknd ), dimension(ngrdcol,nzt,nforcings), intent(in) :: &
      forcings_array ! Forcing profiles (column, level, variable) [units vary]

    !--------------------- Output Variables ---------------------
    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(inout) :: &
      thlm_f, & ! Potential Temperature forcing     [K/s]
      rtm_f,  & ! Total Water Mixing Ration forcing [kg/kg/s]
      um_ref, & ! um reference                      [m/s]
      vm_ref, & ! vm reference                      [m/s]
      um_f,   & ! um tendency                       [m/s/s]
      vm_f,   & ! vm tendency                       [m/s/s]
      wm_zt,  & ! subsidence on zt grid             [m/s]
      ug,     & ! u geostrophic wind                [m/s]
      vg        ! v geostrophic wind                [m/s]

    real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(inout) :: &
      wm_zm     ! subsidense on zm grid             [m/s]

    real( kind = core_rknd ), dimension(ngrdcol,nzt, sclr_dim), intent(inout) :: &
      sclrm_forcing ! Scalar forcing [-]

    real( kind = core_rknd ), dimension(ngrdcol,nzt, edsclr_dim), intent(inout) :: &
      edsclrm_forcing ! Edscalar forcing [-]

    !--------------------- Local Variables ---------------------
    logical :: l_valid(ngrdcol) ! Columns with a complete forcing profile [-]
    logical :: l_column_valid ! Current column has no missing input levels [-]
    real( kind = core_rknd ) :: wm_zm_new(ngrdcol,nzm) ! Interpolated subsidence [m/s]

    integer :: i, k, n ! Column, level and forcing indices [-]

    real( kind = core_rknd ), dimension(ngrdcol,nzt) :: temp_array

    !--------------------- Begin Code ---------------------

    !$acc data copyin( forcings_array ) create( temp_array, l_valid, wm_zm_new )

    ! Apply each named forcing field to its corresponding model array.
    do n=2, nforcings

      ! Each device thread copies and checks one column's complete profile.
      !$acc parallel loop gang vector default(present) private(l_column_valid)
      do i = 1, ngrdcol

        l_column_valid = .true.
        !$acc loop seq
        do k = 1, nzt

          temp_array(i,k) = forcings_array(i,k,n)
          if ( abs(temp_array(i,k) - (-999.9_core_rknd)) < &
               abs(temp_array(i,k) + (-999.9_core_rknd)) / 2 * eps ) then
            l_column_valid = .false.
          end if

        end do
        l_valid(i) = l_column_valid

      end do

      select case (forcing_names(n))
      case(temperature_f_name, theta_f_name, thetal_f_name)

        select case(forcing_names(n))
        case(temperature_f_name)

          !$acc parallel loop gang vector collapse(2) default(present)
          do k = 1, nzt
            do i = 1, ngrdcol
              if ( .not. l_valid(i) ) cycle
              thlm_f(i,k) = temp_array(i,k) / exner(i,k)
            end do
          end do

        case(theta_f_name)

          !$acc parallel loop gang vector collapse(2) default(present)
          do k = 1, nzt
            do i = 1, ngrdcol
              if ( .not. l_valid(i) ) cycle
              thlm_f(i,k) = temp_array(i,k) ! n am not sure on the conversion of this
            end do
          end do

        case(thetal_f_name)

          !$acc parallel loop gang vector collapse(2) default(present)
          do k = 1, nzt
            do i = 1, ngrdcol
              if ( .not. l_valid(i) ) cycle
              thlm_f(i,k) = temp_array(i,k)
            end do
          end do

        end select

        if ( sclr_idx%iisclr_thl > 0 ) then

          !$acc parallel loop gang vector collapse(2) default(present)
          do k = 1, nzt
            do i = 1, ngrdcol
              if ( .not. l_valid(i) ) cycle
              sclrm_forcing(i,k,sclr_idx%iisclr_thl) = thlm_f(i,k)
            end do
          end do

        end if

        if ( sclr_idx%iiedsclr_thl > 0 ) then

          ! With l_do_expldiff_rtm_thlm, advance_windm_edsclrm puts a copy of thlm
          ! in this eddy scalar. thlm is already forced with thlm_f in
          ! advance_xm_wpxp, so the copy must not be forced again.
          !$acc parallel loop gang vector collapse(2) default(present)
          do k = 1, nzt
            do i = 1, ngrdcol
              if ( .not. l_valid(i) ) cycle
              if ( l_do_expldiff_rtm_thlm .and. edsclr_dim > 1 ) then
                edsclrm_forcing(i,k,sclr_idx%iiedsclr_thl) = zero
              else
                edsclrm_forcing(i,k,sclr_idx%iiedsclr_thl) = thlm_f(i,k)
              end if
            end do
          end do

        end if

      case(rt_f_name, sp_humidity_f_name)

        select case(forcing_names(n))
        case(sp_humidity_f_name)

          !$acc parallel loop gang vector collapse(2) default(present)
          do k = 1, nzt
            do i = 1, ngrdcol
              if ( .not. l_valid(i) ) cycle
              rtm_f(i,k) = temp_array(i,k) * ( 1._core_rknd + rtm(i,k) )**2
            end do
          end do

        case(rt_f_name )

          !$acc parallel loop gang vector collapse(2) default(present)
          do k = 1, nzt
            do i = 1, ngrdcol
              if ( .not. l_valid(i) ) cycle
              rtm_f(i,k) = temp_array(i,k)
            end do
          end do

        end select

        if ( sclr_idx%iisclr_rt  > 0 ) then

          !$acc parallel loop gang vector collapse(2) default(present)
          do k = 1, nzt
            do i = 1, ngrdcol
              if ( .not. l_valid(i) ) cycle
              sclrm_forcing(i,k,sclr_idx%iisclr_rt)  = rtm_f(i,k)
            end do
          end do

        end if

        if ( sclr_idx%iiedsclr_rt  > 0 ) then

          ! Copy of rtm with l_do_expldiff_rtm_thlm; see the thlm case above.
          !$acc parallel loop gang vector collapse(2) default(present)
          do k = 1, nzt
            do i = 1, ngrdcol
              if ( .not. l_valid(i) ) cycle
              if ( l_do_expldiff_rtm_thlm .and. edsclr_dim > 1 ) then
                edsclrm_forcing(i,k,sclr_idx%iiedsclr_rt)  = zero
              else
                edsclrm_forcing(i,k,sclr_idx%iiedsclr_rt)  = rtm_f(i,k)
              end if
            end do
          end do

        end if

      case(um_ref_name)

        !$acc parallel loop gang vector collapse(2) default(present)
        do k = 1, nzt
          do i = 1, ngrdcol
            if ( .not. l_valid(i) ) cycle
            um_ref(i,k) = temp_array(i,k)
          end do
        end do

      case(vm_ref_name)

        !$acc parallel loop gang vector collapse(2) default(present)
        do k = 1, nzt
          do i = 1, ngrdcol
            if ( .not. l_valid(i) ) cycle
            vm_ref(i,k) = temp_array(i,k)
          end do
        end do

      case(um_f_name)

        !$acc parallel loop gang vector collapse(2) default(present)
        do k = 1, nzt
          do i = 1, ngrdcol
            if ( .not. l_valid(i) ) cycle
            um_f(i,k) = temp_array(i,k)
          end do
        end do

      case(vm_f_name)

        !$acc parallel loop gang vector collapse(2) default(present)
        do k = 1, nzt
          do i = 1, ngrdcol
            if ( .not. l_valid(i) ) cycle
            vm_f(i,k) = temp_array(i,k)
          end do
        end do

      case(wm_name, omega_name, omega_mb_hr_name)

        select case(forcing_names(n))
        case(wm_name)

          !$acc parallel loop gang vector collapse(2) default(present)
          do k = 1, nzt
            do i = 1, ngrdcol
              if ( .not. l_valid(i) ) cycle
              wm_zt(i,k) = temp_array(i,k)
            end do
          end do

        case(omega_name)

          !$acc parallel loop gang vector collapse(2) default(present)
          do k = 1, nzt
            do i = 1, ngrdcol
              if ( .not. l_valid(i) ) cycle
              wm_zt(i,k) = - temp_array(i,k) / (grav * rho(i,k))
            end do
          end do

        case(omega_mb_hr_name)

          !$acc parallel loop gang vector collapse(2) default(present)
          do k = 1, nzt
            do i = 1, ngrdcol
              if ( .not. l_valid(i) ) cycle
              temp_array(i,k) = temp_array(i,k) * pascal_per_mb / sec_per_hr
            end do
          end do

          !$acc parallel loop gang vector collapse(2) default(present)
          do k = 1, nzt
            do i = 1, ngrdcol
              if ( .not. l_valid(i) ) cycle

              wm_zt(i,k) = - temp_array(i,k) / (grav * rho(i,k))

            end do
          end do

        end select

        wm_zm_new(:,:) = zt2zm_api( nzm, nzt, ngrdcol, gr, wm_zt(:,:) )

        !$acc parallel loop gang vector collapse(2) default(present)
        do k = 1, nzm
          do i = 1, ngrdcol
            if ( l_valid(i) ) wm_zm(i,k) = wm_zm_new(i,k)
          end do
        end do

      case(ug_name)

        !$acc parallel loop gang vector collapse(2) default(present)
        do k = 1, nzt
          do i = 1, ngrdcol
            if ( .not. l_valid(i) ) cycle
            ug(i,k) = temp_array(i,k)
          end do
        end do

      case(vg_name)

        !$acc parallel loop gang vector collapse(2) default(present)
        do k = 1, nzt
          do i = 1, ngrdcol
            if ( .not. l_valid(i) ) cycle
            vg(i,k) = temp_array(i,k)
          end do
        end do

      case default

        ! Unknown fields are errors only when at least one profile is valid.
        !$acc update host( l_valid )
        if ( any(l_valid) ) then
          write(fstderr, *) "Incompatable forcing type: "//forcing_names(n)
          error stop
        end if

      end select


    end do ! 1 .. nforcings

    !$acc end data

    return

  end subroutine apply_time_dependent_forcings_from_array

  !================================================================================================
  subroutine apply_time_dependent_forcings( &
              ngrdcol, nzm, nzt, &
              sclr_dim, edsclr_dim, sclr_idx, &
              l_do_expldiff_rtm_thlm, &
              gr, time, rtm, rho, exner,  &
              thlm_f, rtm_f, um_ref, vm_ref, um_f, vm_f, &
              wm_zt, wm_zm,  ug, vg, &
              sclrm_forcing, edsclrm_forcing )
    !
    !  Description: This subroutine converts the time dependent information stored in
    !  memory (time_dependent_forcing_data) into the format used by CLUBB.
    !
    !---------------------------------------------------------------------------------

    use interpolation, only: &
      linear_interp_factor ! Procedure(s)

    use grid_class, only : &
      grid,     & ! Type
      zt2zm_api   ! Procedure(s)

    use clubb_precision, only: &
      time_precision, &
      core_rknd ! Variable(s)

    use array_index, only: &
      sclr_idx_type

    implicit none

    !--------------------- Input Variables ---------------------
    integer, intent(in) :: &
      ngrdcol, &
      nzm, &
      nzt, &
      sclr_dim, & 
      edsclr_dim

    type (sclr_idx_type), intent(in) :: &
      sclr_idx

    logical, intent(in) :: &
      l_do_expldiff_rtm_thlm ! thlm/rtm eddy scalars are copies diffused in
                             ! advance_windm_edsclrm

    type (grid), intent(in) :: &
      gr

    real(kind=time_precision), intent(in) :: &
      time ! Model Time [s]

    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
      exner,   & ! Exner Function                             [-]
      rho,     & ! Air Density                                [kg/m^3]
      rtm        ! Total Water Mixing Ratio                   [kg/kg]

    !--------------------- Output Variables ---------------------
    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(inout) :: &
      thlm_f, & ! Potential Temperature forcing     [K/s]
      rtm_f,  & ! Total Water Mixing Ration forcing [kg/kg/s]
      um_ref, & ! um reference                      [m/s]
      vm_ref, & ! vm reference                      [m/s]
      um_f,   & ! um tendency                       [m/s/s]
      vm_f,   & ! vm tendency                       [m/s/s]
      wm_zt,  & ! subsidence on zt grid             [m/s]
      ug,     & ! u geostrophic wind                [m/s]
      vg        ! v geostrophic wind                [m/s]

    real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(inout) :: &
      wm_zm     ! subsidense on zm grid             [m/s]

    real( kind = core_rknd ), dimension(ngrdcol,nzt, sclr_dim), intent(inout) :: &
      sclrm_forcing ! Scalar forcing [-]

    real( kind = core_rknd ), dimension(ngrdcol,nzt, edsclr_dim), intent(inout) :: &
      edsclrm_forcing ! Edscalar forcing [-]

    !--------------------- Local Variables ---------------------
    integer :: i, k, n ! Column, level and forcing indices [-]
    integer :: before_time, after_time ! Bracketing forcing time indices [-]

    real( kind = core_rknd ), dimension(ngrdcol,nzt,nforcings) :: &
      forcings_array ! Time-interpolated fields (column, level, variable) [units vary]

    real( kind = core_rknd ) :: time_frac ! Fraction between bracketing input times [-]

    !--------------------- Begin Code ---------------------

    time_frac = -one ! Default initialization

    call time_select( time, size(time_forcing_given), time_forcing_given, &
                                 before_time, after_time, time_frac )

    !$acc data create( forcings_array )

    !$acc parallel loop gang vector collapse(3) default(present)
    do n = 2, nforcings
      do k = 1, nzt
        do i = 1, ngrdcol
          forcings_array(i,k,n) = linear_interp_factor &
            ( time_frac, t_dependent_forcing_data(i,k,n,after_time), &
              t_dependent_forcing_data(i,k,n,before_time) )
        end do
      end do
    end do

    call apply_time_dependent_forcings_from_array( &
                                                   ngrdcol, nzm, nzt, &
                                                   sclr_dim, edsclr_dim, sclr_idx, &
                                                   l_do_expldiff_rtm_thlm, &
                                                   gr, rtm, rho, exner,  &
                                                   forcings_array, &
                                                   thlm_f, rtm_f, um_ref, vm_ref, um_f, vm_f, &
                                                   wm_zt, wm_zm,  ug, vg, &
                                                   sclrm_forcing, edsclrm_forcing )

    !$acc end data

    return

  end subroutine apply_time_dependent_forcings

  !================================================================================================
  subroutine apply_time_dependent_forcings_from_dycore( &
              ngrdcol, nzm, nzt, &
              sclr_dim, edsclr_dim, sclr_idx, &
              l_do_expldiff_rtm_thlm, &
              gr, gr_dycore, time, rtm, rho, exner, &
              grid_remap_method, &
              total_idx_rho_lin_spline, rho_lin_spline_vals, &
              rho_lin_spline_levels, &
              p_sfc, &
              thlm_f, rtm_f, um_ref, vm_ref, um_f, vm_f, &
              wm_zt, wm_zm,  ug, vg, &
              sclrm_forcing, edsclrm_forcing )
    !
    !  Description: This subroutine converts the time dependent information stored in
    !  memory (time_dependent_forcing_data) which was read in to the dycore grid into the format
    !  used by CLUBB and remaps it to the current physics grid.
    !
    !---------------------------------------------------------------------------------

    use interpolation, only: &
      linear_interp_factor ! Procedure(s)

    use grid_class, only : &
        grid,     & ! Type
        zt2zm_api   ! Procedure(s)

    use clubb_precision, only: &
      time_precision, &
      core_rknd ! Variable(s)

    use array_index, only: &
        sclr_idx_type

    use remapping_module, only: &
        remap_vals_to_target

    implicit none

    !--------------------- Input Variables ---------------------
    integer, intent(in) :: &
      ngrdcol, &
      nzm, &
      nzt, &
      sclr_dim, & 
      edsclr_dim

    type (sclr_idx_type), intent(in) :: &
      sclr_idx

    logical, intent(in) :: &
      l_do_expldiff_rtm_thlm ! thlm/rtm eddy scalars are copies diffused in
                             ! advance_windm_edsclrm

    type (grid), intent(in) :: &
      gr, gr_dycore

    real(kind=time_precision), intent(in) :: &
      time ! Model Time [s]

    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
      exner,   & ! Exner Function                             [-]
      rho,     & ! Air Density                                [kg/m^3]
      rtm        ! Total Water Mixing Ratio                   [kg/kg]

    integer, intent(in) :: &
      total_idx_rho_lin_spline, & ! number of indices for the linear spline definition arrays
      grid_remap_method  ! Integer that specifies what remapping method should be used

    real( kind = core_rknd ), dimension(ngrdcol,total_idx_rho_lin_spline), intent(in) :: &
      rho_lin_spline_vals, & ! rho values at the given altitudes
      rho_lin_spline_levels  ! altitudes for the given rho values
    ! Note: both these arrays need to be sorted from low to high altitude

    real( kind = core_rknd ), dimension(ngrdcol), intent(in) :: &
      p_sfc

    !--------------------- Output Variables ---------------------
    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(inout) :: &
      thlm_f, & ! Potential Temperature forcing     [K/s]
      rtm_f,  & ! Total Water Mixing Ration forcing [kg/kg/s]
      um_ref, & ! um reference                      [m/s]
      vm_ref, & ! vm reference                      [m/s]
      um_f,   & ! um tendency                       [m/s/s]
      vm_f,   & ! vm tendency                       [m/s/s]
      wm_zt,  & ! subsidence on zt grid             [m/s]
      ug,     & ! u geostrophic wind                [m/s]
      vg        ! v geostrophic wind                [m/s]

    real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(inout) :: &
      wm_zm     ! subsidense on zm grid             [m/s]

    real( kind = core_rknd ), dimension(ngrdcol,nzt, sclr_dim), intent(inout) :: &
      sclrm_forcing ! Scalar forcing [-]

    real( kind = core_rknd ), dimension(ngrdcol,nzt, edsclr_dim), intent(inout) :: &
      edsclrm_forcing ! Edscalar forcing [-]

    !--------------------- Local Variables ---------------------
    integer :: n ! Forcing field index [-]
    integer :: before_time, after_time ! Bracketing forcing time indices [-]
    integer :: iv ! Index of the field passed to remapping [-]

    integer :: i, k ! Column and level indices [-]

    real( kind = core_rknd ), dimension(ngrdcol,gr_dycore%nzt) :: &
      temp_array_dycore ! Time-interpolated forcing on the dycore grid [units vary]

    real( kind = core_rknd ), dimension(ngrdcol,nzt,nforcings) :: &
      forcings_array ! Time-interpolated fields (column, level, variable) [units vary]

    real( kind = core_rknd ) :: time_frac ! Fraction between bracketing input times [-]

    logical :: &
      l_zt_variable ! Remap a field on thermodynamic levels [-]

    !--------------------- Begin Code ---------------------
    time_frac = -one ! Default initialization
    iv = 1
    l_zt_variable = .true.

    call time_select( time, size(time_forcing_given), time_forcing_given, &
                                 before_time, after_time, time_frac )

    do n = 2, nforcings

      do k = 1, gr_dycore%nzt
        do i = 1, ngrdcol
          temp_array_dycore(i,k) = linear_interp_factor &
            ( time_frac, t_dependent_forcing_data(i,k,n,after_time), &
              t_dependent_forcing_data(i,k,n,before_time) )
        end do
      end do

      forcings_array(:,:,n) = remap_vals_to_target( ngrdcol, &
                                                    gr_dycore, gr, &
                                                    gr_dycore%nzt, &
                                                    temp_array_dycore, &
                                                    gr%nzt, &
                                                    total_idx_rho_lin_spline, &
                                                    rho_lin_spline_vals, &
                                                    rho_lin_spline_levels, &
                                                    iv, p_sfc, &
                                                    grid_remap_method, &
                                                    l_zt_variable )
    end do

    call apply_time_dependent_forcings_from_array( &
                                                   ngrdcol, nzm, nzt, &
                                                   sclr_dim, edsclr_dim, sclr_idx, &
                                                   l_do_expldiff_rtm_thlm, &
                                                   gr, rtm, rho, exner,  &
                                                   forcings_array, &
                                                   thlm_f, rtm_f, um_ref, vm_ref, um_f, vm_f, &
                                                   wm_zt, wm_zm,  ug, vg, &
                                                   sclrm_forcing, edsclrm_forcing )

    return

  end subroutine apply_time_dependent_forcings_from_dycore

  !================================================================================================
  subroutine time_select( time, nvar, time_array, &
                          before_time, after_time, time_frac )
    !
    ! Description: 
    !   This subroutine determines which indexes of the given
    !   time_array should be used when interpolating a value
    !   at the specified time and the location of time between
    !   these indexes.
    !
    ! References:
    !   None
    !---------------------------------------------------------------------------------

    use clubb_precision, only: time_precision, core_rknd ! Variable(s)

    use constants_clubb, only: fstderr ! Constant(s)

    implicit none

    ! External
    intrinsic :: real

    ! Input Variable(s)

    integer, intent(in) :: nvar                     ! Number of array elements [-]

    real(kind=time_precision), intent(in) :: time   ! Target time              [s]

    real( kind = core_rknd ), dimension(nvar), intent(in) :: time_array ! Array of times [s]

    ! Output Variable(s)

    integer, intent(out) :: &
      after_time, &  ! Index of a time later than the target time [-]
      before_time      ! Index of time before the target time       [-]

    real( kind = core_rknd ), intent(out) :: &
      time_frac      ! The fraction representing the point where time
                     ! is located between after_time and before_time [-]

    ! Local Variable(s)

    integer :: k

    ! ----------------- Begin Code --------------------

    ! Default initialization
    before_time = -999
    after_time = -999

    ! convert time to a real so it has the same precision as the values
    ! in time_array   
    if( real( time, kind = core_rknd ) < time_array(1) ) then
      
      ! If time is less than the lowest value in time_array, an invalid
      ! time has been provided. Stop execution.
      write(fstderr,*) "In subroutine time_select:"
      write(fstderr,*) "Selected time is before the first (begin) time"
      write(fstderr,*) "at which data are available.  Cannot interpolate."
      error stop

    else if ( real( time, kind = core_rknd ) > time_array(nvar) ) then

      ! If time is greater than the highest value in time_array, an invalid
      ! time has been provided. Stop execution.
      write(fstderr,*) "In subroutine time_select:"
      write(fstderr,*) "Selected time is after the last (end) time"
      write(fstderr,*) "at which data are available.  Cannot interpolate."
      error stop

    else

      do k=1,nvar-1

        if( time_array(k) >= time_array(k+1) ) then

            ! If times are not increasing then they aren't sorted properly.
            write(fstderr,*) "In subroutine time_select:"
            write(fstderr,*) "Times are not sorted. Check the configured surface and forcing files."
            error stop

        end if 

        if ( before_time == -999 .and. &
             (real( time, kind = core_rknd ) >= time_array(k)) .and. &
             (real( time, kind = core_rknd ) <= time_array(k+1)) ) then

          before_time = k
          after_time = k+1

        end if

      end do

    end if

    ! Compute the position of time between before_time and after_time
    ! as a fraction.
    time_frac = real( ( real( time, kind = core_rknd ) - time_array(before_time) ) / &
                ( time_array(after_time) - time_array(before_time) ), kind = core_rknd )

    return

  end subroutine time_select

!===========================================================================================
end module time_dependent_input
