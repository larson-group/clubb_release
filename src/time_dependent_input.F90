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

  use input_reader, only: &
    two_dim_read_var, &
    one_dim_read_var

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

  type(two_dim_read_var), private, dimension(nforcings) :: &
    t_dependent_forcing_data ! Data structure that defines the change in input
                             ! files over time
!$omp threadprivate( t_dependent_forcing_data )

  type(two_dim_read_var), private, dimension(nforcings) :: &
    t_dependent_forcing_data_f_grid ! Data structure that defines the change in input
                                    ! files over time but stores the raw file input.
                                    ! Is only set if we want to adapt grid but dont want to
                                    ! simulate forcings input from the host model (dycore grid)
!$omp threadprivate( t_dependent_forcing_data_f_grid )

  type(one_dim_read_var), private :: dimension_var ! Data structure that describes other 
                                                   ! dimension of the two_dim_read_var

!$omp threadprivate( dimension_var )

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

  ! File path constants
  character(len=*), private, parameter :: input_path = "../input/case_setups/"

  character(len=*), private, parameter :: forcings_path = "_forcings.in"

  character(len=*), private, parameter :: sfc_path = "_sfc.in"

  private

  contains

  !================================================================================================
  subroutine initialize_t_dependent_input( iunit, runtype, nz, grid, p_in_Pa, &
                                           l_add_dycore_grid, &
                                           grid_adapt_in_time_method )
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
    integer, intent(in) :: iunit ! File I/O

    character(len=*), intent(in) :: runtype ! Runtype

    integer, intent(in) :: nz ! Size of the model grid

    real( kind = core_rknd ), dimension(nz), intent(in) :: grid ! Model grid

    real( kind = core_rknd ), dimension(nz), intent(in) :: p_in_Pa ! Pressure[Pa]

    logical, intent(in) :: &
      l_add_dycore_grid     ! Flag to set remapping from dycore on or off

    integer, intent(in) :: &
      grid_adapt_in_time_method   ! flag that says whether the grid should be adapted over
                                  ! time and if so, how the grid density function should
                                  ! be formed

    ! ----------------- Begin Code --------------------

    if ( .not. l_ignore_forcings ) then
      call initialize_t_dependent_forcings &
                   ( iunit, input_path//trim(runtype)//forcings_path, nz, grid, p_in_Pa, &
                     l_add_dycore_grid, grid_adapt_in_time_method )
    end if

    if ( .not. l_sfc_already_initialized ) then
      call initialize_t_dependent_sfc &
                     ( iunit, input_path//trim(runtype)//sfc_path )
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
  subroutine initialize_t_dependent_sfc( iunit, input_file )
    !
    !  Description: This subroutine reads in a file that details time dependent
    !  input values that vary in one dimension.
    !-----------------------------------------------------------------------------

    use input_reader, only: &
      read_one_dim_file, one_dim_read_var, & ! Procedure(s)
      fill_blanks_one_dim_vars, read_x_profile, &
      get_target_index, deallocate_one_dim_vars, &
      count_columns

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

    character(len=*), intent(in) :: input_file ! Path to surface.in file

    ! Local Variable(s)

    type(one_dim_read_var), allocatable, dimension(:) :: &
      retVars ! retVars stores the name of a variable (e.g. pressure),
              ! the name of the dimension the variable varies along (e.g. time), and
              ! the time-dependent values of the variable to be input into CLUBB.

    integer ::  &
      dim_size, & ! Number of time-dependent values of a variable to be input into CLUBB 
      nforcings       ! Number of variables with time-dependent input data


    ! ----------------- Begin Code --------------------

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

    l_sfc_already_initialized = .true.

    return 
  end subroutine initialize_t_dependent_sfc

  !================================================================================================
  subroutine initialize_t_dependent_forcings( iunit, input_file, nz, grid, p_in_Pa, &
                                              l_add_dycore_grid, &
                                              grid_adapt_in_time_method )
    !
    !  Description: This subroutine reads in a file that details time dependent
    !  input values that vary in two dimensions.
    !
    !-------------------------------------------------------------------------------------

    use input_reader, only: read_two_dim_file, two_dim_read_var, fill_blanks_two_dim_vars

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
    integer, intent(in) :: iunit ! File I/O

    character(len=*), intent(in) :: input_file ! Path to the input file

    integer, intent(in) :: nz  ! Size of Model Grid [-]

    real( kind = core_rknd ), intent(in), dimension(nz) :: grid ! Altitudes of Grid [m]

    real( kind = core_rknd ), intent(in), dimension(nz) :: p_in_Pa ! Pressure [Pa]

    logical, intent(in) :: &
      l_add_dycore_grid     ! Flag to set remapping from dycore on or off

    integer, intent(in) :: &
      grid_adapt_in_time_method   ! flag that says whether the grid should be adapted over
                                  ! time and if so, how the grid density function should
                                  ! be formed

    ! Local Variables

    integer :: i, n_f_grid_z, n_f_grid_t

    ! ----------------- Begin Code --------------------

    if ( .not. allocated(t_dependent_forcing_data_f_grid(1)%values) ) then
      ! Only read in the forcing data from the input file if they aren't already
      ! stored in the global variable
      call read_two_dim_file( iunit, nforcings, input_file, &
                              t_dependent_forcing_data_f_grid, dimension_var )
    end if

    n_f_grid_z = size( t_dependent_forcing_data_f_grid(1)%values, 1 )

    n_f_grid_t = size( dimension_var%values )

    ! Fill in blanks with linear interpolation. Whole profiles of -999.9 will
    ! remain that way thus marking them blank.
    call fill_blanks_two_dim_vars( nforcings, dimension_var, t_dependent_forcing_data_f_grid )

    do i=1, nforcings
      ! Only allocate t_dependent_forcing_data, if it isn't already stored in the global variable
      if ( .not. allocated(t_dependent_forcing_data(i)%values) ) then
        allocate( t_dependent_forcing_data(i)%values(1:nz,1:n_f_grid_t) )
      end if
    end do

    select case( t_dependent_forcing_data_f_grid(1)%name )
    case( z_name )

      t_dependent_forcing_data(1)%name = z_name
      t_dependent_forcing_data(1)%values(:,1:n_f_grid_t) = spread(grid,2,n_f_grid_t)

    case( pressure_name )

      t_dependent_forcing_data(1)%name = pressure_name
      t_dependent_forcing_data(1)%values(:,1:n_f_grid_t) = spread(-p_in_Pa, 2, n_f_grid_t )
      t_dependent_forcing_data_f_grid(1)%values = -t_dependent_forcing_data_f_grid(1)%values

    case default
      error stop "Incompatible grid type in first element of t_dependent_forcings."
    end select

    ! Interpolate the time dependent input data to the appropriate grid.

    do i=2, nforcings

      t_dependent_forcing_data(i)%name = t_dependent_forcing_data_f_grid(i)%name
      t_dependent_forcing_data(i)%values = read_to_grid( nforcings, n_f_grid_z, n_f_grid_t, &
                         nz,t_dependent_forcing_data(1)%values(:,1), &
                         t_dependent_forcing_data_f_grid, t_dependent_forcing_data(i)%name )

    end do

    if ( l_add_dycore_grid .or. grid_adapt_in_time_method == no_grid_adaptation ) then
      ! the array should only be kept in memory, if we want to use grid adaptation
      ! (grid_adapt_in_time_method > no_grid_adaptation) and no simulating forcings input from the
      ! host model on the dycore grid (l_add_dycore_grid == .true.), since then we
      ! need the raw input from the file during the CLUBB run if the grid changes
      do i = 1, nforcings
        if ( allocated( t_dependent_forcing_data_f_grid(i)%values ) ) then
          deallocate( t_dependent_forcing_data_f_grid(i)%values )
        end if
      end do
    end if

    return
  end subroutine initialize_t_dependent_forcings


  !================================================================================================
  subroutine finalize_t_dependent_forcings()
    !
    !   Description: Clears memory initialized in initialize_t_dependent_forcings.
    !   This should be called at the end of the model
    !----------------------------------------------------------

    implicit none

    ! Local Variable
    integer :: i

    ! ----------------- Begin Code --------------------

    ! for the case that we use grid adaptation and do not want to simulate the forcings
    ! input from the host model on the dycore grid, we keep the raw input also in memory
    do i = 1, nforcings
      if ( allocated( t_dependent_forcing_data_f_grid(i)%values ) ) then
        deallocate( t_dependent_forcing_data_f_grid(i)%values )
      end if
    end do

    do i=1, nforcings
      deallocate( t_dependent_forcing_data(i)%values )
    end do

    deallocate( dimension_var%values ) 

    return
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
  function read_to_grid( ntwo_dim_vars, dim_size, other_dim_size, &
                         nz, grid, two_dim_vars, target_name ) result(var)
    !
    !  Description: This is a helper function for doing the translation from the
    !  forcing grid to the model grid.
    !
    !------------------------------------------------------------------------------------

    use input_reader, only: read_x_table, two_dim_read_var

    use interpolation, only: zlinterp_fnc

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    integer, intent(in) :: &
      ntwo_dim_vars, &
      dim_size, &
      other_dim_size, &
      nz

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      grid

    type(two_dim_read_var), dimension(ntwo_dim_vars), intent(in) :: &
      two_dim_vars

    character(len=*), intent(in) :: &
      target_name

    real( kind = core_rknd ), dimension(dim_size, other_dim_size) :: temp_var

    real( kind = core_rknd ), dimension(nz, other_dim_size) :: var

    integer i

    ! ----------------- Begin Code --------------------

    temp_var = read_x_table( ntwo_dim_vars,  dim_size, other_dim_size, target_name, two_dim_vars )

    do i=1, other_dim_size
      var(:,i) = zlinterp_fnc( nz, dim_size, grid, &
                                    two_dim_vars(1)%values(:,i), temp_var(:,i) )
    end do

    return

  end function read_to_grid

  !================================================================================================
  subroutine apply_time_dependent_forcings_from_array( &
              ngrdcol, nzm, nzt, &
              sclr_dim, edsclr_dim, sclr_idx, &
              gr, rtm, rho, exner,  &
              forcings_array, &
              thlm_f, rtm_f, um_ref, vm_ref, um_f, vm_f, &
              wm_zt, wm_zm,  ug, vg, &
              sclrm_forcing, edsclrm_forcing )
    !
    !  Description: This subroutine is a helper subroutine, that takes the forcings
    !               in a 2D matrix and initializes the forcing variables.
    !
    !---------------------------------------------------------------------------------

    use constants_clubb, only: &
      grav, & ! Variable(s)
      sec_per_hr, &
      pascal_per_mb, &
      fstderr

    use interpolation, only: &
      linear_interp_factor ! Procedure(s)

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
      grid, & ! Type
      zt2zm   ! Procedure(s)

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

    type (grid), intent(in) :: &
      gr

    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
      exner,   & ! Exner Function                             [-]
      rho,     & ! Air Density                                [kg/m^3]
      rtm        ! Total Water Mixing Ratio                   [kg/kg]

    real( kind = core_rknd ), dimension(nforcings,nzt), intent(in) :: forcings_array

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
    integer :: i, k, n

    real( kind = core_rknd ), dimension(nzt) :: temp_array

    !--------------------- Begin Code ---------------------

    !$acc data create( temp_array )

    ! Parse the values in t_dependent_forcing_data for CLUBB compatible forcing
    ! data.
    do n=2, nforcings

      temp_array = forcings_array(n,:)

      ! Check to see if temp_array is an actual profile or a dummy profile
      ! If it is a dummy profile we dont want it to apply itself as it may
      ! overwrite legitimate information from another source.
      if( .not. any( abs(temp_array - (-999.9_core_rknd)) < &
          abs(temp_array + (-999.9_core_rknd)) / 2 * eps ) ) then

        !$acc update device( temp_array )
          
        select case (t_dependent_forcing_data(n)%name)
        case(temperature_f_name, theta_f_name, thetal_f_name)

          select case(t_dependent_forcing_data(n)%name)
          case(temperature_f_name)

            !$acc parallel loop gang vector collapse(2) default(present)
            do k = 1, nzt
              do i = 1, ngrdcol
                thlm_f(i,k) = temp_array(k) / exner(i,k)
              end do
            end do

          case(theta_f_name)

            !$acc parallel loop gang vector collapse(2) default(present)
            do k = 1, nzt
              do i = 1, ngrdcol
                thlm_f(i,k) = temp_array(k) ! n am not sure on the conversion of this
              end do
            end do

          case(thetal_f_name)

            !$acc parallel loop gang vector collapse(2) default(present)
            do k = 1, nzt
              do i = 1, ngrdcol
                thlm_f(i,k) = temp_array(k)
              end do
            end do

          end select

          if ( sclr_idx%iisclr_thl > 0 ) then

            !$acc parallel loop gang vector collapse(2) default(present)
            do k = 1, nzt
              do i = 1, ngrdcol
                sclrm_forcing(i,k,sclr_idx%iisclr_thl) = thlm_f(i,k)
              end do
            end do

          end if

          if ( sclr_idx%iiedsclr_thl > 0 ) then

            !$acc parallel loop gang vector collapse(2) default(present)
            do k = 1, nzt
              do i = 1, ngrdcol
                edsclrm_forcing(i,k,sclr_idx%iiedsclr_thl) = thlm_f(i,k)
              end do
            end do

          end if

        case(rt_f_name, sp_humidity_f_name)

          select case(t_dependent_forcing_data(n)%name)
          case(sp_humidity_f_name)

            !$acc parallel loop gang vector collapse(2) default(present)
            do k = 1, nzt
              do i = 1, ngrdcol
                rtm_f(i,k) = temp_array(k) * ( 1._core_rknd + rtm(i,k) )**2
              end do
            end do

          case(rt_f_name )

            !$acc parallel loop gang vector collapse(2) default(present)
            do k = 1, nzt
              do i = 1, ngrdcol
                rtm_f(i,k) = temp_array(k)
              end do
            end do

          end select

          if ( sclr_idx%iisclr_rt  > 0 ) then

            !$acc parallel loop gang vector collapse(2) default(present)
            do k = 1, nzt
              do i = 1, ngrdcol
                sclrm_forcing(i,k,sclr_idx%iisclr_rt)  = rtm_f(i,k)
              end do
            end do

          end if

          if ( sclr_idx%iiedsclr_rt  > 0 ) then

            !$acc parallel loop gang vector collapse(2) default(present)
            do k = 1, nzt
              do i = 1, ngrdcol
                edsclrm_forcing(i,k,sclr_idx%iiedsclr_rt)  = rtm_f(i,k)
              end do
            end do

          end if

        case(um_ref_name)

          !$acc parallel loop gang vector collapse(2) default(present)
            do k = 1, nzt
            do i = 1, ngrdcol
              um_ref(i,k) = temp_array(k)
            end do
          end do

        case(vm_ref_name)

          !$acc parallel loop gang vector collapse(2) default(present)
          do k = 1, nzt
            do i = 1, ngrdcol
              vm_ref(i,k) = temp_array(k)
            end do
          end do

        case(um_f_name)

          !$acc parallel loop gang vector collapse(2) default(present)
          do k = 1, nzt
            do i = 1, ngrdcol
              um_f(i,k) = temp_array(k)
            end do
          end do

        case(vm_f_name)

          !$acc parallel loop gang vector collapse(2) default(present)
          do k = 1, nzt
            do i = 1, ngrdcol
              vm_f(i,k) = temp_array(k)
            end do
          end do

        case(wm_name, omega_name, omega_mb_hr_name)

          select case(t_dependent_forcing_data(n)%name)
          case(wm_name)

            !$acc parallel loop gang vector collapse(2) default(present)
            do k = 1, nzt
              do i = 1, ngrdcol
                wm_zt(i,k) = temp_array(k)
              end do
            end do

          case(omega_name)

            !$acc parallel loop gang vector collapse(2) default(present)
            do k = 1, nzt
              do i = 1, ngrdcol
                wm_zt(i,k) = - temp_array(k) / (grav * rho(i,k))
              end do
            end do

          case(omega_mb_hr_name)

            !$acc parallel loop gang vector default(present)
            do k = 1, nzt
              temp_array(k) = temp_array(k) * pascal_per_mb / sec_per_hr
            end do

            !$acc parallel loop gang vector collapse(2) default(present)
            do k = 1, nzt
              do i = 1, ngrdcol

                wm_zt(i,k) = - temp_array(k) / (grav * rho(i,k))

              end do
            end do

          end select

          wm_zm = zt2zm( nzm, nzt, ngrdcol, gr, wm_zt )

        case(ug_name)

          !$acc parallel loop gang vector collapse(2) default(present)
          do k = 1, nzt
            do i = 1, ngrdcol
              ug(i,k) = temp_array(k)
            end do
          end do

        case(vg_name)

          !$acc parallel loop gang vector collapse(2) default(present)
          do k = 1, nzt
            do i = 1, ngrdcol
              vg(i,k) = temp_array(k)
            end do
          end do

        case default

          write(fstderr, *) "Incompatable forcing type: "//t_dependent_forcing_data(n)%name
          error stop

        end select

      end if 

    end do ! 1 .. nforcings

    !$acc end data

    return

  end subroutine apply_time_dependent_forcings_from_array

  !================================================================================================
  subroutine apply_time_dependent_forcings( &
              ngrdcol, nzm, nzt, &
              sclr_dim, edsclr_dim, sclr_idx, &
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
      grid, & ! Type
      zt2zm   ! Procedure(s)

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
    integer :: n, before_time, after_time

    real( kind = core_rknd ), dimension(nforcings,nzt) :: forcings_array

    real( kind = core_rknd ) :: time_frac

    !--------------------- Begin Code ---------------------

    time_frac = -one ! Default initialization

    call time_select( time, size(dimension_var%values), dimension_var%values, &
                                 before_time, after_time, time_frac )

    do n = 2, nforcings
      forcings_array(n,:) = linear_interp_factor &
                   ( time_frac, t_dependent_forcing_data(n)%values(:,after_time), &
                     t_dependent_forcing_data(n)%values(:,before_time) )
    end do

    call apply_time_dependent_forcings_from_array( &
                                                   ngrdcol, nzm, nzt, &
                                                   sclr_dim, edsclr_dim, sclr_idx, &
                                                   gr, rtm, rho, exner,  &
                                                   forcings_array, &
                                                   thlm_f, rtm_f, um_ref, vm_ref, um_f, vm_f, &
                                                   wm_zt, wm_zm,  ug, vg, &
                                                   sclrm_forcing, edsclrm_forcing )

    return

  end subroutine apply_time_dependent_forcings

  !================================================================================================
  subroutine apply_time_dependent_forcings_from_dycore( &
              ngrdcol, nzm, nzt, &
              sclr_dim, edsclr_dim, sclr_idx, &
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
      grid, & ! Type
      zt2zm   ! Procedure(s)

    use clubb_precision, only: &
      time_precision, &
      core_rknd ! Variable(s)

    use array_index, only: &
      sclr_idx_type

    use grid_adaptation_module, only: &
      remapping_matrix_zt_values, &
      remap_vals_to_target

    use error_code, only: &
      clubb_at_least_debug_level

    use constants_clubb, only: &
      fstderr

    use model_flags, only: &
      cons_ullrich_remap

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
    integer :: n, before_time, after_time

    real( kind = core_rknd ), dimension(gr_dycore%nzt) :: temp_array_dycore

    real( kind = core_rknd ), dimension(1,nforcings,nzt) :: forcings_array

    real( kind = core_rknd ) :: time_frac

    ! right now only one forcings profile can be applied, so ngrdcol=1
    real( kind = core_rknd ), dimension(1,gr%nzt,gr_dycore%nzt) :: R_ij

    !--------------------- Begin Code ---------------------
    time_frac = -one ! Default initialization

    call time_select( time, size(dimension_var%values), dimension_var%values, &
                                 before_time, after_time, time_frac )

    if ( grid_remap_method == cons_ullrich_remap ) then

      ! right now only one forcings profile can be applied, so ngrdcol=1
      call remapping_matrix_zt_values( 1, &                          ! Intent(in)
                                       gr_dycore, gr, &              ! Intent(in)
                                       total_idx_rho_lin_spline, &   ! Intent(in)
                                       rho_lin_spline_vals(1,:), &   ! Intent(in)
                                       rho_lin_spline_levels(1,:), & ! Intent(in)
                                       R_ij )                        ! Intent(out)
    end if

    do n = 2, nforcings

      temp_array_dycore = linear_interp_factor &
                          ( time_frac, t_dependent_forcing_data(n)%values(:,after_time), &
                            t_dependent_forcing_data(n)%values(:,before_time) )

      if ( grid_remap_method == cons_ullrich_remap ) then

        forcings_array(:,n,:) = remap_vals_to_target( 1, &
                                                      gr_dycore%nzm, &
                                                      gr%nzm, &
                                                      gr_dycore%zm(1,:), &
                                                      gr%zm(1,:), &
                                                      total_idx_rho_lin_spline, &
                                                      rho_lin_spline_vals(1,:), &
                                                      rho_lin_spline_levels(1,:), &
                                                      temp_array_dycore, &
                                                      R_ij(1,:,:), p_sfc(1) )

      else
        write(fstderr,*) 'There is currently no method implemented for grid_remap_method=', &
                         grid_remap_method, '. Set flag to different value.'
        error stop 'Invalid value for flag grid_remap_method.'
      end if
    end do

    call apply_time_dependent_forcings_from_array( &
                                                   ngrdcol, nzm, nzt, &
                                                   sclr_dim, edsclr_dim, sclr_idx, &
                                                   gr, rtm, rho, exner,  &
                                                   forcings_array(1,:,:), &
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
            write(fstderr,*) "times are not sorted. Check (case)_sfc.in "
            write(fstderr,*) "and (case)_forcings.in, located in input/case_setups."
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
