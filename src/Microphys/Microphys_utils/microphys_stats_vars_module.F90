!-------------------------------------------------------------------------------
! $Id$
!===============================================================================
module microphys_stats_vars_module

! Description:
!   This module contains the derived type microphys_stats_vars_type, which is used to
!   feed output variables out of microphysics schemes

  use clubb_precision, only: &
    core_rknd

  implicit none

  private ! Set Default Scope

  integer, parameter :: &
    stats_name_len = 64

  public :: microphys_stats_vars_type, microphys_stats_alloc, microphys_put_var, &
            microphys_get_var_by_name, &
            microphys_stats_accumulate, microphys_stats_cleanup

  ! (no debug flags)

  type microphys_stats_vars_type

    logical :: &
      l_allocated = .false. ! This is set to true when the structure is allocated

    integer :: &
      nz, &       ! Number of vertical levels
      num_vars, & ! Number of output variables from microphysics
      alloc_size  ! Size of allocated var_names and output_values arrays

    ! Names of the output variables (used for stats mapping)
    ! An array of names corresponding to the output variables.
    character( len = stats_name_len ), dimension(:), allocatable :: &
      var_names

    ! Values of the output variables (over the vertical domain)
    real( kind = core_rknd ), dimension(:,:), allocatable :: &
      output_values

  end type microphys_stats_vars_type

  contains

  !-----------------------------------------------------------------------
  subroutine microphys_stats_alloc( nz, num_vars, microphys_stats_vars )

  ! Description:
  !   Allocates a new microphysics stats type

  ! References:
  !   None
  !-----------------------------------------------------------------------

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nz, &                  ! Number of vertical levels for statistics output
      num_vars               ! Number of statistics variables desired
                             ! (actual number output can be less)

    ! Input/Output Variables
    type(microphys_stats_vars_type), intent(out) :: &
      microphys_stats_vars   ! Unallocated microphys_stats_vars_type object

  !-----------------------------------------------------------------------

    !----- Begin Code -----

    if ( num_vars <= 0 .or. nz <= 0 ) then
      microphys_stats_vars%nz = max( nz, 0 )
      microphys_stats_vars%alloc_size = 0
      microphys_stats_vars%num_vars = 0
      microphys_stats_vars%l_allocated = .false.
      return
    end if

    allocate( microphys_stats_vars%var_names(num_vars), &
              microphys_stats_vars%output_values(nz,num_vars) )

    microphys_stats_vars%nz = nz
    microphys_stats_vars%alloc_size = num_vars
    microphys_stats_vars%num_vars = 0 ! Since no variables have been put into structure
    microphys_stats_vars%var_names = ""

    microphys_stats_vars%l_allocated = .true.

    return
  end subroutine microphys_stats_alloc
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  subroutine microphys_put_var( var_name, value, microphys_stats_vars )

  ! Description:
  !   Places a variable in the microphys_stats_vars structure

  ! References:
  !   None
  !-----------------------------------------------------------------------

    implicit none

    ! Input/Output Variables
    type(microphys_stats_vars_type), intent(inout) :: &
      microphys_stats_vars   ! The structure to place the variable into

    ! Input Variables
    character(len=*), intent(in) :: &
      var_name               ! Name of the variable

    real( kind = core_rknd ), dimension(microphys_stats_vars%nz), intent(in) :: &
      value

  !-----------------------------------------------------------------------

    !----- Begin Code -----
    if ( .not. microphys_stats_vars%l_allocated ) return

    if ( microphys_stats_vars%num_vars >= microphys_stats_vars%alloc_size ) then
      ! There is no more room in the structure.
      error stop "Fatal error in microphys_put_var: allocation exhausted"
    end if

    microphys_stats_vars%num_vars = microphys_stats_vars%num_vars + 1
    microphys_stats_vars%var_names(microphys_stats_vars%num_vars) = trim(var_name)
    microphys_stats_vars%output_values(:,microphys_stats_vars%num_vars) = value

    return
  end subroutine microphys_put_var
  !-----------------------------------------------------------------------

  subroutine microphys_stats_accumulate( microphys_stats_vars, &
                                         stats, icol )

  ! Description:
  !   Samples all variables stored in the microphys_stats_vars structure by
  !   using stat_update_var. The variables are sampled on the grid given as an
  !   argument to the subroutine. Therefore, there are no global variables used!

  ! References:
  !   None
  !-----------------------------------------------------------------------

    use stats_netcdf, only: &
      stats_type, &
      stats_update

    implicit none

    ! Input Variables
    type(microphys_stats_vars_type), intent(in) :: &
      microphys_stats_vars

    type(stats_type), intent(inout) :: &
      stats

    integer, intent(in) :: &
      icol

    ! Local Variablesbit
    integer :: i ! Loop variable

  !-----------------------------------------------------------------------

    !----- Begin Code -----
    if ( stats%l_sample ) then
      do i = 1, microphys_stats_vars%num_vars
        if ( microphys_stats_vars%nz == 1 ) then
          call stats_update( trim(microphys_stats_vars%var_names(i)), &
                            microphys_stats_vars%output_values(1,i), stats, icol)
        else
          call stats_update( trim(microphys_stats_vars%var_names(i)), &
                            microphys_stats_vars%output_values(:,i), stats, icol)
        end if
      end do
    end if

  end subroutine microphys_stats_accumulate
  !-----------------------------------------------------------------------

  function microphys_get_var_by_name( var_name, microphys_stats_vars ) result( stats_var )

  ! Description:
  !   Gets the specified statistics variable from the input structure by name.
  !   Get an index to the arrays in the microphys_stats_vars structure
  !   associated with the given stats name.

  ! References:
  !   None
  !-----------------------------------------------------------------------
    implicit none

    ! Input Variables
    character(len=*), intent(in) :: &
      var_name                ! Name of the variable

    type(microphys_stats_vars_type), intent(in) :: &
      microphys_stats_vars    ! The statistics structure

    ! Output Variable
    real( kind = core_rknd ), dimension(microphys_stats_vars%nz) :: &
      stats_var               ! The output statistics variable

    ! Local Variables
    integer :: ivar
    logical :: l_found

  !-----------------------------------------------------------------------

    ! Initialize variables
    stats_var = 0._core_rknd
    l_found = .false.

    !----- Begin Code -----
    do ivar = 1, microphys_stats_vars%num_vars
      if ( trim(microphys_stats_vars%var_names(ivar)) == trim(var_name) ) then
        stats_var = microphys_stats_vars%output_values(:,ivar)
        l_found = .true.
        exit
      end if
    end do

    if ( .not. l_found ) then
      error stop "Variable not found in microphys_get_var_by_name"
    end if

    return
  end function microphys_get_var_by_name
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  subroutine microphys_stats_cleanup( microphys_stats_vars )

  ! Description:
  !   Deallocates all (dynamic) memory associated with the
  !   microphys_stats_vars object

  ! References:
  !   None
  !-----------------------------------------------------------------------

    implicit none

    ! Input/Output Variables
    type(microphys_stats_vars_type), intent(inout) :: &
      microphys_stats_vars   ! Unallocated microphys_stats_vars_type object

  !-----------------------------------------------------------------------

    !----- Begin Code -----

    if ( .not. microphys_stats_vars%l_allocated ) return

    if ( allocated(microphys_stats_vars%var_names) ) then
      deallocate( microphys_stats_vars%var_names )
    end if
    if ( allocated(microphys_stats_vars%output_values) ) then
      deallocate( microphys_stats_vars%output_values )
    end if

    microphys_stats_vars%l_allocated = .false.
    return
  end subroutine microphys_stats_cleanup
  !-----------------------------------------------------------------------

end module microphys_stats_vars_module
