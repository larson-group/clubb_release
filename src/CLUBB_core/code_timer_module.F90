!-------------------------------------------------------------------------------
! code_timer_module.F90
!
! Centralized named timing calls for CLUBB.
!
! Usage:
!
!   Call once before timed regions:
!     - timer_initialize(): initialize the active timing backend
!
!
!   Call around each timed region:
!     - timer_start(): begin a named timer region
!
!     - timer_stop(): end the most recent matching named timer region
!
!
!   Call once after timed regions:
!     - timer_finalize(): write timer output and finalize the active backend
!
!
! Understanding:
!
!   Backends:
!     - When CLUBB is compiled with GPTL, timer_initialize(), timer_start(),
!       timer_stop(), and timer_finalize() forward to the matching GPTL calls.
!
!     - When CLUBB is compiled without GPTL, this module uses a small cpu_time
!       backend that records a flat table of named timers.
!
!   Timer names:
!     - Timer regions are identified by the character string passed to
!       timer_start() and timer_stop().  The cpu_time backend creates timer
!       registry entries dynamically when new names are first started.
!
!     - Timer names must not be empty.  The same name string must be passed to
!       timer_stop() that was passed to the matching timer_start().
!
!   Nesting:
!     - Timer calls must be explicitly nested.  For example, if timer B starts
!       inside timer A, timer B must stop before timer A stops.
!
!     - The cpu_time backend enforces this stack order and reports an error if
!       timer_stop() is called for the wrong active timer.
!
!   cpu_time output:
!     - The cpu_time backend writes a summary table to the file passed to
!       timer_finalize().  The table includes call count, inclusive time,
!       exclusive time, average inclusive time per call, and exclusive percent.
!
!     - Inclusive time is the full elapsed time in a timer region.  Exclusive
!       time subtracts time spent in nested timer regions.
!
!   GPTL output:
!     - The GPTL backend writes output through GPTLpr_file() using the file
!       passed to timer_finalize().  GPTL controls the exact output format.
!
!-------------------------------------------------------------------------------
module code_timer_module

  use iso_fortran_env, only: &
    output_unit

#ifndef GPTL
  use clubb_precision, only: &
    core_rknd
#endif

  use constants_clubb, only: &
    fstderr

#ifdef GPTL
  use gptl, only: &
    GPTLsetoption, &
    GPTLprint_method, &
    GPTLfull_tree, &
    GPTLoverhead, &
    GPTLabort_on_error, &
    GPTLinitialize, &
    GPTLstart, &
    GPTLstop, &
    GPTLpr_file, &
    GPTLfinalize
#endif

  implicit none

  private ! Set default scope

  logical :: &
    l_timer_initialized = .false., &
    l_timer_enabled = .true.

  !$omp threadprivate( l_timer_initialized, l_timer_enabled )

#ifndef GPTL
  integer, parameter :: &
    registry_chunk_size = 16, & ! Number of timer entries to add when the registry grows.
    stack_chunk_size    = 32, & ! Number of stack entries to add when the stack grows.
    display_name_length = 48    ! Width of the timer name in the cpu_time printout.

  type timer_entry_t
    character(len=:), allocatable :: name
    real( kind = core_rknd ) :: inclusive_time = 0.0_core_rknd
    real( kind = core_rknd ) :: exclusive_time = 0.0_core_rknd
    integer :: call_count = 0
  end type timer_entry_t

  type timer_stack_entry_t
    integer :: timer_id = 0
    real( kind = core_rknd ) :: start_time = 0.0_core_rknd
    real( kind = core_rknd ) :: child_time = 0.0_core_rknd
  end type timer_stack_entry_t

  type(timer_entry_t), allocatable :: timer_entries(:)
  type(timer_stack_entry_t), allocatable :: timer_stack(:)

  integer :: &
    num_timers = 0, &
    stack_depth = 0

  !$omp threadprivate( timer_entries, timer_stack )
  !$omp threadprivate( num_timers, stack_depth )
#endif

  public :: &
    timer_initialize, &
    timer_disable, &
    timer_start, &
    timer_stop, &
    timer_finalize

  contains

  !-----------------------------------------------------------------------
  subroutine timer_initialize()

  ! Description:
  !   Initializes the active timing backend.
  !-----------------------------------------------------------------------

    implicit none

#ifdef GPTL
    !--------------------------- Local Variables ---------------------------
    integer :: &
      ret_code    ! Return code from GPTL calls. [-]
#endif

    !--------------------------- Begin Code ---------------------------

    l_timer_enabled = .true.

#ifdef GPTL
    if ( l_timer_initialized ) then
      ret_code = GPTLfinalize()
    end if

    ret_code = GPTLsetoption( GPTLprint_method, GPTLfull_tree )
    ret_code = GPTLsetoption( GPTLabort_on_error, 1 )
    ret_code = GPTLsetoption( GPTLoverhead, 0 )
    ret_code = GPTLinitialize()
    if ( ret_code /= 0 ) call timer_fail( "GPTLinitialize failed." )
#else
    if ( l_timer_initialized ) call require_empty_stack( "timer_initialize" )

    call reset_cpu_timer_state()
#endif

    l_timer_initialized = .true.

    return
  end subroutine timer_initialize
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  subroutine timer_disable()

  ! Description:
  !   Disables the active timing backend and discards any recorded timer state.
  !-----------------------------------------------------------------------

    implicit none

#ifdef GPTL
    !--------------------------- Local Variables ---------------------------
    integer :: &
      ret_code    ! Return code from GPTL calls. [-]
#endif

    !--------------------------- Begin Code ---------------------------

    if ( l_timer_initialized ) then
#ifdef GPTL
      ret_code = GPTLfinalize()
#else
      call reset_cpu_timer_state()
#endif
    end if

    l_timer_initialized = .false.
    l_timer_enabled = .false.

    return
  end subroutine timer_disable
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  subroutine timer_start( timer_name )

  ! Description:
  !   Starts or enters a named timer region.
  !-----------------------------------------------------------------------

    implicit none

    !--------------------------- Input Variables ---------------------------
    character(len=*), intent(in) :: &
      timer_name    ! Name of the timer region. [-]

#ifdef GPTL
    !--------------------------- Local Variables ---------------------------
    integer :: &
      ret_code    ! Return code from GPTL calls. [-]
#endif

    !--------------------------- Begin Code ---------------------------

    if ( .not. l_timer_enabled ) return

    if ( .not. l_timer_initialized ) call timer_initialize()

#ifdef GPTL
    ret_code = GPTLstart( trim( timer_name ) )
#else
    call cpu_timer_start( timer_name )
#endif

    return
  end subroutine timer_start
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  subroutine timer_stop( timer_name )

  ! Description:
  !   Stops or exits a named timer region.
  !-----------------------------------------------------------------------

    implicit none

    !--------------------------- Input Variables ---------------------------
    character(len=*), intent(in) :: &
      timer_name    ! Name of the timer region. [-]

#ifdef GPTL
    !--------------------------- Local Variables ---------------------------
    integer :: &
      ret_code    ! Return code from GPTL calls. [-]
#endif

    !--------------------------- Begin Code ---------------------------

    if ( .not. l_timer_enabled ) return

    if ( .not. l_timer_initialized ) then
      call timer_fail( "timer_stop called before timer_initialize or timer_start." )
    end if

#ifdef GPTL
    ret_code = GPTLstop( trim( timer_name ) )
#else
    call cpu_timer_stop( timer_name )
#endif

    return
  end subroutine timer_stop
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  subroutine timer_finalize( output_file )

  ! Description:
  !   Prints timer output and finalizes the active timing backend.
  !-----------------------------------------------------------------------

    implicit none

    !--------------------------- Input Variables ---------------------------
    character(len=*), intent(in) :: &
      output_file    ! File used for timer output. [-]

#ifdef GPTL
    !--------------------------- Local Variables ---------------------------
    integer :: &
      ret_code    ! Return code from GPTL calls. [-]
#else
    !--------------------------- Local Variables ---------------------------
    integer :: &
      file_unit, & ! Unit used for the cpu_time summary file. [-]
      io_status   ! Status from opening the cpu_time summary file. [-]
#endif

    character(len=1024) :: &
      output_path    ! Absolute path used for timer output. [-]

    !--------------------------- Begin Code ---------------------------

    if ( .not. l_timer_enabled ) return

    if ( .not. l_timer_initialized ) return

    if ( len_trim( output_file ) < 1 ) then
      call timer_fail( "Timer output file must not be empty." )
    end if

    output_path = timer_absolute_path( output_file )
    write( output_unit, '(a)' ) "Writing timing results to: " // trim( output_path )

#ifdef GPTL
    ret_code = GPTLpr_file( trim( output_path ) )
    ret_code = GPTLfinalize()
#else
    call require_empty_stack( "timer_finalize" )

    open( newunit = file_unit, file = trim( output_path ), status = "replace", &
          action = "write", iostat = io_status )
    if ( io_status /= 0 ) then
      call timer_fail( "Unable to open timer output file: " // trim( output_path ) )
    end if
    call print_cpu_timer_summary( file_unit )
    close( file_unit )

    call reset_cpu_timer_state()
#endif

    l_timer_initialized = .false.

    return
  end subroutine timer_finalize
  !-----------------------------------------------------------------------

#ifndef GPTL
  !-----------------------------------------------------------------------
  subroutine cpu_timer_start( timer_name )

  ! Description:
  !   Starts or enters a named timer region using cpu_time.
  !-----------------------------------------------------------------------

    implicit none

    !--------------------------- Input Variables ---------------------------
    character(len=*), intent(in) :: &
      timer_name    ! Name of the timer region. [-]

    !--------------------------- Local Variables ---------------------------
    integer :: &
      timer_id    ! Registry index for timer_name. [-]

    real( kind = core_rknd ) :: &
      start_time    ! cpu_time value at timer start. [s]

    !--------------------------- Begin Code ---------------------------

    call cpu_time( start_time )

    timer_id = get_timer_id( timer_name, l_create = .true. )

    call ensure_stack_capacity( stack_depth + 1 )

    stack_depth = stack_depth + 1
    timer_stack(stack_depth)%timer_id = timer_id

    timer_entries(timer_id)%call_count = timer_entries(timer_id)%call_count + 1

    timer_stack(stack_depth)%start_time = start_time
    timer_stack(stack_depth)%child_time = 0.0_core_rknd

    return
  end subroutine cpu_timer_start
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  subroutine cpu_timer_stop( timer_name )

  ! Description:
  !   Stops or exits the current named cpu_time timer region.  Timer stops
  !   must match the most recent active timer start.
  !-----------------------------------------------------------------------

    implicit none

    !--------------------------- Input Variables ---------------------------
    character(len=*), intent(in) :: &
      timer_name    ! Name of the timer region. [-]

    !--------------------------- Local Variables ---------------------------
    integer :: &
      timer_id    ! Registry index for timer_name. [-]

    real( kind = core_rknd ) :: &
      stop_time, & ! cpu_time value at timer stop. [s]
      elapsed      ! Elapsed inclusive timer duration. [s]

    !--------------------------- Begin Code ---------------------------

    call cpu_time( stop_time )

    timer_id = get_timer_id( timer_name, l_create = .false. )

    if ( stack_depth < 1 ) then
      call timer_fail( "timer_stop called with no active timers: " // trim( timer_name ) )
    end if

    if ( timer_stack(stack_depth)%timer_id /= timer_id ) then
      call timer_fail( "timer_stop nesting mismatch. Expected " &
                       // trim( timer_entries(timer_stack(stack_depth)%timer_id)%name ) &
                       // " but received " // trim( timer_name ) // "." )
    end if

    elapsed = stop_time - timer_stack(stack_depth)%start_time

    timer_entries(timer_id)%inclusive_time = timer_entries(timer_id)%inclusive_time + elapsed
    timer_entries(timer_id)%exclusive_time = timer_entries(timer_id)%exclusive_time &
                                             + elapsed - timer_stack(stack_depth)%child_time

    stack_depth = stack_depth - 1

    if ( stack_depth > 0 ) then
      timer_stack(stack_depth)%child_time = timer_stack(stack_depth)%child_time + elapsed
    end if

    return
  end subroutine cpu_timer_stop
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  function get_timer_id( timer_name, l_create ) result( timer_id )

  ! Description:
  !   Returns the registry index for timer_name, optionally creating a new
  !   entry when the name has not been seen before.
  !-----------------------------------------------------------------------

    implicit none

    !--------------------------- Input Variables ---------------------------
    character(len=*), intent(in) :: &
      timer_name    ! Name of the timer region. [-]

    logical, intent(in) :: &
      l_create    ! Whether to create a missing timer entry. [-]

    !--------------------------- Output Variables ---------------------------
    integer :: &
      timer_id    ! Registry index for timer_name, or 0 if not found. [-]

    !--------------------------- Local Variables ---------------------------
    character(len=:), allocatable :: &
      trimmed_name    ! Timer name without trailing blanks. [-]

    integer :: &
      i    ! Loop index. [-]

    !--------------------------- Begin Code ---------------------------

    trimmed_name = trim( timer_name )

    if ( len( trimmed_name ) < 1 ) then
      call timer_fail( "Timer name must not be empty." )
    end if

    timer_id = 0

    do i = 1, num_timers
      if ( timer_entries(i)%name == trimmed_name ) then
        timer_id = i
        return
      end if
    end do

    if ( .not. l_create ) then
      call timer_fail( "Unknown timer name: " // trimmed_name )
    end if

    call ensure_registry_capacity( num_timers + 1 )

    num_timers = num_timers + 1
    timer_entries(num_timers)%name = trimmed_name
    timer_entries(num_timers)%inclusive_time = 0.0_core_rknd
    timer_entries(num_timers)%exclusive_time = 0.0_core_rknd
    timer_entries(num_timers)%call_count = 0

    timer_id = num_timers

    return
  end function get_timer_id
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  subroutine ensure_registry_capacity( required_capacity )

  ! Description:
  !   Grows the timer registry if required_capacity exceeds its current size.
  !-----------------------------------------------------------------------

    implicit none

    !--------------------------- Input Variables ---------------------------
    integer, intent(in) :: &
      required_capacity    ! Minimum required number of registry slots. [-]

    !--------------------------- Local Variables ---------------------------
    type(timer_entry_t), allocatable :: &
      expanded_entries(:)    ! Temporary expanded timer registry. [-]

    integer :: &
      old_capacity, & ! Existing registry capacity. [-]
      new_capacity    ! Expanded registry capacity. [-]

    !--------------------------- Begin Code ---------------------------

    if ( .not. allocated( timer_entries ) ) then
      allocate( timer_entries( max( registry_chunk_size, required_capacity ) ) )
      return
    end if

    old_capacity = size( timer_entries )

    if ( required_capacity <= old_capacity ) return

    new_capacity = max( required_capacity, old_capacity + registry_chunk_size )
    allocate( expanded_entries( new_capacity ) )

    expanded_entries(1:old_capacity) = timer_entries(1:old_capacity)

    call move_alloc( expanded_entries, timer_entries )

    return
  end subroutine ensure_registry_capacity
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  subroutine ensure_stack_capacity( required_capacity )

  ! Description:
  !   Grows the active timer stack if required_capacity exceeds its current size.
  !-----------------------------------------------------------------------

    implicit none

    !--------------------------- Input Variables ---------------------------
    integer, intent(in) :: &
      required_capacity    ! Minimum required number of stack slots. [-]

    !--------------------------- Local Variables ---------------------------
    type(timer_stack_entry_t), allocatable :: &
      expanded_stack(:)    ! Temporary expanded timer stack. [-]

    integer :: &
      old_capacity, & ! Existing stack capacity. [-]
      new_capacity    ! Expanded stack capacity. [-]

    !--------------------------- Begin Code ---------------------------

    if ( .not. allocated( timer_stack ) ) then
      allocate( timer_stack( max( stack_chunk_size, required_capacity ) ) )
      return
    end if

    old_capacity = size( timer_stack )

    if ( required_capacity <= old_capacity ) return

    new_capacity = max( required_capacity, old_capacity + stack_chunk_size )
    allocate( expanded_stack( new_capacity ) )

    expanded_stack(1:old_capacity) = timer_stack(1:old_capacity)

    call move_alloc( expanded_stack, timer_stack )

    return
  end subroutine ensure_stack_capacity
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  subroutine print_cpu_timer_summary( timer_unit )

  ! Description:
  !   Prints a flat summary of the cpu_time backend timer registry.
  !-----------------------------------------------------------------------

    implicit none

    !--------------------------- Input Variables ---------------------------
    integer, intent(in) :: &
      timer_unit    ! Unit used for timer output. [-]

    !--------------------------- Local Variables ---------------------------
    real( kind = core_rknd ) :: &
      total_exclusive_time, & ! Sum of exclusive time over all timers. [s]
      average_time, &         ! Inclusive time per timer call. [s]
      exclusive_percent       ! Percentage of total exclusive time. [%]

    integer :: &
      i    ! Loop index. [-]

    !--------------------------- Begin Code ---------------------------

    write( unit = timer_unit, fmt = '(a)' ) &
      "========================= CLUBB TIMER SUMMARY ========================="
    write( unit = timer_unit, fmt = '(a)' ) &
      "Backend: cpu_time"

    if ( num_timers < 1 ) then
      write( unit = timer_unit, fmt = '(a)' ) "No timers were recorded."
      write( unit = timer_unit, fmt = '(a)' ) &
        "======================================================================="
      return
    end if

    total_exclusive_time = sum( timer_entries(1:num_timers)%exclusive_time )

    write( unit = timer_unit, fmt = '(a48,1x,a10,1x,a14,1x,a14,1x,a14,1x,a10)' ) &
      "Timer", "Calls", "Inclusive(s)", "Exclusive(s)", "Avg Inc(s)", "Excl Pct"
    write( unit = timer_unit, fmt = '(a48,1x,a10,1x,a14,1x,a14,1x,a14,1x,a10)' ) &
      repeat( "-", display_name_length ), repeat( "-", 10 ), repeat( "-", 14 ), &
      repeat( "-", 14 ), repeat( "-", 14 ), repeat( "-", 10 )

    do i = 1, num_timers

      average_time = timer_entries(i)%inclusive_time &
                     / real( timer_entries(i)%call_count, kind = core_rknd )

      if ( total_exclusive_time > 0.0_core_rknd ) then
        exclusive_percent = 100.0_core_rknd * timer_entries(i)%exclusive_time &
                            / total_exclusive_time
      else
        exclusive_percent = 0.0_core_rknd
      end if

      write( unit = timer_unit, fmt = '(a48,1x,i10,1x,f14.6,1x,f14.6,1x,f14.6,1x,f10.2)' ) &
        timer_display_name( timer_entries(i)%name ), &
        timer_entries(i)%call_count, &
        timer_entries(i)%inclusive_time, &
        timer_entries(i)%exclusive_time, &
        average_time, &
        exclusive_percent

    end do

    write( unit = timer_unit, fmt = '(a)' ) &
      "----------------------------------------------------------------------"
    write( unit = timer_unit, fmt = '(a,f14.6)' ) &
      "Total exclusive timed cpu_time seconds: ", total_exclusive_time
    write( unit = timer_unit, fmt = '(a)' ) &
      "======================================================================="

    return
  end subroutine print_cpu_timer_summary
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  function timer_display_name( timer_name ) result( display_name )

  ! Description:
  !   Returns a fixed-width timer name for the cpu_time summary table.
  !-----------------------------------------------------------------------

    implicit none

    !--------------------------- Input Variables ---------------------------
    character(len=*), intent(in) :: &
      timer_name    ! Full timer name. [-]

    !--------------------------- Output Variables ---------------------------
    character(len=display_name_length) :: &
      display_name    ! Fixed-width display name. [-]

    !--------------------------- Begin Code ---------------------------

    display_name = " "

    if ( len_trim( timer_name ) <= display_name_length ) then
      display_name(1:len_trim( timer_name )) = trim( timer_name )
    else
      display_name(1:display_name_length-3) = timer_name(1:display_name_length-3)
      display_name(display_name_length-2:display_name_length) = "..."
    end if

    return
  end function timer_display_name
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  subroutine require_empty_stack( caller_name )

  ! Description:
  !   Verifies that no cpu_time timer regions are active.
  !-----------------------------------------------------------------------

    implicit none

    !--------------------------- Input Variables ---------------------------
    character(len=*), intent(in) :: &
      caller_name    ! Name of the routine checking the stack. [-]

    !--------------------------- Begin Code ---------------------------

    if ( stack_depth > 0 ) then
      call timer_fail( trim( caller_name ) // " called while timer " &
                       // trim( timer_entries(timer_stack(stack_depth)%timer_id)%name ) &
                       // " is still active." )
    end if

    return
  end subroutine require_empty_stack
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  subroutine reset_cpu_timer_state()

  ! Description:
  !   Clears the cpu_time timer registry and active stack.
  !-----------------------------------------------------------------------

    implicit none

    !--------------------------- Begin Code ---------------------------

    if ( allocated( timer_entries ) ) deallocate( timer_entries )
    if ( allocated( timer_stack ) ) deallocate( timer_stack )

    num_timers = 0
    stack_depth = 0

    return
  end subroutine reset_cpu_timer_state
  !-----------------------------------------------------------------------
#endif

  !-----------------------------------------------------------------------
  function timer_absolute_path( file_path ) result( absolute_path )

  ! Description:
  !   Returns an absolute timer output path, using PWD for relative paths.
  !-----------------------------------------------------------------------

    implicit none

    !--------------------------- Input Variables ---------------------------
    character(len=*), intent(in) :: &
      file_path    ! Timer output path, absolute or relative. [-]

    !--------------------------- Output Variables ---------------------------
    character(len=1024) :: &
      absolute_path    ! Absolute timer output path. [-]

    !--------------------------- Local Variables ---------------------------
    character(len=1024) :: &
      cwd    ! Current working directory from the environment. [-]

    integer :: &
      env_status    ! Status from get_environment_variable. [-]

    !--------------------------- Begin Code ---------------------------

    absolute_path = " "

    if ( file_path(1:1) == "/" ) then
      absolute_path = trim( file_path )
    else
      call get_environment_variable( "PWD", cwd, status = env_status )
      if ( env_status == 0 .and. len_trim( cwd ) > 0 ) then
        absolute_path = trim( cwd ) // "/" // trim( file_path )
      else
        absolute_path = trim( file_path )
      end if
    end if

    return
  end function timer_absolute_path
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  subroutine timer_fail( message )

  ! Description:
  !   Reports an unrecoverable timer API error and stops execution.
  !-----------------------------------------------------------------------

    implicit none

    !--------------------------- Input Variables ---------------------------
    character(len=*), intent(in) :: &
      message    ! Error message. [-]

    !--------------------------- Begin Code ---------------------------

    write( unit = fstderr, fmt = '(a)' ) "CLUBB timer error: " // trim( message )
    error stop "CLUBB timer error"

  end subroutine timer_fail
  !-----------------------------------------------------------------------

end module code_timer_module
