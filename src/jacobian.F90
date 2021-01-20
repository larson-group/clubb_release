!-----------------------------------------------------------------------
! $Id$

program jacobian


! Description:
!   Generates a matrix based on variation between parameter
!   constants (C1,C2...etc) and variables (cloud_frac,rcm,thlm...etc)

! References:
!   None
!-----------------------------------------------------------------------

  use clubb_driver, only:  & 
      run_clubb ! Procedure(s)

  use parameter_indices, only:  & 
      nparams ! Variable(s)

  use parameters_tunable, only:  & 
      params_list,  & ! Variable(s)
      read_parameters ! Procedure(s)

  use constants_clubb, only:  & 
      fstdout,  & ! Constant(s) 
      fstderr

  use stat_file_utils, only:  & 
      stat_file_average_interval,  & ! Procedure(s) 
      stat_file_num_vertical_levels, &
      stat_file_vertical_levels

  use stats_variables, only:  & 
      fname_zt,  & ! Variable(s) 
      fname_zm

  use error_code, only: &
        clubb_at_least_debug_level,  & ! Procedure
        err_code,                    & ! Error Indicator
        clubb_fatal_error              ! Constant

  use parameters_model, only: &
    PosInf ! Variable(s)

  use clubb_precision, only: &
    core_rknd ! Variable(s)

  implicit none

  !----------------------------------------------------------------------------
  ! Variable derived types
  !----------------------------------------------------------------------------
  type param_array

    integer :: entries ! Total tunable parameters

    character(len=32), pointer :: name(:)

    real( kind = core_rknd ), pointer :: value(:)

  end type param_array
  !----------------------------------------------------------------------------
  type variable_array

    integer ::  & 
      nz,   & ! Z dimension [grid boxes] 
      entries ! Total variables

    real( kind = core_rknd ), pointer, dimension(:) :: z

    character(len=12), pointer :: name(:)

    real( kind = core_rknd ), pointer, dimension(:,:) :: value ! (1:nz, entries)

  end type variable_array
  !-----------------------------------------------------------------------------

  ! External
  intrinsic sum, transfer, abs, int, trim

  ! Constant Parameters
  integer, parameter :: & 
    nvarzt = 15, & 
    nvarzm = 40

! character, parameter :: delta   = 'Δ' ! Only works on unicode terminals

  character, parameter :: delta   = 'D'

  ! Local Variables
  integer, dimension(10) :: & 
    times ! Times to read in [GraDS output file units]

  ! Types to hold GrADS variables and parameter constants
  type (param_array) :: clubb_params

  type (variable_array) ::  & 
    var1zt,  & ! Thermo grid GrADS results   [units vary]
    var2zt,  & ! Thermo grid GrADS results   [units vary]
    var1zm,  & ! Momentum grid GrADS results [units vary]
    var2zm     ! Momentum grid GrADS results [units vary]


  real( kind = core_rknd ), dimension(nparams, nvarzt+nvarzm) :: &
    jmatrix, &           ! Jacobian matrix
    impact_matrix, &     ! Impact matrix
    fc_impact_matrix     !

  ! XX must be changed to be equal to nparams
  character(len=12), save :: write_format = "(XX(E18.10))"

  integer :: &
    nzt, &        ! Thermo grid levels
    nzm, &        ! Momentum grid levels
    alloc_stat, & ! Det. whether array allocation worked
    i, j, k       ! loop variables

  real( kind = core_rknd ) ::  & 
    delta_factor, & ! Factor that tunable parameters are multiplied by
    tmp_param       ! Temporary variable

  logical :: & 
    l_use_standard_vars ! Whether to use the standard tunable parameters

  ! Namelists
  namelist /jcbn_nml/  & 
    times, delta_factor, l_use_standard_vars

!-----------------------------------------------------------------------

  ! ---- Begin Code ----

  ! Use an internal file write to specify the write format for the jacobian_matrix.txt
  ! and impact_matrix.txt files.
  write(unit=write_format(2:3), fmt='(i2.2)') nparams

  times(1:10) = 0

  allocate( clubb_params%value( nparams ),  & 
            clubb_params%name( nparams ), & 
            stat=alloc_stat )
  if (alloc_stat /= 0 ) stop "allocate failed"

  clubb_params%entries = nparams

  ! Read namelists
  open( unit=10, file='jacobian.in', status='old')
  read( unit=10, nml=jcbn_nml )

  close( unit=10 )

  if ( .not. l_use_standard_vars ) then
    call read_parameters( 10, 'jacobian.in', clubb_params%value )

  else
    call read_parameters & 
         ( 10, "", clubb_params%value )

  end if

  clubb_params%name(1:nparams) = params_list(1:nparams)

  write(unit=fstdout,fmt='(a27,2a12)')  & 
    "Parameter", "Initial", "Varied"

  do i = 1, clubb_params%entries, 1
    write(unit=*,fmt='(a27,2f12.5)') trim( clubb_params%name(i) ),  & 
      clubb_params%value(i), clubb_params%value(i) * delta_factor
  end do

  call run_clubb  & 
       ( clubb_params%value(:), 'jacobian.in', .false. )

  if ( clubb_at_least_debug_level( 0 ) ) then
    if ( err_code == clubb_fatal_error ) then
      stop "The initial set of parameters caused a fatal error."
    end if
  end if

  ! Obtain number of vertical levels from the generated GrADS files

  nzt = stat_file_num_vertical_levels( "thlm", "../output/"//trim( fname_zt )//".ctl" )
  nzm = stat_file_num_vertical_levels( "thlm", "../output/"//trim( fname_zm )//".ctl" )

  ! Initialize the structures holding the variables

  allocate( var1zt%value(nzt, nvarzt), & 
            var2zt%value(nzt, nvarzt), & 
            var1zt%name(nvarzt), & 
            var2zt%name(nvarzt), & 
            var1zt%z(nzt), & 
            var2zt%z(nzt), & 
            stat=alloc_stat )

  if (alloc_stat /= 0 ) stop "allocate failed"

  var1zt%entries = nvarzt
  var1zt%nz      = nzt
  var2zt%entries = nvarzt
  var2zt%nz      = nzt

  var1zt%z = stat_file_vertical_levels &
    ( var1zt%name(1), "../output/"//trim( fname_zt )//".ctl", nzt )
  var2zt%z = stat_file_vertical_levels &
    ( var1zt%name(1), "../output/"//trim( fname_zt )//".ctl", nzt )

  allocate( var1zm%value(nzm, nvarzm), & 
            var2zm%value(nzm, nvarzm), & 
            var1zm%name(nvarzm), & 
            var2zm%name(nvarzm), & 
            var1zm%z(nzm), & 
            var2zm%z(nzm), & 
            stat = alloc_stat )

  if (alloc_stat /= 0 ) stop "allocate failed"

  var1zm%entries = nvarzm
  var1zm%nz      = nzm
  var2zm%entries = nvarzm
  var2zm%nz      = nzm

  var1zm%z = stat_file_vertical_levels &
    ( var1zm%name(1), "../output/"//trim( fname_zm )//".ctl", nzm )
  var2zm%z = stat_file_vertical_levels &
    ( var1zm%name(1), "../output/"//trim( fname_zm )//".ctl", nzm )

  var1zt%name(1:nvarzt) =  & 
  (/"cloud_frac  ", "rcm         ", "rtm         ", & 
    "thlm        ", "um          ", "vm          ", & 
    "wp3         ", "wp3_ta      ", "wp3_tp      ", & 
    "wp3_bp1     ", "wp3_pr_turb ", "wp3_cl      ", & 
    "wp3_pr1     ", "wp3_pr2     ", "wp3_dp1     "/)

  var1zm%name(1:nvarzm) =  & 
  (/"wp2         ", "rtp2        ", "thlp2       ", & 
    "rtpthlp     ", "wprtp       ", "wpthlp      ", & 
    "wp2_ta      ", "wp2_bp      ", "wp2_pr2     ", & 
    "wp2_pr3     ", "wp2_dp1     ", "wp2_dp2     ", & 
    "wp2_cl      ",  & 
    "wprtp_ta    ", "wprtp_tp    ", "wprtp_bp    ",  & 
    "wprtp_pr1   ", "wprtp_pr2   ",  & 
    "wprtp_pr3   ", "wprtp_dp1   ", & 
    "wpthlp_ta   ", "wpthlp_tp   ", "wpthlp_bp   ",  & 
    "wpthlp_pr1  ", "wpthlp_pr2  ",  & 
    "wpthlp_pr3  ", "wpthlp_dp1  ", & 
    "rtp2_ta     ", "rtp2_tp     ",  & 
    "rtp2_dp1    ", "rtp2_dp2    ",  & 
    "thlp2_ta    ", "thlp2_tp    ",  & 
    "thlp2_dp1   ", "thlp2_dp2   ", & 
    "rtpthlp_ta  ", "rtpthlp_tp1 ", "rtpthlp_tp2 ",  & 
    "rtpthlp_dp1 ", "rtpthlp_dp2 "/)

  var2zt%name(1:nvarzt) = var1zt%name(1:nvarzt)

  var2zm%name(1:nvarzm) = var1zm%name(1:nvarzm)

  ! Set var1 fields with initial run results

  call getvariables( var1zt, trim( fname_zt )//".ctl" )

  call getvariables( var1zm, trim( fname_zm )//".ctl" )

  do i = 1, clubb_params%entries
    tmp_param = clubb_params%value(i)
    clubb_params%value(i) = clubb_params%value(i) * delta_factor

    call run_clubb & 
      ( clubb_params%value(:), 'jacobian.in', .false. )

    ! Print a period so the user knows something is happening
    write(unit=fstdout, fmt='(a1)', advance='no') "."

    if ( err_code == clubb_fatal_error ) then

      ! Pos. Infinity bit pattern
      jmatrix(i,:) = real(PosInf, kind = core_rknd)
      clubb_params%value(i) = tmp_param
      cycle
    end if

    ! Set var2 results with results from altering the constants

    call getvariables( var2zt, trim( fname_zt )//".ctl" )
    call getvariables( var2zm, trim( fname_zm )//".ctl" )

    do j = 1, nvarzt
      impact_matrix(i, j) =  & 
      sum(  var2zt%value(1:nzt,j) - var1zt%value(1:nzt,j) ) & 
      / real( nzt, kind = core_rknd )

      fc_impact_matrix(i, j) =  & 
      sum(  var2zt%value(1:nzt, j) - var1zt%value(1:nzt, j) ) & 
      / sum( abs( var1zt%value(1:nzt, j) ) )
    end do

    do j = 1, nvarzt
      jmatrix(i, j) = impact_matrix(i, j)  & 
                    / ( clubb_params%value(i) - tmp_param )
    end do

    do j = 1, nvarzm
      impact_matrix(i, j+nvarzt)  & 
      = sum( var2zm%value(1:nzm, j) - var1zm%value(1:nzm, j) ) & 
      / real( nzm, kind = core_rknd )

      fc_impact_matrix(i, j+nvarzt)  & 
      = sum( var2zm%value(1:nzm, j) - var1zm%value(1:nzm, j) ) & 
      / sum( abs( var1zm%value(1:nzm, j) ) )
    end do

    do j = 1, nvarzm
      jmatrix(i, j+nvarzt) =  & 
      impact_matrix(i, j+nvarzt)  & 
       / (clubb_params%value(i) - tmp_param)
    end do

    clubb_params%value(i) = tmp_param ! Set parameter back

  end do !i = 1..clubb_params%entries

  ! Write a newline
  write(6,*) ""

  ! Output Results to the terminal
  do i = 1, clubb_params%entries

    do j = 1, var2zt%entries
      write(unit=*,fmt='(3(a,e11.4))')  & 
      delta//var2zt%name(j)//"/" & 
      //delta//clubb_params%name(i)//" = ", jmatrix(i, j), & 
      " impact: ", impact_matrix(i, j), & 
      " fc imp: ", fc_impact_matrix(i, j)

    end do

    do j = 1, var2zm%entries
      write(unit=*,fmt='(3(a,e11.4))')  & 
      delta//var2zm%name(j)//"/" & 
      //delta//clubb_params%name(i)//" = ", jmatrix(i, j+nvarzt), & 
      " impact: ", impact_matrix(i, j+nvarzt), & 
      " fc imp: ", fc_impact_matrix(i, j+nvarzt)

    end do

    write(unit=fstderr,fmt=*) ""

  end do ! 1..clubb_params%entries

  ! Output results to ASCII text files

  ! Jacobian Matrix
  open(unit=20, file="jacobian_matrix.txt", action="write" )

  do j = 1, nvarzm + nvarzt
    write(unit=20, fmt=write_format) jmatrix(:,j)
  end do

  close(unit=20)

  ! Impact Matrix
  open(unit=20, file="impact_matrix.txt", action="write" )

  do j = 1, nvarzm + nvarzt
    write(unit=20, fmt=write_format) impact_matrix(:,j)
  end do

  close(unit=20)

  ! Deallocate memory
  deallocate( var1zt%value, var2zt%value, & 
              var1zt%name, var2zt%name, & 
              var1zt%z, var2zt%z )

  deallocate( var1zm%value, var2zm%value, & 
              var1zm%name, var2zm%name, & 
              var1zm%z, var2zm%z )

  stop "Program exited normally"

  contains
!-----------------------------------------------------------------------
  subroutine getvariables( varray, fname_zx )

! Description:
!   Returns a variable_array structure over namelist defined
!   time intervals

! References:
!   None
!-----------------------------------------------------------------------

    implicit none

    ! Input Variables
    type (variable_array), intent(inout) :: varray

    character(len=*), intent(in) :: fname_zx

    ! Local Variable
    logical :: l_error

    ! ---- Begin Code ----

    do k = 1, varray%entries, 1

      varray%value(1:varray%nz, k) =  & 
      stat_file_average_interval & 
      ( "../output/"//fname_zx, varray%nz, times(:), varray%name(k), &
        varray%z, 1, l_error )

      if ( l_error ) then
        write(unit=fstderr,fmt=*) "Error in reading"//varray%name(i)
        stop
      end if

    end do

    return
  end subroutine getvariables

!-----------------------------------------------------------------------
end program jacobian
