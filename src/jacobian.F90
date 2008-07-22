!-----------------------------------------------------------------------
! $Id: jacobian.F90,v 1.1 2008-07-22 16:04:14 faschinj Exp $

      program jacobian
#ifdef STATS

!     Description:
!     Generates a matrix based on variation between parameter 
!     constants (C1,C2...etc) and variables (cf,rcm,thlm...etc)

!     References:
!     None
!-----------------------------------------------------------------------

      use hoc, only:  & 
          hoc_model ! Procedure(s)
      use param_index, only:  & 
          nparams ! Variable(s)
      use parameters, only:  & 
          params_list,  & ! Variable(s)
          read_parameters ! Procedure(s) 
      use constants, only:  & 
          fstdout,  & ! Variable(s) 
          fstderr
      use grads_common, only:  & 
          grads_average_interval,  & ! Procedure(s) 
          grads_zlvl
      use stats_variables, only:  & 
          fname_zt,  & ! Variable(s) 
          fname_zm
      use error_code, only:  & 
          clubb_no_error,  & ! Variable(s)
          fatal_error ! Procedure(s)

      implicit none

      ! Variable derived types
!----------------------------------------------------------------------
        type param_array

          integer :: entries ! Total constant parameters

          character(len=11), pointer :: name(:)

          real, pointer :: value(:)

        end type param_array
!----------------------------------------------------------------------
        type variable_array

          integer ::  & 
          nz,      & ! Z dimension [grid boxes] 
          entries ! Total variables

          character(len=12), pointer :: name(:)

          real, pointer :: value(:,:) ! (1:nz, entries)

        end type variable_array
!----------------------------------------------------------------------

        ! External
        intrinsic sum, transfer, abs, int, trim

        ! Local Constants
        integer, parameter :: & 
        nvarzt = 14, & 
        nvarzm = 40

        ! 42 must be changed to be equal to nparams 
        character(len=12), parameter :: write_format = "(42(e18.10))"

!       character, parameter :: delta   = greek_'Δ' ! Doesn't work
        character, parameter :: delta   = 'D'

!-----------------------------------------------------------------------

        ! Other namelist Variables
        real ::  & 
        delta_factor ! factor constant parameters are multiplied by

        integer, dimension(10) :: & 
        times ! Times to read in [GraDS output file units]

        character(len=50) ::  & 
        run_file ! namelist for case being run

        logical :: & 
        use_standard_vars ! use standard constant parameters


!     ! Types to hold GrADS variables and parameter constants
        type (param_array) :: hoc_params

        type (variable_array) ::  & 
        var1zt,  & ! thermo grid GRaDS results [units vary]
        var2zt,  & ! thermo grid GRaDS results [units vary]
        var1zm,  & ! momentum grid GRaDS results [units vary]
        var2zm  ! momentum grid GRaDS results [units vary]


        real, dimension(nparams, nvarzt+nvarzm) :: jmatrix
        real, dimension(nparams, nvarzt+nvarzm) :: impact_matrix
        real, dimension(nparams, nvarzt+nvarzm) :: fc_impact_matrix

        integer :: nzt ! Thermo grid levels
        integer :: nzm ! Momentum grid levels
        integer :: alloc_stat ! Det. whether array allocation worked

        real :: tmp_param ! Temp variable

        integer :: i, j, k ! loop variables
        integer :: nanbits ! Holds a NaN value for purposes of output

        integer :: err_code ! Determines whether a run went unstable

        ! Namelists
        namelist /jcbn_nml/  & 
          run_file, times, delta_factor, use_standard_vars

!-----------------------------------------------------------------------

        ! Initialize data
        data nanbits /Z"7F800000"/

        times(1:10) = 0

        err_code = clubb_no_error

        allocate( hoc_params%value( nparams ),  & 
                  hoc_params%name( nparams ), & 
                  stat=alloc_stat )
        if (alloc_stat /= 0 ) stop "allocate failed"

        hoc_params%entries = nparams

        ! Read namelists
        open( unit=10, file='jacobian.in', status='old')
        read( unit=10, nml=jcbn_nml )

        close( unit=10 )

        if ( .not. use_standard_vars ) then
          call read_parameters( 10, 'jacobian.in', hoc_params%value )

        else
          call read_parameters & 
               ( 10, "../cmp_stats/std_const/std_const_nml.in", & 
                 hoc_params%value )

        end if

        hoc_params%name(1:nparams) = params_list(1:nparams)
         
        write(unit=fstdout,fmt='(a11,2a12)')  & 
          "Parameter  ", "Initial     ", "Varied      "

        do i = 1, hoc_params%entries, 1
          write(unit=*,fmt='(a11,2f12.5)') hoc_params%name(i),  & 
            hoc_params%value(i), hoc_params%value(i) * delta_factor
        end do

        call hoc_model  & 
             ( hoc_params%value(:), run_file, err_code,  & 
              .false., .false. )

        if ( fatal_error(err_code) ) then
        
           stop "Initial run wasn't valid"
        endif
        ! Obtain nz from the generated GrADS files

        nzt = grads_zlvl( trim( fname_zt )//".ctl" )
        nzm = grads_zlvl( trim( fname_zm )//".ctl" )

        ! Initialize the structures holding the variables

        allocate( var1zt%value(nzt, nvarzt), & 
                  var2zt%value(nzt, nvarzt), & 
                  var1zt%name(nvarzt), & 
                  var2zt%name(nvarzt), & 
                  stat=alloc_stat )

        if (alloc_stat /= 0 ) stop "allocate failed"

        var1zt%entries = nvarzt 
        var1zt%nz      = nzt
        var2zt%entries = nvarzt 
        var2zt%nz      = nzt

        allocate( var1zm%value(nzm, nvarzm), & 
                  var2zm%value(nzm, nvarzm), & 
                  var1zm%name(nvarzm), & 
                  var2zm%name(nvarzm), & 
                  stat = alloc_stat )

        if (alloc_stat /= 0 ) stop "allocate failed"

        var1zm%entries = nvarzm 
        var1zm%nz      = nzm
        var2zm%entries = nvarzm 
        var2zm%nz      = nzm



        var1zt%name(1:nvarzt) =  & 
        (/"cf          ", "rcm         ", "rtm         ",  & 
          "thlm        ", "um          ", "vm          ", & 
          "wp3         ", "wp3_ta      ", "wp3_tp      ", & 
          "wp3_bp      ", "wp3_cl      ", & 
          "wp3_pr1     ", "wp3_pr2     ", "wp3_dp1     "/)

        var1zm%name(1:nvarzm) =  & 
        (/"wp2         ", "rtp2        ", "thlp2       ",  & 
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

        do i = 1, hoc_params%entries
          tmp_param = hoc_params%value(i)
          hoc_params%value(i) = hoc_params%value(i) * delta_factor

          call hoc_model & 
            ( hoc_params%value(:), trim( run_file ), err_code, & 
               .false., .false. )

          ! Print a period so the user knows something is happening
          write(unit=fstdout, fmt='(a1)', advance='no') "."

          if ( fatal_error(err_code) ) then
         
            ! Pos. Infinity bit pattern 
            jmatrix(i, :) & 
            = transfer( int( nanbits ), jmatrix(1,1) )
            hoc_params%value(i) = tmp_param
            cycle
          end if

          ! Set var2 results with results from altering the constants

          call getvariables( var2zt, trim( fname_zt )//".ctl" )
          call getvariables( var2zm, trim( fname_zm )//".ctl" )

          do j = 1, nvarzt
            impact_matrix(i, j) =  & 
            sum(  var2zt%value(1:nzt,j) - var1zt%value(1:nzt,j) ) & 
            / nzt

            fc_impact_matrix(i, j) =  & 
            sum(  var2zt%value(1:nzt, j) - var1zt%value(1:nzt, j) ) & 
            / sum( abs( var1zt%value(1:nzt, j) ) )
          end do

          do j = 1, nvarzt
            jmatrix(i, j) = impact_matrix(i, j)  & 
                          / ( hoc_params%value(i) - tmp_param )
          end do

          do j = 1, nvarzm
            impact_matrix(i, j+nvarzt)  & 
            = sum( var2zm%value(1:nzm, j) - var1zm%value(1:nzm, j) ) & 
            / nzm

            fc_impact_matrix(i, j+nvarzt)  & 
            = sum( var2zm%value(1:nzm, j) - var1zm%value(1:nzm, j) ) & 
            / sum( abs( var1zm%value(1:nzm, j) ) )
          end do

          do j = 1, nvarzm
            jmatrix(i, j+nvarzt) =  & 
            impact_matrix(i, j+nvarzt)  & 
             / (hoc_params%value(i) - tmp_param)
          end do

          hoc_params%value(i) = tmp_param ! Set parameter back

        end do !i = 1..hoc_params%entries

        ! Output Results to the terminal
        do i = 1, hoc_params%entries

            do j = 1, var2zt%entries
              write(unit=*,fmt='(3(a,e10.4))')  & 
              delta//var2zt%name(j)//"/" & 
              //delta//hoc_params%name(i)//" = ", jmatrix(i, j), & 
              " impact: ", impact_matrix(i, j), & 
              " fc imp: ", fc_impact_matrix(i, j)

            end do

            do j = 1, var2zm%entries
              write(unit=*,fmt='(3(a,e10.4))')  & 
              delta//var2zm%name(j)//"/" & 
              //delta//hoc_params%name(i)//" = ", jmatrix(i, j+nvarzt), & 
              " impact: ", impact_matrix(i, j+nvarzt), & 
              " fc imp: ", fc_impact_matrix(i, j+nvarzt)

            end do

          write(unit=fstderr,fmt=*) ""

        end do ! 1..hoc_params%entries

        ! Output results to an ASCII file

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

        stop "Program exited normally"

        contains
!-----------------------------------------------------------------------
        subroutine getvariables( varray, fname_zx ) 

!       Description:
!       Returns a variable_array structure over namelist defined
!       time intervals

!       References:
!       None
!-----------------------------------------------------------------------
 
        implicit none

        ! Input
        type (variable_array), intent(inout) :: varray

        character(len=*), intent(in) :: fname_zx

        ! Local Variable
        logical :: error


    
        do k = 1, varray%entries, 1

          varray%value(1:varray%nz, k) =  & 
          grads_average_interval & 
          ( fname_zx, varray%nz, times(:), varray%name(k), 1, error )

          if ( error ) then
            write(unit=fstderr,fmt=*) "Error in reading"//varray%name(i)
            stop
          end if

        end do

        return
        end subroutine getvariables 

!-----------------------------------------------------------------------
#else
        stop "Compile all files with -DSTATS to use this program."
#endif /*STATS*/

        end program jacobian
