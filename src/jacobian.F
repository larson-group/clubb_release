!-----------------------------------------------------------------------
! $Id: jacobian.F,v 1.31 2006-07-27 01:46:48 dschanen Exp $

      program jacobian

!     Description:
!     Generates a matrix based on variation between parameter 
!     constants (C1,C2...etc) and variables (cf,rcm,thlm...etc)

!     References:
!     None
!-----------------------------------------------------------------------

      use hoc, only: hoc_model
      use grads_common, only: grads_average_interval, grads_zlvl
      use statistics, only: fname_zt, fname_zm

      implicit none

      ! Variable derived types
!----------------------------------------------------------------------
        type constant_array

          integer entries ! Total constant parameters

          character(len=11), pointer :: name(:)

          real, pointer :: constant(:)

        end type constant_array
!----------------------------------------------------------------------
        type variable_array

          integer :: 
     .    nz,     ! Z dimension [grid boxes] 
     .    entries ! Total variables

          character(len=12), pointer :: name(:)

          real, pointer :: variable(:,:) ! (1:nz, entries)

        end type variable_array
!----------------------------------------------------------------------

        ! External
        intrinsic sum, transfer, abs, int, trim

        ! Local Constants
        integer, parameter ::
     .  nvarzt = 14,
     .  nvarzm = 40,
     .  nconst = 35

        ! 35 must be changed to be equal to nconst
        character(len=12), parameter :: write_format = "(35(e18.10))"

!       character, parameter :: delta   = greek_'Δ' ! Doesn't work
        character, parameter :: delta   = 'D'

!-----------------------------------------------------------------------

        ! Model constant parameters 
        real ::
     .  C1,         ! Low Skewness in C1 Skewness Function.
     .  C1b,        ! High Skewness in C1 Skewness Function.
     .  C1c,        ! Degree of Slope of C1 Skewness Function.
     .  C2,         ! Low Skewness in C2 Skewness Function.
     .  C2b,        ! High Skewness in C2 Skewness Function.  
     .  C2c,        ! Degree of Slope of C2 Skewness Function.
     .  C4,         ! Currently Not Used in the Model.
     .  C5,         ! Constant from Chris Golaz.
     .  C6rt,       ! Low Skewness in C6rt Skewness Function.
     .  C6rtb,      ! High Skewness in C6rt Skewness Function.
     .  C6rtc,      ! Degree of Slope of C6rt Skewness Function.
     .  C6thl,      ! Low Skewness in C6thl Skewness Function.
     .  C6thlb,     ! High Skewness in C6thl Skewness Function.
     .  C6thlc,     ! Degree of Slope of C6thl Skewness Function.
     .  C7,         ! Low Skewness in C7 Skewness Function.
     .  C7b,        ! High Skewness in C7 Skewness Function.
     .  C7c,        ! Degree of Slope of C7 Skewness Function.
     .  C8,         ! Coefficient #1 in C8 Skewness Equation.
     .  C8b,        ! Coefficient #2 in C8 Skewness Equation.  From Golaz.
     .  C10,        ! Currently Not Used in the Model.
     .  C11,        ! Low Skewness in C11 Skewness Function.
     .  C11b,       ! High Skewness in C11 Skewness Function.  From Golaz.
     .  C11c,       ! Degree of Slope of C11 Skewness Function.  From Golaz.
     .  C12,        ! Constant from Chris Golaz.
     .  nu1,        ! Previously determined eddy viscosity.
     .  nu2,        ! Previously determined eddy viscosity.
     .  nu6,        ! Previously determined eddy viscosity.
     .  nu8,        ! Previously determined eddy viscosity.
     .  nu_r,       ! Tuned to LES from Surface Value of DYCOMS2 RF02 DS.
     .  beta,       ! Constant from Chris Golaz.
     .  gamma_coef, ! Low Skewness in gamma coefficient Skewness Function.
     .  gamma_coefb,! High Skewness in gamma coefficient Skewness Function.
     .  gamma_coefc,! Degree of Slope of gamma coefficient Skewness Function.
     .  c_K,        ! Constant from Chris Golaz.
     .  mu          ! Constant from Chris Golaz.

      ! These aren't varied in matrix generation
        real ::
     .  taumin,    ! Previously determined value.
     .  taumax,    ! Previously determined value.
     .  lmin_coef  ! Previously determined value.

        ! Other namelist Variables
        real :: 
     .  delta_factor ! factor constant parameters are multiplied by

        integer, dimension(10) ::
     .  times ! Times to read in [GraDS output file units]

        character(len=50) :: 
     .  run_file ! namelist for case being run

        logical ::
     .  use_standard_vars ! use standard constant parameters


!     ! Types to hold GRaDS variables and parameter constants
        type (constant_array) :: ca1

        type (variable_array) :: 
     .  var1zt, ! thermo grid GRaDS results [units vary]
     .  var2zt, ! thermo grid GRaDS results [units vary]
     .  var1zm, ! momentum grid GRaDS results [units vary]
     .  var2zm  ! momentum grid GRaDS results [units vary]


        real, dimension(nconst, nvarzt+nvarzm) :: jmatrix
        real, dimension(nconst, nvarzt+nvarzm) :: impact_matrix
        real, dimension(nconst, nvarzt+nvarzm) :: fc_impact_matrix

        integer :: nzt ! Thermo grid levels
        integer :: nzm ! Momentum grid levels
        integer :: alloc_stat ! Det. whether array allocation worked

        real :: tmp_constant ! Temp variable

        integer :: i, j, k ! loop variables
        integer :: nanbits ! Holds a NaN value for purposes of output

        logical :: isValid ! Determines whether a run went unstable

        ! Namelists
        namelist /jcbn_nml/ 
     .  run_file, times, delta_factor, use_standard_vars

        namelist /initvars/ 
     .  C1, C1b, C1c, C2, C2b, C2c, C4, C5,
     .  C6rt, C6rtb, C6rtc, C6thl, C6thlb, C6thlc,
     .  C7, C7b, C7c, C8, C8b, C10, 
     .  C11, C11b, C11c, C12,
     .  nu1, nu2, nu6, nu8, nu_r, beta, gamma_coef, 
     .  gamma_coefb, gamma_coefc,
     .  c_K, lmin_coef, taumin, taumax, mu

!-----------------------------------------------------------------------

        ! Initialize data
        data nanbits /Z"7F800000"/

        times(1:10) = 0
        isValid     = .true.


        allocate( ca1%constant( nconst ), ca1%name( nconst ),
     .            stat=alloc_stat )
        if (alloc_stat /= 0 ) stop "allocate failed"

        ca1%entries = nconst

        ! Read namelists
        open( unit=10, file='jacobian.in', status='old')
        read( unit=10, nml=jcbn_nml )

        if ( .not. use_standard_vars ) then
          read( unit=10, nml=initvars )
        else
          open( unit=20, file="../cmp_stats/std_const/std_const_nml.in",
     .          status="old" )

          read( unit=20, nml=initvars )
          close( unit=20 )
        end if
        close( unit=10 )
      
        ca1%constant(1:nconst) = 
     .  (/C1, C1b, C1c, C2, C2b, C2c,
     .    C4, C5, C6rt, C6rtb, C6rtc, C6thl,
     .    C6thlb, C6thlc, C7, C7b, C7c, C8, C8b,
     .    C10, C11, C11b, C11c, C12,
     .    nu1, nu2, nu6, nu8, nu_r, beta,
     .    gamma_coef, gamma_coefb, gamma_coefc, 
     .    c_K, mu/)     
        ca1%name(1:nconst) =
     .  (/"C1         ", "C1b        ", "C1c        ",
     .    "C2         ", "C2b        ", "C2c        ",
     .    "C4         ", "C5         ",
     .    "C6rt       ", "C6rtb      ", "C6rtc      ",
     .    "C6thl      ", "C6thlb     ", "C6thlc     ",
     .    "C7         ", "C7b        ", "C7c        ",
     .    "C8         ", "C8b        ",
     .    "C10        ", "C11        ", "C11b       ",
     .    "C11c       ", "C12        ",
     .    "nu1        ", "nu2        ", 
     .    "nu6        ", "nu8        ", "nu_r       ",
     .    "beta       ", "gamma_coef ", "gamma_coefb", "gamma_coefc",
     .    "c_K        ", "mu         "/)

        write(unit=*,fmt='(3a12)') "Parameter", "Initial", "Varied"

        do i = 1, ca1%entries
          write(unit=*,fmt='(a12,2f12.5)') ca1%name(i), 
     .    ca1%constant(i), ca1%constant(i) * delta_factor
        end do

        call hoc_model 
     .  ( ca1%constant(1), ca1%constant(2), 
     .    ca1%constant(3), ca1%constant(4), 
     .    ca1%constant(5), ca1%constant(6), 
     .    ca1%constant(7), ca1%constant(8), 
     .    ca1%constant(9), ca1%constant(10), 
     .    ca1%constant(11), ca1%constant(12), 
     .    ca1%constant(13), ca1%constant(14), 
     .    ca1%constant(15), ca1%constant(16), 
     .    ca1%constant(17), ca1%constant(18), 
     .    ca1%constant(19), ca1%constant(20),
     .    ca1%constant(21), ca1%constant(22),
     .    ca1%constant(23), ca1%constant(24),
     .    ca1%constant(25), ca1%constant(26),
     .    ca1%constant(27), ca1%constant(28),
     .    ca1%constant(29), ca1%constant(30),
     .    ca1%constant(31), ca1%constant(32),
     .    ca1%constant(33), ca1%constant(34),
     .    lmin_coef, taumin, taumax, 
     .    ca1%constant(35), trim( run_file ),
     .    isValid, .false. )

        if ( .not. isValid ) stop "Initial run wasn't valid"

        ! Obtain nz from the generated GrADS files

        nzt = grads_zlvl( trim(fname_zt)//".ctl" )
        nzm = grads_zlvl( trim(fname_zm)//".ctl" )


        ! Initialize the structures holding the variables

        allocate( var1zt%variable( nzt, nvarzt ),
     .            var2zt%variable( nzt, nvarzt ),
     .            var1zt%name( nvarzt ),
     .            var2zt%name( nvarzt ),
     .            stat=alloc_stat )

        if (alloc_stat /= 0 ) stop "allocate failed"

        var1zt%entries = nvarzt 
        var1zt%nz      = nzt
        var2zt%entries = nvarzt 
        var2zt%nz      = nzt

        allocate( var1zm%variable( nzm, nvarzm ),
     .            var2zm%variable( nzm, nvarzm ),
     .            var1zm%name( nvarzm ),
     .            var2zm%name( nvarzm ),
     .            stat = alloc_stat )

        if (alloc_stat /= 0 ) stop "allocate failed"

        var1zm%entries = nvarzm 
        var1zm%nz      = nzm
        var2zm%entries = nvarzm 
        var2zm%nz      = nzm



        var1zt%name(1:nvarzt) = 
     .  (/"cf          ", "rcm         ", "rtm         ", 
     .    "thlm        ", "um          ", "vm          ",
     .    "wp3         ", "wp3_ta      ", "wp3_tp      ",
     .    "wp3_bp      ", "wp3_cl      ",
     .    "wp3_pr1     ", "wp3_pr2     ", "wp3_dp1     "/)

        var1zm%name(1:nvarzm) = 
     .  (/"wp2         ", "rtp2        ", "thlp2       ", 
     .    "rtpthlp     ", "wprtp       ", "wpthlp      ",
     .    "wp2_ta      ", "wp2_bp      ", "wp2_pr2     ",
     .    "wp2_pr3     ", "wp2_dp1     ", "wp2_dp2     ",
     .    "wp2_cl      ", 
     .    "wprtp_ta    ", "wprtp_tp    ", "wprtp_bp    ", 
     .    "wprtp_pr1   ", "wprtp_pr2   ", 
     .    "wprtp_pr3   ", "wprtp_dp1   ",
     .    "wpthlp_ta   ", "wpthlp_tp   ", "wpthlp_bp   ", 
     .    "wpthlp_pr1  ", "wpthlp_pr2  ", 
     .    "wpthlp_pr3  ", "wpthlp_dp1  ",
     .    "rtp2_ta     ", "rtp2_tp     ", 
     .    "rtp2_dp1    ", "rtp2_dp2    ", 
     .    "thlp2_ta    ", "thlp2_tp    ", 
     .    "thlp2_dp1   ", "thlp2_dp2   ",
     .    "rtpthlp_ta  ", "rtpthlp_tp1 ", "rtpthlp_tp2 ", 
     .    "rtpthlp_dp1 ", "rtpthlp_dp2 "/)

        var2zt%name(1:nvarzt) = var1zt%name(1:nvarzt) 

        var2zm%name(1:nvarzm) = var1zm%name(1:nvarzm) 

        ! Set var1 fields with initial run results

        call getvariables( var1zt, trim(fname_zt)//".ctl" )

        call getvariables( var1zm, trim(fname_zm)//".ctl" )

        do i = 1, ca1%entries
          tmp_constant    = ca1%constant(i)
          ca1%constant(i) = ca1%constant(i) * delta_factor

          call hoc_model
     .         ( ca1%constant(1),  ca1%constant(2), 
     .           ca1%constant(3),  ca1%constant(4), 
     .           ca1%constant(5),  ca1%constant(6), 
     .           ca1%constant(7),  ca1%constant(8), 
     .           ca1%constant(9),  ca1%constant(10), 
     .           ca1%constant(11), ca1%constant(12), 
     .           ca1%constant(13), ca1%constant(14), 
     .           ca1%constant(15), ca1%constant(16), 
     .           ca1%constant(17), ca1%constant(18), 
     .           ca1%constant(19), ca1%constant(20),
     .           ca1%constant(21), ca1%constant(22),
     .           ca1%constant(23), ca1%constant(24),
     .           ca1%constant(25), ca1%constant(26),
     .           ca1%constant(27), ca1%constant(28),
     .           ca1%constant(29), ca1%constant(30),
     .           ca1%constant(31), ca1%constant(32),
     .           ca1%constant(33), ca1%constant(34),
     .           lmin_coef, taumin, taumax, 
     .           ca1%constant(35), trim( run_file ),
     .           isValid, .false. )

          if ( .not. isValid ) then
            ! Pos. Infinity bit pattern 
            jmatrix(i, :) = 
     .      transfer( int( nanbits ), jmatrix(1,1) )
            ca1%constant(i) = tmp_constant
            cycle
          end if

          ! Set var2 results with results from altering the constants

          call getvariables( var2zt, trim(fname_zt)//".ctl" )
          call getvariables( var2zm, trim(fname_zm)//".ctl" )

          do j = 1, nvarzt
            impact_matrix(i, j) = 
     .      sum(  var2zt%variable(1:nzt,j) - var1zt%variable(1:nzt,j) )
     .      / nzt

            fc_impact_matrix(i, j) = 
     .      sum(  var2zt%variable(1:nzt,j) - var1zt%variable(1:nzt,j) )
     .      / sum( abs( var1zt%variable(1:nzt,j) ) )
          end do

          do j = 1, nvarzt
            jmatrix(i, j) = impact_matrix(i, j) 
     .                    / ( ca1%constant(i) - tmp_constant )
          end do

          do j = 1, nvarzm
            impact_matrix(i, j+nvarzt) 
     .      = sum( var2zm%variable(1:nzm,j) - var1zm%variable(1:nzm,j) )
     .      / nzm

            fc_impact_matrix(i, j+nvarzt) 
     .      = sum( var2zm%variable(1:nzm,j) - var1zm%variable(1:nzm,j) )
     .      / sum( abs( var1zm%variable(1:nzm,j) ) )
          end do

          do j = 1, nvarzm
            jmatrix(i, j+nvarzt) = 
     .      impact_matrix(i,j+nvarzt) / (ca1%constant(i) - tmp_constant)
          end do

        ca1%constant(i) = tmp_constant ! Set parameter back

        end do !i = 1..ca1%entries

       ! Calculate impact matrix
!      do j = 1, var2zt%entries
!        impact_matrix(:,j) = jmatrix(:,j)
!     .                     * ca1%constant(:) * (delta_factor - 1.0)
!      end do 
!
!      do j = 1, var2zm%entries
!        impact_matrix(:,j+nvarzt) = jmatrix(:,j+nvarzt)
!     .                            * ca1%constant(:) 
!     .                            * (delta_factor - 1.0)
!      end do 


        ! Output Results to the terminal
        do i = 1, ca1%entries

            do j = 1, var2zt%entries
              write(unit=*,fmt='(3(a,e10.4))') 
     .        delta//var2zt%name(j)//"/"
     .        //delta//ca1%name(i)//" = ", jmatrix(i, j),
     .        " impact: ", impact_matrix(i, j),
     .        " fc imp: ", fc_impact_matrix(i, j)

            end do

            do j = 1, var2zm%entries
              write(unit=*,fmt='(3(a,e10.4))') 
     .        delta//var2zm%name(j)//"/"
     .        //delta//ca1%name(i)//" = ", jmatrix(i, j+nvarzt),
     .        " impact: ", impact_matrix(i, j+nvarzt),
     .        " fc imp: ", fc_impact_matrix(i, j+nvarzt)

            end do

          write(unit=*,fmt=*) ""

        end do ! 1..ca1%entries

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
        logical error


    
        do k = 1, varray%entries, 1

          varray%variable(1:varray%nz, k) = 
     .    grads_average_interval
     .    ( fname_zx, varray%nz, times(:), varray%name(k), 1, error )

          if ( error ) then
            write(unit=*,fmt=*) "Error in reading"//varray%name(i)
            stop
          end if

        end do

        return
        end subroutine getvariables 

!-----------------------------------------------------------------------
        end program jacobian
