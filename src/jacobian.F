!-----------------------------------------------------------------------
! $Id: jacobian.F,v 1.24 2005-11-05 18:23:25 griffinb Exp $

! Program Jacobian

! Generates a matrix based on variation between constants (C1,C2...etc)
! and variables (cf,rcm,thlm...etc)

!-----------------------------------------------------------------------

      program jacobian

      use hoc
      use grads_common

      implicit none

!-----------------------------------------------------------------------
!  Variable structures
!-----------------------------------------------------------------------
      type constant_array

        integer entries

        character(len=10), pointer :: name(:)

        real, pointer :: constant(:)

      end type constant_array
!----------------------------------------------------------------------
      type variable_array

        integer nz
        integer entries

        character(len=12), pointer :: name(:)

        real, pointer :: variable(:,:) ! (1:nz, entries)

      end type variable_array

!----- Parameters ------------------------------------------------------
      integer, parameter :: nvarzt = 14
      integer, parameter :: nvarzm = 40
      integer, parameter :: nconst = 31

!     31 must be changed to be equal to nconst

      character(len=12), parameter :: write_format = "(27(e18.10))"

!     character, parameter :: delta   = greek_'Δ' ! Doesn't work
      character, parameter :: delta   = 'D'

!----- Model constants -------------------------------------------------
! Notes:  C4 and C10 do not currently play a role in the model 5/31/05
!-----------------------------------------------------------------------
      real       :: C1, C1b, C1c, C2rt, C2thl, C2rtthl
      real       :: C4, C5, C6rt, C6rtb, C6rtc, C6thl
      real       :: C7, C7b, C7c, C8, C8b, C10, C11, C11b, C11c, C12
      real       :: nu1, nu2, nu6, nu8, nu_r
      real       :: beta
      real       :: gamma_coef, c_K
      real       :: mu

!     These aren't varied in matrix generation
      real       :: lmin_coef  
      real       :: taumin    
      real       :: taumax

      real       :: tmp_constant 

!     Namelist Variables
      real                      :: delta_factor
      integer, dimension(10)    :: times
      character(len=50)         :: run_name
      logical                   :: use_standard_vars

!     Internal
      integer                   :: nzt
      integer                   :: nzm
      integer                   :: alloc_stat

!----------------------------------------------------------------------
!  Variable and constant data structures
!----------------------------------------------------------------------
      type (constant_array) :: ca1

      type (variable_array) :: var1zt
      type (variable_array) :: var2zt

      type (variable_array) :: var1zm
      type (variable_array) :: var2zm

      real, dimension(nconst, nvarzt+nvarzm) :: jmatrix
      real, dimension(nconst, nvarzt+nvarzm) :: impact_matrix
      real, dimension(nconst, nvarzt+nvarzm) :: fc_impact_matrix

      integer i, j, k        ! loop variables
      logical isValid

!----- Namelists -------------------------------------------------------
      namelist /model/    run_name, times, delta_factor, 
     .                    use_standard_vars

      namelist /initvars/ C1, C1b, C1c, C2rt, C2thl, C2rtthl,
     .                    C4, C5, C6rt, C6rtb, C6rtc, C6thl,
     .                    C7, C7b, C7c, C8, C8b, C10, 
     .                    C11, C11b, C11c, C12,
     .                    nu1, nu2, nu6, nu8, nu_r, beta, gamma_coef,
     .                    c_K, lmin_coef, taumin, taumax, mu

      times      = 0
      isValid    = .true.

      allocate( ca1%constant( nconst ), ca1%name( nconst ),
     .          stat = alloc_stat )
      ca1%entries = nconst

      if (alloc_stat /= 0 ) stop "allocate failed"


      open(10, file = 'jacobian.in', status='old')
      read(10, nml = model)
      if ( .not. use_standard_vars ) then
        read(10, nml = initvars)
      else
        open(20, file = "../standalone/standalone_"
     .                //trim(run_name)//".in", status="old")

        read(20, nml = initvars)
        close(20)
      endif
      close(10)
      
      ca1%constant(1:nconst) = (/C1, C1b, C1c, C2rt, C2thl, C2rtthl,
     .                           C4, C5, C6rt, C6rtb, C6rtc, C6thl,
     .                           C7, C7b, C7c, C8, C8b,
     .                           C10, C11, C11b, C11c, C12,
     .                           nu1, nu2, nu6, nu8, nu_r, beta,
     .                           gamma_coef, c_K, mu/)     
      ca1%name(1:nconst) =
     .  (/"C1        ", "C1b       ", "C1c       ",
     .    "C2rt      ", "C2thl     ", "C2rtthl   ",
     .    "C4        ", "C5        ",
     .    "C6rt      ", "C6rtb     ", "C6rtc     ",
     .    "C6thl     ", 
     .    "C7        ", "C7b       ", "C7c       ",
     .    "C8        ", "C8b       ", 
     .    "C10       ", "C11       ", "C11b      ",
     .    "C11c      ", "C12       ",
     .    "nu1       ", "nu2       ", 
     .    "nu6       ", "nu8       ", "nu_r      ",
     .    "beta      ", "gamma_coef",
     .    "c_K       ", "mu        "/)

      write(*, '(3a12)') "Parameter", "Initial", "Varied"

      do i = 1, ca1%entries
        write(*,'(a12,2f12.5)') ca1%name(i), ca1%constant(i), 
     .             ca1%constant(i) * delta_factor
      enddo

      call hoc_model( ca1%constant(1), ca1%constant(2), 
     .                ca1%constant(3), ca1%constant(4), 
     .                ca1%constant(5), ca1%constant(6), 
     .                ca1%constant(7), ca1%constant(8), 
     .                ca1%constant(9), ca1%constant(10), 
     .                ca1%constant(11), ca1%constant(12), 
     .                ca1%constant(13), ca1%constant(14), 
     .                ca1%constant(15), ca1%constant(16), 
     .                ca1%constant(17), ca1%constant(18), 
     .                ca1%constant(19), ca1%constant(20),
     .                ca1%constant(21), ca1%constant(22),
     .                ca1%constant(23), ca1%constant(24),
     .                ca1%constant(25), ca1%constant(26),
     .                ca1%constant(27), ca1%constant(28),
     .                ca1%constant(29), ca1%constant(30),
     .                lmin_coef, taumin, taumax, 
     .                ca1%constant(31),
     .                "../tune/"//trim(run_name)//"_hoc.in", 
     .                isValid, .false. )

      if ( .not. isValid ) stop "Initial run wasn't valid"

!     Obtain nz from the generated GrADS files

      nzt = grads_zlvl( trim(run_name)//"_zt.ctl" )
      nzm = grads_zlvl( trim(run_name)//"_zm.ctl" )


!     Initialize the structures holding the variables

      allocate( var1zt%variable( nzt, nvarzt ),
     .          var2zt%variable( nzt, nvarzt ),
     .          var1zt%name( nvarzt ),
     .          var2zt%name( nvarzt ),
     .          stat = alloc_stat )

      if (alloc_stat /= 0 ) stop "allocate failed"

      var1zt%entries = nvarzt 
      var1zt%nz      = nzt
      var2zt%entries = nvarzt 
      var2zt%nz      = nzt

      allocate( var1zm%variable( nzm, nvarzm ),
     .          var2zm%variable( nzm, nvarzm ),
     .          var1zm%name( nvarzm ),
     .          var2zm%name( nvarzm ),
     .          stat = alloc_stat )

      if (alloc_stat /= 0 ) stop "allocate failed"

      var1zm%entries = nvarzm 
      var1zm%nz      = nzm
      var2zm%entries = nvarzm 
      var2zm%nz      = nzm



      var1zt%name(1:nvarzt) = 
     .           (/"cf          ", "rcm         ", "rtm         ", 
     .             "thlm        ", "um          ", "vm          ",
     .             "wp3         ", "wp3_ta      ", "wp3_tp      ",
     .             "wp3_bp      ", "wp3_cl      ",
     .             "wp3_pr1     ", "wp3_pr2     ", "wp3_dp1     "/)

      var1zm%name(1:nvarzm) = 
     .           (/"wp2         ", "rtp2        ", "thlp2       ", 
     .             "rtpthlp     ", "wprtp       ", "wpthlp      ",
     .             "wp2_ta      ", "wp2_bp      ", "wp2_pr2     ",
     .             "wp2_pr3     ", "wp2_dp1     ", "wp2_dp2     ",
     .             "wp2_cl      ", 
     .             "wprtp_ta    ", "wprtp_tp    ", "wprtp_bp    ", 
     .             "wprtp_pr1   ", "wprtp_pr2   ", 
     .             "wprtp_pr3   ", "wprtp_dp1   ",
     .             "wpthlp_ta   ", "wpthlp_tp   ", "wpthlp_bp   ", 
     .             "wpthlp_pr1  ", "wpthlp_pr2  ", 
     .             "wpthlp_pr3  ", "wpthlp_dp1  ",
     .             "rtp2_ta     ", "rtp2_tp     ", 
     .             "rtp2_dp1    ", "rtp2_dp2    ", 
     .             "thlp2_ta    ", "thlp2_tp    ", 
     .             "thlp2_dp1   ", "thlp2_dp2   ",
     .             "rtpthlp_ta  ", "rtpthlp_tp1 ", "rtpthlp_tp2 ", 
     .             "rtpthlp_dp1 ", "rtpthlp_dp2 "/)

      var2zt%name(1:nvarzt) = var1zt%name(1:nvarzt) 

      var2zm%name(1:nvarzm) = var1zm%name(1:nvarzm) 

!     Set var1 fields with initial run results

      call getvariables( var1zt, "zt" )

      call getvariables( var1zm, "zm" )

      do i = 1, ca1%entries
        tmp_constant    = ca1%constant(i)
        ca1%constant(i) = ca1%constant(i) * delta_factor

        call hoc_model( ca1%constant(1),  ca1%constant(2), 
     .                  ca1%constant(3),  ca1%constant(4), 
     .                  ca1%constant(5),  ca1%constant(6), 
     .                  ca1%constant(7),  ca1%constant(8), 
     .                  ca1%constant(9),  ca1%constant(10), 
     .                  ca1%constant(11), ca1%constant(12), 
     .                  ca1%constant(13), ca1%constant(14), 
     .                  ca1%constant(15), ca1%constant(16), 
     .                  ca1%constant(17), ca1%constant(18), 
     .                  ca1%constant(19), ca1%constant(20),
     .                  ca1%constant(21), ca1%constant(22),
     .                  ca1%constant(23), ca1%constant(24),
     .                  ca1%constant(25), ca1%constant(26),
     .                  ca1%constant(27), ca1%constant(28),
     .                  ca1%constant(29), ca1%constant(30),
     .                  lmin_coef, taumin, taumax, 
     .                  ca1%constant(31),
     .                  "../tune/"//trim(run_name)//"_hoc.in", 
     .                  isValid, .false. )

        if ( .not. isValid ) then
          ! Infinity bit pattern 
          jmatrix(i, :)   = transfer( int( Z"7F800000" ), 
     .                              jmatrix(1,1) )
          ca1%constant(i) = tmp_constant
          cycle
        endif

!       Set var2 results with results from altering the constants

        call getvariables( var2zt, "zt" )
        call getvariables( var2zm, "zm" )

        do j = 1, nvarzt
          impact_matrix(i, j) = sum(  var2zt%variable(1:nzt,j)   
     .                              - var1zt%variable(1:nzt,j) )
     .                        / nzt
          fc_impact_matrix(i, j) = sum(  var2zt%variable(1:nzt,j)   
     .                                 - var1zt%variable(1:nzt,j) )
     .                        / sum( abs( var1zt%variable(1:nzt,j) ) )
        enddo

        do j = 1, nvarzt
          jmatrix(i, j) = impact_matrix(i, j) 
     .                  / ( ca1%constant(i) - tmp_constant )
        enddo

        do j = 1, nvarzm
          impact_matrix(i, j+nvarzt) 
     .    = sum(  var2zm%variable(1:nzm,j)   
     .          - var1zm%variable(1:nzm,j) )
     .    / nzm
          fc_impact_matrix(i, j+nvarzt) 
     .    = sum(  var2zm%variable(1:nzm,j)   
     .          - var1zm%variable(1:nzm,j) )
     .    / sum( abs( var1zm%variable(1:nzm,j) ) )
        enddo

        do j = 1, nvarzm
          jmatrix(i, j+nvarzt) = impact_matrix(i, j+nvarzt)
     .                         / ( ca1%constant(i) - tmp_constant )
        enddo

        ca1%constant(i) = tmp_constant
      enddo

!---- Calculate impact matrix ------------------------------------------
!      do j = 1, var2zt%entries
!        impact_matrix(:,j) = jmatrix(:,j)
!     .                     * ca1%constant(:) * (delta_factor - 1.0)
!      enddo 
!
!      do j = 1, var2zm%entries
!        impact_matrix(:,j+nvarzt) = jmatrix(:,j+nvarzt)
!     .                            * ca1%constant(:) 
!     .                            * (delta_factor - 1.0)
!      enddo 


!---- Output Results to terminal ---------------------------------------
      do i = 1, ca1%entries

          do j = 1, var2zt%entries
            write(*,'(3(a,e10.4))') 
     .      delta//var2zt%name(j)//"/"
     .      //delta//ca1%name(i)//" = ", jmatrix(i, j),
     .      " impact: ", impact_matrix(i, j),
     .      " fc imp: ", fc_impact_matrix(i, j)
          enddo

          do j = 1, var2zm%entries
            write(*,'(3(a,e10.4))') 
     .      delta//var2zm%name(j)//"/"
     .      //delta//ca1%name(i)//" = ", jmatrix(i, j+nvarzt),
     .      " impact: ", impact_matrix(i, j+nvarzt),
     .      " fc imp: ", fc_impact_matrix(i, j+nvarzt)

          enddo

          write(*,*) ""

      enddo

!---- Output results to an ASCII file ----------------------------------

!     Jacobian Matrix
      open(20, file="jacobian_matrix.txt", action="write" )

      do j = 1, nvarzm + nvarzt
        write(20, write_format) jmatrix(:,j)
      enddo
      close(20)

!     Impact Matrix
      open(20, file="impact_matrix.txt", action="write" )

      do j =1, nvarzm + nvarzt
        write(20, write_format) impact_matrix(:,j)
      enddo

      contains
!-----------------------------------------------------------------------
!     SUBROUTINE getvariables( ) 
!     returns a variable_array structure over namelist defined
!     time intervals
!-----------------------------------------------------------------------
      subroutine getvariables( varray, grid ) 
 
      implicit none

      logical error

!     Input

      type (variable_array), intent(inout) :: varray
      character(len=2), intent(in)         :: grid

    
      do k = 1, varray%entries
        varray%variable(1:varray%nz, k) = 
     .    grads_average_interval( trim(run_name)//"_"//grid//".ctl",
     .                            varray%nz,
     .                            times(:), varray%name(k), 1, error )
        if ( error ) stop "getvariables(): error"
      enddo

      return
      end subroutine getvariables 

!-----------------------------------------------------------------------
      end program jacobian
