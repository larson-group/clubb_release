!-----------------------------------------------------------------------
! $Id: jacobian.F,v 1.17 2005-08-31 22:50:48 dschanen Exp $

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

        character(len=10), pointer :: name(:)

        real, pointer :: variable(:,:) ! (1:nz, entries)

      end type variable_array

!----- Parameters ------------------------------------------------------
      integer, parameter :: nvar      = 13
      integer, parameter :: nconst    = 21

!      character, parameter :: delta   = greek_'Δ'
      character, parameter :: delta   = 'D'

!----- Model constants -------------------------------------------------
! Notes:  C4 and C10 do not currently play a role in the model 5/31/05
!-----------------------------------------------------------------------
!      real       :: C1, C2, C4, C5, C6, C7, C8, C10, C11
      real       :: C1, C2rt, C2thl, C2rtthl, C4, C5, C6rt, C6thl
      real       :: C7, C8, C10, C11
      real       :: nu1, nu2, nu6, nu8, nu_r
      real       :: beta
      real       :: gamma_coef, c_K
      real       :: mu

!     These aren't varied in matrix generation
      real       :: lmin_coef  
      real       :: taumin    
      real       :: taumax

      real       :: tmp_constant 

!----- Namelist Variables ----------------------------------------------
      real                      :: delta_factor
      integer, dimension(10)    :: times
      character(len=50)         :: run_name
      logical                   :: use_standard_vars

      integer                   :: nz
      integer                   :: alloc_stat

!----------------------------------------------------------------------
!  Variable and constant data structures
!----------------------------------------------------------------------
      type (constant_array)                :: ca1

      type (variable_array)                :: var1zt
      type (variable_array)                :: var2zt

      real, dimension(nconst, nvar)        :: jmatrix
      real, dimension(nconst, nvar)        :: impact_matrix

      integer i, j, k        ! loop variables
      logical isValid

!----- Namelists -------------------------------------------------------
      namelist /model/    run_name, times, delta_factor, 
     .                    use_standard_vars

      namelist /initvars/ C1, C2rt, C2thl, C2rtthl, C4, C5, 
     .                    C6rt, C6thl, C7, C8, C10, C11, 
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
      
      ca1%constant(1:nconst) = (/C1, C2rt, C2thl, C2rtthl, C4, C5, 
     .                           C6rt, C6thl, C7, C8, C10, C11,
     .                           nu1, nu2, nu6, nu8, nu_r, beta,
     .                           gamma_coef, c_K, mu/)     
      ca1%name(1:nconst) =
     .  (/"C1        ", "C2rt      ", "C2thl     ", 
     .    "C2rtthl   ", "C4        ", 
     .    "C5        ", "C6rt      ", "C6thl     ", 
     .    "C7        ", 
     .    "C8        ", "C10       ", "C11       ",
     .    "nu1       ", "nu2       ", 
     .    "nu6       ", "nu8       ", "nu_r      ",
     .    "beta      ", "gamma_coef",
     .    "c_K       ", "mu        "/)

      write(*, '(3a12)') "constant   ", "initial    ", "varied     "

      do i = 1, ca1%entries
        write(*,*) ca1%name(i), ca1%constant(i), 
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
     .                lmin_coef, taumin,
     .                taumax, ca1%constant(21),
     .                "../tune/"//trim(run_name)//"_hoc.in", 
     .                isValid, .false. )

      if ( .not. isValid ) stop "Initial run wasn't valid"

!-----------------------------------------------------------------------
      nz = grads_zlvl( trim(run_name)//"_zt.ctl" )

      allocate( var1zt%variable( nz, nvar ),
     .          var2zt%variable( nz, nvar ),
     .          var1zt%name( nvar ),
     .          var2zt%name( nvar ),
     .          stat = alloc_stat )

      if (alloc_stat /= 0 ) stop "allocate failed"

      var1zt%entries = nvar
      var1zt%nz      = nz
      var2zt%entries = nvar
      var2zt%nz      = nz



      var1zt%name(1:nvar) = 
     .                  (/"cf        ", "rcm       ", "rtm       ", 
     .                    "thlm      ", "um        ", "vm        ",
     .                    "wp3       ", "wp2zt     ", 
     .                    "rtp2zt    ", "thlp2zt   ", "rtpthlpzt ",
     .                    "wprtpzt   ", "wpthlpzt  "/)

      var2zt%name(1:nvar) = var1zt%name(1:nvar) 

      call getvariables( var1zt, "zt" )

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
     .                  lmin_coef, taumin,
     .                  taumax, ca1%constant(21),
     .                  "../tune/"//trim(run_name)//"_hoc.in", 
     .                  isValid, .false. )

        if ( .not. isValid ) then
          ! Infinity bit pattern 
          jmatrix(i, :) = transfer( int( Z"7F800000" ), 
     .                              jmatrix(1,1) )
          ca1%constant(i) = tmp_constant
          cycle
        endif

        call getvariables( var2zt, "zt" )

        do j = 1, nvar
          jmatrix(i, j) = sum(  var2zt%variable(1:nz,j)   
     .                        - var1zt%variable(1:nz,j) )
        enddo

        do j = 1, nvar
          jmatrix(i, j) = jmatrix(i, j) / 
     .                       (tmp_constant - ca1%constant(i) )
        enddo

        ca1%constant(i) = tmp_constant
      enddo

!---- Calculate impact matrix ------------------------------------------
      do j = 1, var2zt%entries
        impact_matrix(:,j) = jmatrix(:,j) *
     .                       ca1%constant(:) * (delta_factor - 1.0)
      enddo 

!---- Output Results to terminal ---------------------------------------
      do i = 1, ca1%entries
          do j = 1, var2zt%entries
            write(*,*) delta//var2zt%name(j)//" / "
     .                 //delta//ca1%name(i)//" =",
     .                 jmatrix(i, j),
     .                 "impact: ", impact_matrix(i, j)
          enddo
          write(*,*) ""
      enddo

!---- Output results to an ascii file ----------------------------------
!     Jacobian Matrix
      open(20, file="jacobian_matrix.txt", action="write" )
     
!     Matlab cannot import with this included
!      write(20,'(A)') "% Tuning parameters vary according to columns"
!      write(20,'(A)') "% Output fields vary according to rows"
!      write(20,'(18A)') "% columns:", ca1%name(1:nconst)
!      write(20,'(14A)') "% rows:   ", var2%name(1:nvar)
      do j =1, var2zt%entries
        write(20, '(17e18.10)') jmatrix(:,j)
      enddo
      close(20)

!     Impact Matrix
      open(20, file="impact_matrix.txt", action="write" )

!     Matlab cannot import with this included
!      write(20,'(A)') "% Tuning parameters vary according to columns"
!      write(20,'(A)') "% Output fields vary according to rows"
!      write(20,'(18A)') "% columns:", ca1%name(1:nconst)
!      write(20,'(13A)') "% rows:   ", var2%name(1:nvar)
      do j =1, var2zt%entries
        write(20, '(17e18.10)') impact_matrix(:,j)
      enddo

      contains
!-----------------------------------------------------------------------
!      subroutine getvariables( ) 
!      returns a variable_array structure over time intervals
!-----------------------------------------------------------------------
      subroutine getvariables( varray, grid ) 
 
      implicit none

      logical error

!     Input

      type (variable_array), intent(inout) :: varray

      character(len=2) grid

    
      do k = 1, varray%entries
        varray%variable(1:varray%nz, k) = 
     .    grads_average_interval( trim(run_name)//"_"//grid//".ctl",
     .                            varray%nz,
     .                            times(:), varray%name(k), error )
        if ( error ) stop "getvariables(): error"
      enddo

      return
      end subroutine getvariables 

!-----------------------------------------------------------------------
      end program jacobian
