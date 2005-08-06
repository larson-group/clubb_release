!-----------------------------------------------------------------------
! $Id: jacobian.F,v 1.16 2005-08-06 13:37:53 griffinb Exp $

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
      type avariable
        real, dimension(160) :: value ! this would better allocatable, but
                                      ! fortran is dumb
        character(len=10)    :: name 
      end type avariable
!----------------------------------------------------------------------
      type constant
        real              :: value
        character(len=10) :: name 
      end type constant
!----------------------------------------------------------------------
      type constant_array
        integer entries
        type (constant), allocatable, dimension(:) :: element
      end type constant_array
!----------------------------------------------------------------------
      type variable_array
        integer nz
        integer entries
        type (avariable), allocatable, dimension(:) :: element
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

      type (variable_array)                :: var1
      type (variable_array)                :: var2

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

      allocate( ca1%element( nconst ), 
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
      
      ca1%element(1:nconst)%value =
     .                     (/C1, C2rt, C2thl, C2rtthl, C4, C5, 
     .                       C6rt, C6thl, C7, C8, C10, C11,
     .                       nu1, nu2, nu6, nu8, nu_r, beta,
     .                       gamma_coef, c_K, mu/)     
      ca1%element(1:nconst)%name =
     .  (/"C1        ", "C2rt      ", "C2thl     ", 
     .    "C2rtthl   ", "C4        ", 
     .    "C5        ", "C6rt      ", "C6thl     ", 
     .    "C7        ", 
     .    "C8        ", "C10       ", "C11       ",
     .    "nu1       ", "nu2       ", 
     .    "nu6       ", "nu8       ", "nu_r      ",
     .    "beta      ", "gamma_coef",
     .    "c_K       ", "mu        "/)

      write(*, '(a12,a12,a12)') "constant   ", "initial    ", 
     .                          "varied     "
      do i = 1, ca1%entries
        write(*,*) ca1%element(i)%name, ca1%element(i)%value, 
     .             ca1%element(i)%value * delta_factor
      enddo

      call hoc_model( ca1%element(1)%value, ca1%element(2)%value, 
     .                ca1%element(3)%value, ca1%element(4)%value, 
     .                ca1%element(5)%value, ca1%element(6)%value, 
     .                ca1%element(7)%value, ca1%element(8)%value, 
     .                ca1%element(9)%value, ca1%element(10)%value, 
     .                ca1%element(11)%value, ca1%element(12)%value, 
     .                ca1%element(13)%value, ca1%element(14)%value, 
     .                ca1%element(15)%value, ca1%element(16)%value, 
     .                ca1%element(17)%value, ca1%element(18)%value, 
     .                ca1%element(19)%value, ca1%element(20)%value,
     .                lmin_coef, taumin,
     .                taumax, ca1%element(21)%value,
     .                "../tune/"//trim(run_name)//"_hoc.in", 
     .                isValid, .false. )

      if ( .not. isValid ) stop "Initial run wasn't valid"

!-----------------------------------------------------------------------
      nz = grads_zlvl( trim(run_name)//"_zt.ctl" )

      allocate( var1%element( nvar ),
     .          var2%element( nvar ),
     .          stat = alloc_stat )

      if (alloc_stat /= 0 ) stop "allocate failed"

      var1%entries = nvar
      var1%nz      = nz
      var2%entries = nvar
      var2%nz      = nz



      var1%element(1:nvar)%name = 
     .                     (/"cf        ", "rcm       ", "rtm       ", 
     .                       "thlm      ", "um        ", "vm        ",
     .                       "wp3       ", "wp2zt     ", 
     .                       "rtp2zt    ", "thlp2zt   ", "rtpthlpzt ",
     .                       "wprtpzt   ", "wpthlpzt  "/)

      var2%element(1:nvar)%name = var1%element(1:nvar)%name 

      call getvariables( var1 )

      do i = 1, ca1%entries
        tmp_constant         = ca1%element(i)%value
        ca1%element(i)%value = ca1%element(i)%value * delta_factor

        call hoc_model( ca1%element(1)%value,  ca1%element(2)%value, 
     .                  ca1%element(3)%value,  ca1%element(4)%value, 
     .                  ca1%element(5)%value,  ca1%element(6)%value, 
     .                  ca1%element(7)%value,  ca1%element(8)%value, 
     .                  ca1%element(9)%value,  ca1%element(10)%value, 
     .                  ca1%element(11)%value, ca1%element(12)%value, 
     .                  ca1%element(13)%value, ca1%element(14)%value, 
     .                  ca1%element(15)%value, ca1%element(16)%value, 
     .                  ca1%element(17)%value, ca1%element(18)%value, 
     .                  ca1%element(19)%value, ca1%element(20)%value,
     .                  lmin_coef, taumin,
     .                  taumax, ca1%element(21)%value,
     .                  "../tune/"//trim(run_name)//"_hoc.in", 
     .                  isValid, .false. )

        if ( .not. isValid ) then
          ! Infinity bit pattern 
          jmatrix(i, :) = transfer( int( Z"7F800000" ), real( 1.0 ) )
          ca1%element(i)%value = tmp_constant
          cycle
        endif

        call getvariables( var2 )

        do j = 1, nvar
          jmatrix(i, j) =    sum(   var2%element(j)%value(1:nz)   
     .                            - var1%element(j)%value(1:nz)  )
        enddo

        do j = 1, nvar
          jmatrix(i, j) = jmatrix(i, j) / 
     .                       (tmp_constant - ca1%element(i)%value)
        enddo

        ca1%element(i)%value = tmp_constant
      enddo

!---- Calculate impact matrix ------------------------------------------
      do j = 1, var2%entries
        impact_matrix(:,j) = jmatrix(:,j) *
     .                       ca1%element(:)%value * (delta_factor - 1.)
      enddo 

!---- Output Results to terminal ---------------------------------------
      do i = 1, ca1%entries
          do j = 1, var2%entries
            write(*,*) delta//var2%element(j)%name//" / "
     .                 //delta//ca1%element(i)%name//" =",
     .                 jmatrix(i, j),
     .                 "impact: ",
     .                 impact_matrix(i, j)
          enddo
          write(*,*) ""
      enddo

!---- Output results to an ascii file ----------------------------------
!     Jacobian Matrix
      open(20, file="jacobian_matrix.txt", action="write" )
     
!     Matlab cannot import with this included
!      write(20,'(A)') "% Tuning parameters vary according to columns"
!      write(20,'(A)') "% Output fields vary according to rows"
!      write(20,'(18A)') "% columns:", ca1%element(1:nconst)%name
!      write(20,'(14A)') "% rows:   ", var2%var(1:nvar)%name 
      do j =1, var2%entries
        write(20, '(17e18.10)') jmatrix(:,j)
      enddo
      close(20)

!     Impact Matrix
      open(20, file="impact_matrix.txt", action="write" )

!     Matlab cannot import with this included
!      write(20,'(A)') "% Tuning parameters vary according to columns"
!      write(20,'(A)') "% Output fields vary according to rows"
!      write(20,'(18A)') "% columns:", ca1%element(1:nconst)%name
!      write(20,'(13A)') "% rows:   ", var2%element(1:nvar)%name
      do j =1, var2%entries
        write(20, '(17e18.10)') impact_matrix(:,j)
      enddo

      contains
!-----------------------------------------------------------------------
!      subroutine getvariables( ) 
!      returns a variable_array structure over time intervals
!-----------------------------------------------------------------------
      subroutine getvariables( varray ) 
 
      implicit none

      logical error

      type (variable_array), intent(inout) :: varray
    
      do k = 1, varray%entries
        varray%element(k)%value(1:varray%nz) = 
     .                            grads_average_interval( 
     .                            trim(run_name)//"_zt.ctl", varray%nz,
     .                            times(:), varray%element(k)%name, 
     .                            error )
      enddo

      return
      end subroutine getvariables 

!-----------------------------------------------------------------------
      end program jacobian
