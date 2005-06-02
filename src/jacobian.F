!-----------------------------------------------------------------------
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
      type variable
        real, dimension(150) :: value ! this would better allocatable, but
                                      ! fortran is dumb
        character(len=6)     :: name 
      end type variable
!----------------------------------------------------------------------
      type constant
        real              :: value
        character(len=10) :: name 
      end type constant
!----------------------------------------------------------------------
      type constant_array
        integer entries
        type (constant), allocatable, dimension(:) :: const
      end type constant_array
!----------------------------------------------------------------------
      type variable_array
        integer nz
        integer entries
        type (variable), allocatable, dimension(:) :: var
      end type variable_array

!----- Parameters ------------------------------------------------------
      integer, parameter :: nvar      = 9
      integer, parameter :: nconst    = 17
      real, parameter    :: lmin_coef = 0.5
      real, parameter    :: taumin    = 90.0
      real, parameter    :: taumax    = 3600.0

!      character, parameter :: delta   = greek_'Δ'
      character, parameter :: delta   = 'D'
!----- Model constants -------------------------------------------------
! Notes:  C4 and C10 do not currently play a role in the model 5/31/05
!-----------------------------------------------------------------------
      real       :: C1, C2, C4, C5, C6, C7, C8, C10, C11
      real       :: nu1, nu2, nu6, nu8
      real       :: beta
      real       :: gamma_coef, c_K
      real       :: mu

!----- Model constants multiplier --------------------------------------
      real       :: delta_factor
      real       :: tmp_constant 

      integer, dimension(10)    :: times
      integer                   :: nz
      integer                   :: alloc_stat

      character(len=50)         :: run_name
      logical                   :: use_standard_vars
!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
      type (constant_array)                :: ca1

!      type (constant_array)                :: c_var

      type (variable_array)                :: var1
      type (variable_array)                :: var2

      real, dimension(nconst, nvar)        :: jmatrix
      real, dimension(nconst, nvar)        :: impact_matrix

      integer i, j, k        ! loop variables
      logical isValid

!----- Namelists -------------------------------------------------------
      namelist /model/    run_name, times, delta_factor, 
     .                    use_standard_vars

      namelist /initvars/ C1, C2, C4, C5, C6, C7, C8, C10, C11, 
     .                    nu1, nu2, nu6, nu8, beta, gamma_coef,
     .                    c_K, mu

      times      = 0
      isValid    = .true.

      allocate( ca1%const( nconst ), 
     .          stat = alloc_stat )
      ca1%entries = nconst

      if (alloc_stat /= 0 ) stop "allocate failed"


      open(10, file = 'jacobian.in', status='old')
      read(10, nml = model)
      if ( .not. use_standard_vars ) then
        read(10, nml = initvars)
      else
!       Original values (hoc_v2, v2.1)
        C1  = 1.7
        C2  = 0.8 * 1.3
        C4  = 4.5
        C5  = 0.25
        C6  = 4.85 * 2.0
        C7  = 0.8
        C8  = 0.7 * 0.6 * 6.5
        C10 = 3. * C2
        C11 = 0.2

        nu1 = 20.
        nu2 = 5.
        nu6 = 5.
        nu8 = 20.

        gamma_coef = 0.45
        c_K        = 0.548
        mu         = 6.e-4
      endif
      close(10)
      
      ca1%const(1:nconst)%value =
     .                     (/C1, C2, C4, C5, C6, C7, C8, C10, C11,
     .                       nu1, nu2, nu6, nu8, beta, gamma_coef,
     .                       c_K, mu/)     
      ca1%const(1:nconst)%name =
     .                    (/"C1         ", "C2", "C4", "C5", "C6", 
     .                      "C7", "C8", "C10", "C11",
     .                      "nu1", "nu2", "nu6", "nu8", "beta", 
     .                      "gamma_coef",
     .                      "c_K", "mu"/)

      write(*, 'a12,a12,a12'), "constant", "initial", "varied"
      do i = 1, ca1%entries
        write(*,*) ca1%const(i)%name, ca1%const(i)%value, 
     .             ca1%const(i)%value * delta_factor
      enddo

      call hoc_sub( ca1%const(1)%value, ca1%const(2)%value, 
     .              ca1%const(3)%value, ca1%const(4)%value, 
     .              ca1%const(5)%value, ca1%const(6)%value, 
     .              ca1%const(7)%value, ca1%const(8)%value, 
     .              ca1%const(9)%value, ca1%const(10)%value, 
     .              ca1%const(11)%value, ca1%const(12)%value, 
     .              ca1%const(13)%value, ca1%const(14)%value, 
     .              ca1%const(15)%value, ca1%const(16)%value, 
     .              lmin_coef, taumin,
     .              taumax, ca1%const(17)%value,
     .                "../tune/"//trim(run_name)//"_hoc.in", 
     .              trim(run_name), 
     .              isValid, .false. )

      if ( .not. isValid ) stop "Initial run wasn't valid"

!-----------------------------------------------------------------------
      nz = grads_zlvl( trim(run_name)//"_zt.ctl" )

      allocate( var1%var( nvar ),
     .          var2%var( nvar ),
     .          stat = alloc_stat )

      if (alloc_stat /= 0 ) stop "allocate failed"

      var1%entries = nvar
      var1%nz      = nz
      var2%entries = nvar
      var2%nz      = nz



      var1%var(1:nvar)%name = 
     .                     (/"rcm   ", "rtm   ", "cf    ", "thlm  ",
     .                       "um    ", "vm    ", "wp3   ", "wp2rtp", 
     .                       "wp2rcp"/)

      var2%var(1:nvar)%name = var1%var(1:nvar)%name 

      call getvariables( var1 )

      do i = 1, ca1%entries
!        print *, i
        tmp_constant       = ca1%const(i)%value
        ca1%const(i)%value = ca1%const(i)%value * delta_factor

        call hoc_sub( ca1%const(1)%value, ca1%const(2)%value, 
     .                ca1%const(3)%value, ca1%const(4)%value, 
     .                ca1%const(5)%value, ca1%const(6)%value, 
     .                ca1%const(7)%value, ca1%const(8)%value, 
     .                ca1%const(9)%value, ca1%const(10)%value, 
     .                ca1%const(11)%value, ca1%const(12)%value, 
     .                ca1%const(13)%value, ca1%const(14)%value, 
     .                ca1%const(15)%value, ca1%const(16)%value, 
     .                lmin_coef, taumin,
     .                taumax, ca1%const(17)%value,
     .                "../tune/"//trim(run_name)//"_hoc.in", 
     .                trim(run_name), 
     .                isValid, .false. )

        if ( .not. isValid ) then
          jmatrix(i, :) = transfer(Z"7F800000", 1.0)
          ca1%const(i)%value = tmp_constant
          cycle
        endif

        call getvariables( var2 )

        do j = 1, nvar
          jmatrix(i, j) =    sum( abs(   var1%var(j)%value(1:nz)   
     .                                 - var2%var(j)%value(1:nz)  ) )
        enddo

        do j = 1, nvar
          jmatrix(i, j) = jmatrix(i, j) / 
     .                    abs( (ca1%const(i)%value - tmp_constant) )
        enddo

        ca1%const(i)%value = tmp_constant
      enddo

!---- Calculate impact matrix ----
      do j = 1, var2%entries
        impact_matrix(:,j) = jmatrix(:,j) *
     .                       ( ca1%const(:)%value * delta_factor )
      enddo 

!---- Output Results to terminal ----
      do i = 1, ca1%entries
          do j = 1, var2%entries
            write(*,*) delta//var2%var(j)%name//" / "
     .                 //delta//ca1%const(i)%name,
     .                 jmatrix(i, j),
     .                 "impact: ",
     .                 impact_matrix(i, j)
          enddo
          write(*,*) ""
      enddo

!---- Output results to an ascii file ----
      open(20, file="results_jmatrix.txt", action="write" )
      do j =1, var2%entries
        write(20, '(17e)') jmatrix(:,j)
      enddo
      close(20)

      open(20, file="results_imatrix.txt", action="write" )
      do j =1, var2%entries
        write(20, '(17e)') impact_matrix(:,j)
      enddo

      contains
!-----------------------------------------------------------------------
!      subroutine getvariables( ) 
!      returns an array of variables over times
!-----------------------------------------------------------------------
      subroutine getvariables( varray ) 
      use grads_common
 
      implicit none

      logical error

      type (variable_array), intent(inout) :: varray
    
      do k = 1, varray%entries
      varray%var(k)%value(:) =  grads_average_interval( 
     .                          trim(run_name)//"_zt.ctl", varray%nz,
     .                          times(:), varray%var(k)%name, 
     .                          error )
      enddo

      return
      end subroutine getvariables 

      end program jacobian
