!----------------------------------------------------------------------
!  COMPARE_RUNS.F
!    Tests GRaDS file variables for low variation
!----------------------------------------------------------------------
        program compare_runs 

        use inputfile_class
        use grads_common        

        implicit none
 
!        integer, parameter               :: nvar = 14
        integer, parameter               :: nvar = 9
        integer, parameter               :: max_times = 10


        real(kind=4), allocatable, dimension(:)      :: fvar1
        real(kind=4), allocatable, dimension(:)      :: fvar2
        real(kind=4), allocatable, dimension(:)      :: fvar3

        character(len=10), dimension(nvar)           :: hoc_variables
        character(len=10), dimension(nvar)           :: les_variables

        real(kind=8) :: percent_mean_sqr
        real(kind=8) :: reference

        real(kind=4) :: minmax

        character(len=80)     :: file1 ! Always a HOC file
        character(len=80)     :: file2 ! Either a HOC or LES file
        character(len=80)     :: file3 ! Always the 2.1 reference files

        integer, dimension(max_times) :: t        ! times to analyze

        integer                ::       nz, nz_file1,
     .                            nz_file2, nz_file3    ! z-levels
        integer                :: n               ! loop
        integer                :: i               ! loop
        integer                :: tmax            ! max time to analyze 

        logical                :: les_comparison

        logical                :: error           ! set, never used

        namelist /compare/ les_comparison, file1, file2, file3, t

        open ( 10, file="compare_runs.in", status='old' )
        read ( 10, nml = compare )
        
        hoc_variables(1)  = "rcm"
        hoc_variables(2)  = "rtm" 
        hoc_variables(3)  = "cf" 
        hoc_variables(4)  = "thlm"
        hoc_variables(5)  = "um"
        hoc_variables(6)  = "vm"
        hoc_variables(7)  = "wp3"
        hoc_variables(8)  = "wp2rtp"
        hoc_variables(9)  = "wp2rcp"
!        hoc_variables(10) = "wp2zt"
!        hoc_variables(11) = "thlp2zt"
!        hoc_variables(12) = "wpthlpzt"
!        hoc_variables(13) = "wprtpzt"
!        hoc_variables(14) = "rtp2zt"

        if ( les_comparison ) then
          les_variables(1)  = "qcm"
          les_variables(2)  = "qtm" 
          les_variables(3)  = "cf" 
          les_variables(4)  = "thlm"
          les_variables(5)  = "um"
          les_variables(6)  = "vm"
          les_variables(7)  = "wp3"
          les_variables(8)  = "wp2qtp"
          les_variables(9)  = "wp2qlp"
!          les_variables(10) = "wp2"
!          les_variables(11) = "thlp2"
!          les_variables(12) = "wpthlp"
!          les_variables(13) = "wpqtp"
!          les_variables(14) = "qtp2"
        else
          les_variables( 1:nvar ) = hoc_variables( 1:nvar )
        endif

        do i=1, max_times
          if ( t(i) == 0 ) exit
        enddo
        tmax = i - 1
 
        do i=1, tmax, 2
        write(*,*) "time: ", t(i), " to ", t(i+1)
        write(*,*) "Comparison","                (file1 - file2)^2",
     .             "          (file2 - file3)^2"
          do n=1, nvar

!           hoc_v2.1 reference file
            nz_file3 = grads_zlvl( file3 )
            allocate( fvar3(nz_file3) )

            fvar3 = grads_average( file3, nz_file3, t(i), t(i+1), 
     .                             hoc_variables(n), error )

!           second hoc file or the les file
            nz_file2 = grads_zlvl( file2 )
            allocate( fvar2(nz_file2) )

            fvar2 = grads_average( file2, nz_file2, t(i), t(i+1), 
     .                             les_variables(n), error )

!           first hoc file
            nz_file1 = grads_zlvl( file1 )
            allocate( fvar1(nz_file1) )

            fvar1 = grads_average( file1, nz_file1, t(i), t(i+1), 
     .                             hoc_variables(n), error )

!           Choose smallest number of vertical gridpoints.
            nz = min( nz_file3, nz_file2, nz_file1 )

            minmax = maxval( fvar2 ) - minval( fvar2 )
            if ( minmax == 0. ) then
              write(*,*) "minmax = 0. for ", les_variables(n)
              minmax = 1.             ! no normalization
            endif

            if ( .not. les_comparison ) then
            percent_mean_sqr = sum( ( ( ( fvar1(1:nz)-fvar2(1:nz) )
     .                                  / minmax )        )**2, 1 )
     .                                                 / (nz)
            else
              call les_interpolate( )
            endif

            reference = sum( ( ( ( fvar2(1:nz)-fvar3(1:nz) )
     .                                  / minmax )        )**2, 1 )
     .                                                 / (nz)

            write(*,*) "PM^2 for "//hoc_variables(n)// "=", 
     .                  percent_mean_sqr, reference
            deallocate ( fvar1, fvar2 )
          enddo
        enddo

        contains

        subroutine les_interpolate( ) 
        
        implicit none

        select case ( nz )
          case (48, 75, 110, 147, 150)
!----------------------------------------------------------------------
!  Due to hoc's lower starting point, we can only use
!  (total number of z-levels) - 1 (a maximum of 74 for BOMEX).
!  (Need to interpolate hoc to LES' levels.  Right now we just
!   compare adjacent z levels.  Vince Larson 12 Jan 2005)
!----------------------------------------------------------------------
             percent_mean_sqr = sum(((( fvar1(2:nz)-fvar2(1:(nz-1))) 
     .                              / minmax )         ) ** 2, 1 )
     .                                               / (nz-1)
          case (132)  !  the dycoms II RF01 case
             percent_mean_sqr = sum(((( fvar1(3:nz) - fvar2(1:(nz-2))) 
     .                               / minmax )         ) ** 2, 1 )
     .                                               / (nz-2)
          case (50)   !  the wangara case
             percent_mean_sqr = sum(((( fvar1(1:nz)-fvar2(1:nz) )
     .                               / minmax )        )**2, 1 )
     .                                                 / (nz)
          case default !
            stop "Not able to handle specified number of HOC z-levels"
        end select

        return
        end subroutine les_interpolate

        end program compare_runs 
!----------------------------------------------------------------------
