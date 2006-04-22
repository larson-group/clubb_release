!------------------------------------------------------------------------
! $Id: numer_check.F,v 1.2 2006-04-22 20:27:49 dschanen Exp $
        module numerical_check
        implicit none

        contains
!------------------------------------------------------------------------
! Function invalid_model_arrays()
! Checks for invalid floating point values in select model arrays.
!------------------------------------------------------------------------
        logical function invalid_model_arrays( ) 

        use diagnostic_variables
        use prognostic_variables

        implicit none

        invalid_model_arrays = .false.
           
!------------------------------------------------------------------------
!       Check whether any variable array contains a NaN for
!       um, vm, thlm, rtm, rtp2, thlp2, wprtp, wpthlp, rtpthlp, 
!       wp2, & wp3.       
!------------------------------------------------------------------------

        if ( isnan2d( um ) ) then
          print *, "NaN in um model array" 
!         print *, "um= ", um
          invalid_model_arrays = .true.
!         return
        endif

        if ( isnan2d( vm ) ) then
          print *, "NaN in vm model array" 
!         print *, "vm= ", vm
          invalid_model_arrays = .true.
!         return
        endif

        if ( isnan2d( wp2 ) ) then
          print *, "NaN in wp2 model array" 
!         print *, "wp2= ", wp2
          invalid_model_arrays = .true.
!         return
        endif

        if ( isnan2d( wp3 ) ) then
          print *, "NaN in wp3 model array" 
!         print *, "wp3= ", wp3
          invalid_model_arrays = .true.
!         return
        endif

        if ( isnan2d( rtm ) ) then
          print *, "NaN in rtm model array" 
!         print *, "rtm= ", rtm
          invalid_model_arrays = .true.
!         return
        endif

        if ( isnan2d( thlm ) ) then
          print *, "NaN in thlm model array" 
!         print *, "thlm= ", thlm
          invalid_model_arrays = .true.
!         return
        endif

        if ( isnan2d( rtp2 ) ) then
          print *, "NaN in rtp2 model array" 
!         print *, "rtp2= ", rtp2
          invalid_model_arrays = .true.
!         return
        endif

        if ( isnan2d( thlp2 ) ) then
          print *, "NaN in thlp2 model array" 
!         print *, "thlp2= ", thlp2
          invalid_model_arrays = .true.
!         return
        endif

        if ( isnan2d( wprtp ) ) then
          print *, "NaN in wprtp model array" 
!         print *, "wprtp= ", wprtp
          invalid_model_arrays = .true.
!         return
        endif

        if ( isnan2d( wpthlp ) ) then
          print *, "NaN in wpthlp model array" 
!         print *, "wpthlp= ", wpthlp
          invalid_model_arrays = .true.
!         return
        endif

        if ( isnan2d( rtpthlp ) ) then
          print *, "NaN in rtpthlp model array" 
!         print *, "rtpthlp= ", rtpthlp
          invalid_model_arrays = .true.
!         return
        endif

!       if ( isnan2d( Nrm ) ) then
!         print *, "NaN in Nrm model array" 
!         print *, "Nrm= ", Nrm
!         invalid_model_arrays = .true.
!         return
!       endif

!       if ( isnan2d( rrm ) ) then
!         print *, "NaN in rrm model array" 
!         print *, "rrm= ", rrm
!         invalid_model_arrays = .true.
!         return
!       endif

        if ( isnan2d( wmt ) ) then
          print *, "NaN in wmt model array" 
!         print *, "wmt= ", wmt
          invalid_model_arrays = .true.
!         return
        endif

        if ( isnan2d( wp2thvp ) ) then
          print *, "NaN in wp2thvp model array" 
!         print *, "wp2thvp = ", wp2thvp
          invalid_model_arrays = .true.
!         return
        endif

        if ( isnan2d( rtpthvp ) ) then
          print *, "NaN in rtpthvp model array" 
!         print *, "rtpthvp = ", rtpthvp
          invalid_model_arrays = .true.
        endif

        if ( isnan2d( thlpthvp ) ) then
          print *, "NaN in thlpthvp model array" 
!         print *, "thlpthvp = ", thlpthvp
          invalid_model_arrays = .true.
        endif

        return
        end function invalid_model_arrays
!------------------------------------------------------------------------
! function isnan2d( )

! Checks if a given real is NaN, +inf or -inf

! Certain compiler optimizations may cause variables with invalid
! results to flush to zero.  Avoid these!

! Input: x2d: the real type array to be tested.

! Output: isnan2d: true if a variable in tVar is of Not-a-Number type

! Notes: 
! Got a tip from Andy Vaught to use transfer( ) for portability.
! Ideally this would be done with intrinsic module ieee in fortran 2003
!  or with the isnan( ) function provided by some ( not pgf90 ) compilers
!------------------------------------------------------------------------
        logical function isnan2d( x2d )
        
        implicit none

!       External

        intrinsic transfer, any
    
!       Input
        real, intent(in), dimension(:) :: x2d

!       Internal Variables
        integer :: nanbits

        real :: PosInf, NegInf

        data nanbits /Z"7F800000"/

        PosInf = transfer( nanbits, PosInf )
        NegInf = -( transfer( nanbits, NegInf ) )

        isnan2d = .false.

!       This works on compilers with standardized floating point, because
!       the IEEE spec defines that subnormals should not equal themselves.
!       However, all compilers do not seem to follow this. 
        if ( any( ( x2d /= x2d ), dim=1) ) then
          isnan2d = .true. 
          return
        endif

!       This a second check, assuming the above does not work as expected.
        if ( any( ( x2d == PosInf), dim=1) ) then
          isnan2d = .true. 
          return
        endif

!       This may never be needed, it's here just in case.
        if ( any( ( x2d == NegInf), dim=1) ) then
          isnan2d = .true. 
          return
        endif

        return
        end function isnan2d
!------------------------------------------------------------------------
        end module numerical_check
