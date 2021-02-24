! $Id$
!==========================================
! 
! Description:
!   This file tests (mainly) the ESA minimization algorithm that the clubb tuner
!   uses. It is tested with a number of different function and the passing criteria
!   calls for the functions in question to locate the global minimum a certain
!   percent of the time. This percentage was determined experimentally by spending
!   a long time seeing what the best percentage I could get was. The variables that
!   were changed to increase this global minimum finding ratio are intial_temp, 
!   stp_adjst_center, and stp_adjust_spread
! 
! Functions:
! 
!   Goldstein-Price Function (2 vars) - This is the simplest one with 1 global minimum and 
!       few tricky local minimum, this was mostly used for a quick test. The new esa_driver 
!       is capable of finding the global minima over 99% of the time, but the criteria is 
!       90% to avoid false fails.
! 
!   Himmelblau Function (2 vars) - This one has multiple (possible insignificant) local minima 
!       and 4 significantly different global minima. This test is useful to ensure that the 
!       algorithm is capable of finding all the global minima at different times, this gives us 
!       confidence that the algorithm is random enough to yield different results while accurate 
!       enough to still locate the global minimums. The algorithm finds any one of the global 
!       minimums usually 100% of the time, set to 95% to be safe.
! 
!   Schaffer Function #2 (2 vars) - This one has a massive number of local minima that are all 
!       close to and surrounding the global minimum. This test is useful to ensure that the 
!       algorithm isn't tricked by local minima easily. It manages to find the global minimum 
!       over 97% of the time, 90% to pass.
! 
!   Rastrigin Function (n vars) - This function has 11n extreme local minima equally spaced around
!       the global minimum. This is similar to the Schaffer function in use, but allows the 
!       testing of an unlimited number of variables. It finds the global minimum over 97% with 5 
!       variables, set to 90%. The number of vars can be increased/decreased.
! 
!   Eggholder Function (2 vars) - This one is wild. It has many local minima all over the place, 
!       and most are within 0.1 of the global minimum. On top of this the global minimum is 
!       located at the boundary of the domain of the function. So this one is a great test to 
!       check that algorithms ability to find minimums when the function varies so wildly in 
!       smoothness as well as checking for the ability to locate global minimums at boundary 
!       conditions. The best the algorithm could do is find the global minimum around 86% of 
!       the time, test conditions set to 80%
! 
! ============================================================
module tuner_tests

    use clubb_precision, only: &
        core_rknd
    
    use enhanced_simann, only: &
        esa_driver,  & ! Procedure(s)
        esa_driver_siarry

    implicit none

    public :: tuner_tests_driver

    private :: goldstein_price_test, rastrigin_test, himmelblau_test, &
               schaffer_test, eggholder_test, goldstein_price,       &
               rastrigin, himmelblau, schaffer, eggholder

    private

    logical :: &
        l_esa_siarry = .false.,     & ! Turn on to use Siarry's ESA algorithm
        l_print_outs = .false.        ! Print outs of results, useful for testing


    integer, parameter :: &
        rastrigin_vars = 5,     & ! Number of variables for the Rastrigin function
        samples = 10000           ! Number of times the functions are minimized

    contains

    !============================= DRIVER =============================
    function tuner_tests_driver()

        implicit none

        integer :: tuner_tests_driver
            
        ! initial random number generator
        call init_random

        tuner_tests_driver = 0

        print *, "Running tuner tests"
        print *, "-------------------"
        
        if ( goldstein_price_test() < 0.9 ) then
            print *, "goldstein_price fail"
            tuner_tests_driver = 1
        else 
            print *, "goldstein_price pass"
        end if

        if ( rastrigin_test() < 0.9 ) then
            print *, "rastrigin fail"
            tuner_tests_driver = 1
        else 
            print *, "rastrigin pass"
        end if

        if ( himmelblau_test() < 0.95 ) then
            print *, "himmelblau fail"
            tuner_tests_driver = 1
        else 
            print *, "himmelblau pass"
        end if

        if ( eggholder_test() < 0.7 ) then
            print *, "eggholder fail"
            tuner_tests_driver = 1
        else 
            print *, "eggholder pass"
        end if

        if ( schaffer_test() < 0.9 ) then
            print *, "schaffer fail"
            tuner_tests_driver = 1
        else 
            print *, "schaffer pass"
        end if

        if ( l_print_outs ) then

            l_esa_siarry = .false.

            print *, "esa new"
            print '(A5,F6.2)', "RA: ", rastrigin_test()*100.
            print '(A4,F6.2)', "GP: ", goldstein_price_test()*100.
            print '(A4,F6.2)', "HI: ", himmelblau_test()*100.
            print '(A4,F6.2)', "EG: ", eggholder_test()*100.
            print '(A4,F6.2)', "SC: ", schaffer_test()*100.
            print *,""

            l_esa_siarry = .true.

            print *, "esa siarry"
            print '(A5,F6.2)', "RA: ", rastrigin_test()*100.
            print '(A4,F6.2)', "GP: ", goldstein_price_test()*100.
            print '(A4,F6.2)', "HI: ", himmelblau_test()*100.
            print '(A4,F6.2)', "EG: ", eggholder_test()*100.
            print '(A4,F6.2)', "SC: ", schaffer_test()*100.
            print *,""
        end if


    end function tuner_tests_driver

    !============================= Goldstein-Price =============================
    function goldstein_price_test() 

        implicit none

        real( kind = core_rknd ), dimension(2) :: &
                xinit, xmin, xmax, xopt, xrand

        real( kind = core_rknd ) :: goldstein_price_test, nrgy_opt, &
                intial_temp, max_final_temp, stp_adjst_center, stp_adjst_spread, f_tol

        integer :: &
            pass_count, &
            k
            

        ! setup variables
        xmin = -2._core_rknd
        xmax = 2._core_rknd
        f_tol = 1.e-2_core_rknd
        pass_count = 0

        intial_temp = 450._core_rknd
        max_final_temp = 1.e-10_core_rknd
        stp_adjst_center = 0.58_core_rknd
        stp_adjst_spread = 0.28_core_rknd

        do k = 1, samples

            ! initialize random starting values
            call random_number( xrand )
            xinit = xmin + ( xmax - xmin ) * xrand


            if ( l_esa_siarry ) then

                call esa_driver_siarry( xinit, xmin, xmax, intial_temp, goldstein_price, &
                                        xopt, nrgy_opt )
            else 
                
                ! minimize
                call esa_driver( xinit, xmin, xmax,                     & ! intent(in)
                                 intial_temp, max_final_temp,           & ! intent(in)
                                 xopt, nrgy_opt,                        & ! intent(out)
                                 goldstein_price,                       & ! procedure    
                                 stp_adjst_center, stp_adjst_spread,    & ! optional(in)
                                 f_tol                                  ) ! ^

            end if
           
            ! if x(1) ~= 0 and x(2) ~= -1, then it found the global minimum
            if ( abs( xopt(1) ) <= 1e-2 .and. abs( xopt(2) + 1 ) <= 1e-2 ) then
                pass_count = pass_count + 1
            end if

        end do

        goldstein_price_test = real( pass_count, kind = core_rknd ) &
                             / real( samples, kind = core_rknd )

        return

    end function goldstein_price_test


    ! Goldstein-Price function, 2 vars
    real function goldstein_price( x )

        implicit none

        real, dimension(:), intent(in) :: x

        ! Goldstein-Price, 2 vars, [-2,2]
        goldstein_price = ( 1 + (x(1) + x(2) + 1 )**2 * (19 - 14*x(1) + 3 * x(1)**2 - 14 *   &
                          x(2) + 6 * x(1) * x(2) + 3 * x(2)**2) ) * ( ( 30 + ( 2 * x(1) - 3  &
                          * x(2) )**2 * ( 18 - 32 * x(1) + 12 * x(1)**2 + 48 * x(2) - 36 *   &
                          x(1) * x(2) + 27 * x(2)**2 ) ) )

    end function goldstein_price








    !============================= Rastrigin =============================
    function rastrigin_test() 

        implicit none

        real( kind = core_rknd ), dimension(rastrigin_vars) :: &
                xinit, xmin, xmax, xopt, xrand

        real( kind = core_rknd ) :: rastrigin_test, nrgy_opt, &
                intial_temp, max_final_temp, stp_adjst_center, stp_adjst_spread, f_tol

        integer :: &
            pass_count, &
            k

        ! setup variables
        xmin = -5.12_core_rknd
        xmax = 5.12_core_rknd
        f_tol = 1.e-2_core_rknd
        pass_count = 0

        intial_temp = 20.
        max_final_temp = 1.e-10_core_rknd
        stp_adjst_center = 0.694
        stp_adjst_spread = 0.1
        

        do k = 1, samples

            ! initialize random starting values
            call random_number( xrand )
            xinit = xmin + ( xmax - xmin ) * xrand

            if ( l_esa_siarry ) then

                call esa_driver_siarry( xinit, xmin, xmax, intial_temp, rastrigin, xopt, nrgy_opt )
            else 
                ! minimize
                call esa_driver( xinit, xmin, xmax,                     & ! intent(in)
                                 intial_temp, max_final_temp,           & ! intent(in)
                                 xopt, nrgy_opt,                        & ! intent(out)
                                 rastrigin,                             & ! procedure    
                                 stp_adjst_center, stp_adjst_spread,    & ! optional(in)
                                 f_tol                                  ) ! ^
            end if

           
            ! if all x ~= 0 then global minimum found
            if ( all( abs(xopt) < 0.1 )  ) then
                pass_count = pass_count + 1
            end if

        end do

        rastrigin_test = real( pass_count, kind = core_rknd ) &
                       / real( samples, kind = core_rknd )

        return

    end function rastrigin_test

    ! Rastrigin function, n vars
    real function rastrigin( x )

        implicit none

        real, dimension(:), intent(in) :: x

        integer :: k 
    
        rastrigin = 10 * size(x)
        do k = 1, size(x)
            rastrigin = rastrigin + x(k)**2 - 10 * cos( 2 * 3.14159265358979323 * x(k) )       
        end do

    end function rastrigin




    !============================= Himmelblau =============================
    function himmelblau_test( ) 

        implicit none

        real( kind = core_rknd ), dimension(2) :: &
                xinit, xmin, xmax, xopt, xrand

        real( kind = core_rknd ) :: himmelblau_test, nrgy_opt, &
                intial_temp, max_final_temp, stp_adjst_center, stp_adjst_spread, f_tol

        integer :: &
            pass_count, &
            k
            

        ! setup variables
        xmin = -5._core_rknd
        xmax = 5._core_rknd
        f_tol = 1.e-2_core_rknd
        pass_count = 0

        intial_temp = 10._core_rknd
        max_final_temp = 1.e-10_core_rknd
        stp_adjst_center = 0.96_core_rknd
        stp_adjst_spread = 0.8_core_rknd
        
        do k = 1, samples

            ! initialize random starting values
            call random_number( xrand )
            xinit = xmin + ( xmax - xmin ) * xrand

            if ( l_esa_siarry ) then

                call esa_driver_siarry( xinit, xmin, xmax, intial_temp, himmelblau, xopt, nrgy_opt )
            else 
                ! minimize
                call esa_driver( xinit, xmin, xmax,                     & ! intent(in)
                                 intial_temp, max_final_temp,           & ! intent(in)
                                 xopt, nrgy_opt,                        & ! intent(out)
                                 himmelblau,                            & ! procedure    
                                 stp_adjst_center, stp_adjst_spread,    & ! optional(in)
                                 f_tol                                  ) ! ^
            end if
           
            
            if ( nrgy_opt < 1.e-2_core_rknd ) then
                pass_count = pass_count + 1
            end if

        end do

        himmelblau_test = real( pass_count, kind = core_rknd ) &
                        / real( samples, kind = core_rknd )

        return

    end function himmelblau_test

    ! Himmelblau function, 2 var
    real function himmelblau( x )

        implicit none

        real, dimension(:), intent(in) :: x
        
        himmelblau = ( x(1)**2 + x(2) - 11 )**2 + ( x(1) + x(2)**2 - 7 )**2

        return

    end function himmelblau


    





    !============================= Eggholder =============================
    function eggholder_test() 
       
        implicit none

        real( kind = core_rknd ), dimension(2) :: &
                xinit, xmin, xmax, xopt, xrand

        real( kind = core_rknd ) :: eggholder_test, nrgy_opt, &
                intial_temp, max_final_temp, stp_adjst_center, stp_adjst_spread, f_tol

        integer :: &
            pass_count, &
            k
            

        ! setup variables
        xmin = -512._core_rknd
        xmax = 512._core_rknd
        f_tol = 1.e-2_core_rknd
        pass_count = 0

        intial_temp = 85000._core_rknd
        max_final_temp = 1.e-10_core_rknd
        stp_adjst_center = 0.128_core_rknd
        stp_adjst_spread = 0.87_core_rknd
        
        do k = 1, samples

            ! initialize random starting values
            call random_number( xrand )
            xinit = xmin + ( xmax - xmin ) * xrand


            if ( l_esa_siarry ) then

                call esa_driver_siarry( xinit, xmin, xmax, intial_temp, eggholder, xopt, nrgy_opt )
            else 
                ! minimize
                call esa_driver( xinit, xmin, xmax,                     & ! intent(in)
                                 intial_temp, max_final_temp,           & ! intent(in)
                                 xopt, nrgy_opt,                        & ! intent(out)
                                 eggholder,                             & ! procedure    
                                 stp_adjst_center, stp_adjst_spread,    & ! optional(in)
                                 f_tol                                  ) ! ^
                
            end if
           
            if ( abs( xopt(1) - 512. ) < 1e-2 .and. abs( xopt(2) - 404.23 ) < 1 ) then
                pass_count = pass_count + 1
            end if

        end do

        eggholder_test = real( pass_count, kind = core_rknd ) &
                       / real( samples, kind = core_rknd )

        return

    end function eggholder_test

    ! Eggholder function, 2 var
    real function eggholder( x )

        implicit none

        real, dimension(:), intent(in) :: x
    
        eggholder = -( x(2) + 47 ) * sin( sqrt( abs( x(1)/2.0 + x(2) + 47 ) ) ) &
                    - x(1) * sin( sqrt( abs( x(1) - x(2) - 47 ) ) )

    end function eggholder
    
    

    !============================= Schaffer =============================
    function schaffer_test() 
       
        implicit none

        real( kind = core_rknd ), dimension(2) :: &
                xinit, xmin, xmax, xopt, xrand

        real( kind = core_rknd ) :: schaffer_test, nrgy_opt, &
                intial_temp, max_final_temp, stp_adjst_center, stp_adjst_spread, f_tol

        integer :: &
            pass_count, &
            k
            

        ! setup variables
        xmin = -100._core_rknd
        xmax = 100._core_rknd
        f_tol = 1.e-2_core_rknd
        pass_count = 0

        intial_temp = 20._core_rknd
        max_final_temp = 1.e-10_core_rknd
        stp_adjst_center = 0.8_core_rknd
        stp_adjst_spread = 0.4_core_rknd
       

        do k = 1, samples

            ! initialize random starting values
            call random_number( xrand )
            xinit = xmin + ( xmax - xmin ) * xrand

            if ( l_esa_siarry ) then

                call esa_driver_siarry( xinit, xmin, xmax, intial_temp, schaffer, xopt, nrgy_opt )
            else 
                ! minimize
                call esa_driver( xinit, xmin, xmax,                     & ! intent(in)
                                 intial_temp, max_final_temp,           & ! intent(in)
                                 xopt, nrgy_opt,                        & ! intent(out)
                                 schaffer,                              & ! procedure    
                                 stp_adjst_center, stp_adjst_spread,    & ! optional(in)
                                 f_tol                                  ) ! ^
            end if
           
            ! if all x ~= 0 then global minimum found
            if ( nrgy_opt <= 1e-5 ) then
                pass_count = pass_count + 1
            end if

        end do

        schaffer_test = real( pass_count, kind = core_rknd ) &
                       / real( samples, kind = core_rknd )

        return

    end function schaffer_test

    ! Schaffer function, 2 var
    real function schaffer( x )

        implicit none

        real, dimension(:), intent(in) :: x
    
        schaffer = 0.5 + ( sin( x(1)**2 - x(2)**2 )**2 - 0.5 ) &
                         / ( 1 + 0.001 * ( x(1)**2 + x(2)**2 ) )**2

    end function schaffer
            

    subroutine init_random

        implicit none

        integer :: i, n, clock
        integer, dimension(:), allocatable :: seed
      
        call random_seed(size = n)
        allocate(seed(n))
      
        call system_clock(count=clock)
      
        seed = clock + 37 * (/ (i - 1, i = 1, n) /)
        call random_seed(PUT = seed)
      
        deallocate(seed)

    end subroutine init_random

end module tuner_tests
