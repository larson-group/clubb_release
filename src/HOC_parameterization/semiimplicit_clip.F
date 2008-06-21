!-----------------------------------------------------------------------
! $Id: semiimplicit_clip.F,v 1.1 2008-06-21 22:04:06 griffinb Exp $
!===============================================================================
        module semiimplicit_clip

!       Description of the semi-implicit clipping code:
!       The semi-implicit clipping code based on an upper threshold
!       and/or a lower threshold value for variable f.
!
!       The semi-implicit clipping code is used when the value of 
!       variable f should not exceed the designated threshold(s) when it
!       is advanced to timestep index (t+1).
!
!
!       Clipping at an Upper Threshold:
!
!       When there is an upper threshold to be applied, the equation for
!       the clipped value of the variable, f_clipped, is:
!
!       f_clipped(t+1) = min(f_unclipped(t+1),upper_threshold)
!                      = ( f_unclipped(t+1) - upper_threshold )
!                          * H(upper_threshold-f_unclipped(t+1))
!                        + upper_threshold;
!
!       where f_unclipped is the value of the variable without 
!       clipping, and H(upper_threshold-f_unclipped(t+1)) is the 
!       Heaviside Step function.  The clipping term is turned into a 
!       time tendency term, such that:
!
!       (df/dt)_clipping = (1/dt_clip) * f_diff(t+1) 
!                                      * [ 1 - H(f_diff(t+1)) ];
!
!       where dt_clip is the time scale for the clipping term, and 
!       f_diff is defined as:
!
!       f_diff = upper_threshold - f_unclipped.
!
!       In order to solve for f_unclipped (and f_diff) at timestep 
!       index (t+1), the clipping term must be linearized.  A Taylor
!       Series expansion (truncated after the first derivative term) of
!       H(f_diff) around f_diff = f_diff(t) is used to linearize the
!       term.  However, the Heaviside Step function, H(f_diff), is not 
!       differentiable at all points.  Therefore, it is approximated by
!       a sigmoidal function, A_fnc, (which has first derivative B_fnc),
!       resulting in:
!
!       H(f_diff(t+1)) approx.= A_fnc(f_diff(t+1))
!          = A_fnc(f_diff(t)) 
!            + B_fnc(f_diff(t)) * ( f_diff(t+1) - f_diff(t) ).
!       
!       The approximation is substituted into the (df/dt)_clipping 
!       equation.  The rate of change of variable f due to clipping with
!       the upper threshold is:
!
!       (df/dt)_clipping
!       = (1/dt_clip)
!        * {   B_fnc * f_diff(t) * f_diff(t)
!            + [ ( 1 - A_fnc ) - B_fnc * f_diff(t) ] * upper_threshold
!            - [ ( 1 - A_fnc ) - B_fnc * f_diff(t) ] * f_unclipped(t+1) }.
!
!       The implicit (LHS) portion of the equation for clipping with the
!       upper threshold is:
!
!       - (1/dt_clip)
!        * [ ( 1 - A_fnc ) - B_fnc * f_diff(t) ] * f_unclipped(t+1).
!
!       Note:  When the term is brought over to the left-hand side, the
!              sign is reversed and the leading "-" in front of the 
!              term is changed to a "+".
!
!       The explicit (RHS) portion of the equation for clipping with the
!       upper threshold is:
!
!       + (1/dt_clip)
!        * {   B_fnc * f_diff(t) * f_diff(t)
!            + [ ( 1 - A_fnc ) - B_fnc * f_diff(t) ] * upper_threshold }.
!       
!       Timestep index (t) stands for the index of the current timestep,
!       while timestep index (t+1) stands for the index of the next 
!       timestep, which is being advanced to in solving the d(f)/dt 
!       equation.
!
!
!       Clipping at a Lower Threshold:
!
!       When there is a lower threshold to be applied, the equation for
!       the clipped value of the variable, f_clipped, is:
!
!       f_clipped(t+1) = max(f_unclipped(t+1),lower_threshold)
!                      = ( f_unclipped(t+1) - lower_threshold )
!                          * H(f_unclipped(t+1)-lower_threshold)
!                        + lower_threshold.
!
!       The clipping term is turned into a time tendency term, such that:
!
!       (df/dt)_clipping = - (1/dt_clip) * f_diff(t+1) 
!                                        * [ 1 - H(f_diff(t+1)) ];
!
!       where dt_clip is the time scale for the clipping term, and 
!       f_diff is defined as:
!
!       f_diff = f_unclipped - lower_threshold.
!
!       The linearization process is the same for the lower threshold as
!       it is for the upper threshold.  The formulas for A_fnc and B_fnc
!       are the same, but the values (based on f_diff) are different.  
!       The rate of change of variable f due to clipping with the lower 
!       threshold is:
!
!       (df/dt)_clipping
!       = - (1/dt_clip)
!        * {   B_fnc * f_diff(t) * f_diff(t)
!            - [ ( 1 - A_fnc ) - B_fnc * f_diff(t) ] * lower_threshold
!            + [ ( 1 - A_fnc ) - B_fnc * f_diff(t) ] * f_unclipped(t+1) }.
!
!       The implicit (LHS) portion of the equation for clipping with the
!       lower threshold is:
!
!       - (1/dt_clip)
!        * [ ( 1 - A_fnc ) - B_fnc * f_diff(t) ] * f_unclipped(t+1).
!
!       Note:  When the term is brought over to the left-hand side, the
!              sign is reversed and the leading "-" in front of the 
!              term is changed to a "+".
!
!       The explicit (RHS) portion of the equation for clipping with the
!       lower threshold is:
!
!       - (1/dt_clip)
!        * {   B_fnc * f_diff(t) * f_diff(t)
!            - [ ( 1 - A_fnc ) - B_fnc * f_diff(t) ] * lower_threshold }.
!
!       All variables in these equations are on the same levels as the
!       variable f.

!       References:
!-----------------------------------------------------------------------

        implicit none

        private

        public :: semiimp_clip_lhs,
     .            semiimp_clip_rhs

        private :: compute_clip_lhs,
     .             compute_functions_A_B

        contains

!===============================================================================
        function semiimp_clip_lhs( dt_clip, f_unclipped, 
     .                             lupper_thresh, upper_threshold, 
     .                             llower_thresh, lower_threshold )
     .  result( lhs )

!       Description:
!       The implicit portion of the semi-implicit clipping code.
!
!       The implicit (LHS) portion of the equation for clipping with the
!       upper threshold is:
!
!       - (1/dt_clip)
!        * [ ( 1 - A_fnc ) - B_fnc * f_diff(t) ] * f_unclipped(t+1).
!       
!       The implicit (LHS) portion of the equation for clipping with the
!       lower threshold is:
!
!       - (1/dt_clip)
!        * [ ( 1 - A_fnc ) - B_fnc * f_diff(t) ] * f_unclipped(t+1).
!
!       Note:  When either term is brought over to the left-hand side, 
!              the sign is reversed and the leading "-" in front of 
!              either term is changed to a "+".
!
!       Timestep index (t) stands for the index of the current timestep, 
!       while timestep index (t+1) stands for the index of the next 
!       timestep, which is being advanced to in solving the d(f)/dt 
!       equation.
!
!       While the formulas are the same for both the upper threshold
!       and the lower threshold, the values of A_fnc, B_fnc, and f_diff
!       will differ between the two thresholds.
!
!       The overall implicit (LHS) portion for the clipping term is the
!       sum of the implicit portion from the upper threshold and the 
!       implicit portion from the lower threshold.

!       References:
!-----------------------------------------------------------------------

        use stats_prec, only: 
     .      time_prec ! Variable(s)

        implicit none

        ! Input Variables
        real(kind=time_prec), intent(in) :: 
     .  dt_clip            ! Time scale for semi-implicit clipping term.   [s]

        real, intent(in) ::
     .  f_unclipped,       ! The unclipped value of variable f at timestep (t).
     .  upper_threshold,   ! Greatest allowable value of variable f.
     .  lower_threshold    ! Smallest allowable value of variable f.

        logical, intent(in) ::
     .  lupper_thresh,   ! Flag for having an upper threshold value.
     .  llower_thresh    ! Flag for having a lower threshold value.

        ! Return Variable
        real :: lhs

        ! Local Variables
        real ::
     .  f_diff,    ! Difference between the threshold value and f_unclipped.
     .  A_fnc,     ! Sigmoidal function approximation to Heaviside Step function.
     .  B_fnc,     ! Derivative w/ respect to f_diff of sigmoidal function, A_fnc.
     .  lhs_upper, ! Contribution of upper threshold to implicit portion (LHS).
     .  lhs_lower  ! Contribution of lower threshold to implicit portion (LHS).


        ! Upper Threshold
        if ( lupper_thresh ) then

           ! f_diff is the difference between the threshold value and 
           ! f_unclipped.  In regards to the upper threshold, it is 
           ! defined as upper_threshold - f_unclipped.
           f_diff = upper_threshold - f_unclipped

           ! Compute the values of functions A_fnc and B_fnc evaluated 
           ! at f_diff(t) for the upper threshold.
           call compute_functions_A_B( f_diff, A_fnc, B_fnc )

           ! Compute the implicit (LHS) contribution from clipping for 
           ! the upper threshold. 
           lhs_upper = compute_clip_lhs( dt_clip, f_diff, A_fnc, B_fnc )

        else

           lhs_upper = 0.0

        endif


        ! Lower Threshold
        if ( llower_thresh ) then

           ! f_diff is the difference between the threshold value and 
           ! f_unclipped.  In regards to the lower threshold, it is 
           ! defined as f_unclipped - lower_threshold.
           f_diff = f_unclipped - lower_threshold

           ! Compute the values of functions A_fnc and B_fnc evaluated 
           ! at f_diff(t) for the lower threshold.
           call compute_functions_A_B( f_diff, A_fnc, B_fnc )

           ! Compute the implicit (LHS) contribution from clipping for 
           ! the lower threshold. 
           lhs_lower = compute_clip_lhs( dt_clip, f_diff, A_fnc, B_fnc )

        else

           lhs_lower = 0.0

        endif


        ! Total implicit (LHS) contribution to clipping.
        ! Main diagonal: [ x f_unclipped(k,<t+1>) ]
        lhs = lhs_upper + lhs_lower


        end function semiimp_clip_lhs

!===============================================================================
        pure function compute_clip_lhs( dt_clip, f_diff, A_fnc, B_fnc )
     .  result( lhs_contribution )

!       Description:
!       Calculation of the implicit portion of the semi-implicit 
!       clipping term.
!
!       The implicit portion of the semi-implicit clipping term is:
!
!       - (1/dt_clip) * [ ( 1 - A_fnc ) - B_fnc * f_diff(t) ] 
!                     * f_unclipped(t+1).
!
!       Note:  When the term is brought over to the left-hand side, the
!              sign is reversed and the leading "-" in front of the 
!              term is changed to a "+".
!
!       Timestep index (t) stands for the index of the current timestep, 
!       while timestep index (t+1) stands for the index of the next 
!       timestep, which is being advanced to in solving the d(f)/dt 
!       equation.

!       References:
!-----------------------------------------------------------------------

        use stats_prec, only: 
     .      time_prec ! Variable(s)

        implicit none

        ! Input Variables
        real(kind=time_prec), intent(in) :: 
     .  dt_clip ! Time scale for semi-implicit clipping term.   [s]

        real, intent(in) ::
     .  f_diff, ! Difference between the threshold value and f_unclipped.
     .  A_fnc,  ! Sigmoidal function approximation to Heaviside Step function.
     .  B_fnc   ! Derivative w/ respect to f_diff of sigmoidal function, A_fnc.

        ! Return Variable
        real :: lhs_contribution


        ! Main diagonal: [ x f_unclipped(k,<t+1>) ]
        lhs_contribution
     .  = + (1.0/dt_clip) * ( ( 1.0 - A_fnc ) - B_fnc * f_diff )


        end function compute_clip_lhs

!===============================================================================
        function semiimp_clip_rhs( dt_clip, f_unclipped, 
     .                             lupper_thresh, upper_threshold, 
     .                             llower_thresh, lower_threshold )
     .  result( rhs )

!       Description:
!       The explicit portion of the semi-implicit clipping code.
!
!       The explicit (RHS) portion of the equation for clipping with the
!       upper threshold is:
!
!       + (1/dt_clip)
!        * {   B_fnc * f_diff(t) * f_diff(t)
!            + [ ( 1 - A_fnc ) - B_fnc * f_diff(t) ] * upper_threshold }.
!       
!       The explicit (RHS) portion of the equation for clipping with the
!       lower threshold is:
!
!       - (1/dt_clip)
!        * {   B_fnc * f_diff(t) * f_diff(t)
!            - [ ( 1 - A_fnc ) - B_fnc * f_diff(t) ] * lower_threshold }.
!
!       Timestep index (t) stands for the index of the current timestep.
!
!       The overall explicit (RHS) portion for the clipping term is the
!       sum of the explicit portion from the upper threshold and the 
!       explicit portion from the lower threshold.

!       References:
!-----------------------------------------------------------------------

        use stats_prec, only: 
     .      time_prec ! Variable(s)

        implicit none

        ! Input Variables
        real(kind=time_prec), intent(in) :: 
     .  dt_clip            ! Time scale for semi-implicit clipping term.   [s]

        real, intent(in) ::
     .  f_unclipped,       ! The unclipped value of variable f at timestep (t).
     .  upper_threshold,   ! Greatest allowable value of variable f.
     .  lower_threshold    ! Smallest allowable value of variable f.

        logical, intent(in) ::
     .  lupper_thresh,   ! Flag for having an upper threshold value.
     .  llower_thresh    ! Flag for having a lower threshold value.

        ! Return Variable
        real :: rhs

        ! Local Variables
        real ::
     .  f_diff,    ! Difference between the threshold value and f_unclipped.
     .  A_fnc,     ! Sigmoidal function approximation to Heaviside Step function.
     .  B_fnc,     ! Derivative w/ respect to f_diff of sigmoidal function, A_fnc.
     .  rhs_upper, ! Contribution of upper threshold to explicit portion (RHS).
     .  rhs_lower  ! Contribution of lower threshold to explicit portion (RHS).


        ! Upper Threshold
        if ( lupper_thresh ) then

           ! f_diff is the difference between the threshold value and 
           ! f_unclipped.  In regards to the upper threshold, it is 
           ! defined as upper_threshold - f_unclipped.
           f_diff = upper_threshold - f_unclipped

           ! Compute the values of functions A_fnc and B_fnc evaluated 
           ! at f_diff(t) for the upper threshold.
           call compute_functions_A_B( f_diff, A_fnc, B_fnc )

           ! Compute the explicit (RHS) contribution from clipping for 
           ! the upper threshold. 
           rhs_upper 
     .     = + (1.0/dt_clip) 
     .        * (   B_fnc * f_diff * f_diff
     .            + ( ( 1.0 - A_fnc ) - B_fnc * f_diff )
     .              * upper_threshold   )

        else

           rhs_upper = 0.0

        endif


        ! Lower Threshold
        if ( llower_thresh ) then

           ! f_diff is the difference between the threshold value and 
           ! f_unclipped.  In regards to the lower threshold, it is 
           ! defined as f_unclipped - lower_threshold.
           f_diff = f_unclipped - lower_threshold

           ! Compute the values of functions A_fnc and B_fnc evaluated 
           ! at f_diff(t) for the lower threshold.
           call compute_functions_A_B( f_diff, A_fnc, B_fnc )

           ! Compute the explicit (RHS) contribution from clipping for 
           ! the lower threshold. 
           rhs_lower 
     .     = - (1.0/dt_clip) 
     .        * (   B_fnc * f_diff * f_diff
     .            - ( ( 1.0 - A_fnc ) - B_fnc * f_diff )
     .              * lower_threshold   )

        else

           rhs_lower = 0.0

        endif


        ! Total explicit (RHS) contribution to clipping.
        rhs = rhs_upper + rhs_lower


        end function semiimp_clip_rhs

!===============================================================================
        subroutine compute_functions_A_B( f_diff, A_fnc, B_fnc )

!       Description:
!       This subroutine computes the values of two functions used in
!       semi-implicit clipping.  Both of the functions are based on the 
!       values of f_diff(t) and the constant parameter alpha_clip.  One 
!       function is A_fnc, a sigmoidal function that is used to 
!       approximate the Heaviside Step function.  The other function is 
!       B_fnc, the derivative with respect to f_diff of the sigmoidal 
!       function A_fnc.
!
!       The equation for A_fnc is:
!
!       A(alpha_clip,f_diff(t)) = 1 / [ 1 + e^(-alpha_clip*f_diff(t)) ];
!
!       while the equation for B_fnc is:
!
!       B(alpha_clip,f_diff(t)) 
!            = [ alpha_clip * e^(-alpha_clip*f_diff(t)) ]
!              / [ 1 + e^(-alpha_clip*f_diff(t)) ]^2;
!
!       where timestep index (t) stands for the index of the current 
!       timestep.

!       References:
!-----------------------------------------------------------------------

        implicit none

        ! Constant parameter
        real, parameter :: alpha_clip = 10.0

        ! Input Variable
        real, intent(in) ::
     .  f_diff   ! Difference between the threshold value and f_unclipped.

        ! Output Variables
        real, intent(out) ::
     .  A_fnc,  ! Sigmoidal function approximation to Heaviside Step function.
     .  B_fnc   ! Derivative w/ respect to f_diff of sigmoidal function, A_fnc.


        ! A_fnc is a sigmoidal function that approximates the Heaviside 
        ! Step function.  This is needed because the Heaviside Step 
        ! function is not differentiable at point f_diff = 0.  The 
        ! sigmoidal function is differentiable at all points.  It is
        ! evaluated for f_diff at timestep index (t).
        A_fnc = 1.0 / ( 1.0 + exp( -alpha_clip * f_diff ) )

        ! B_fnc is the derivative with respect to f_diff of sigmoidal 
        ! function A_fnc.  It is evaluated for f_diff at timestep 
        ! index (t).
        B_fnc = ( alpha_clip * exp( -alpha_clip * f_diff ) )
     .          / ( 1.0 + exp( -alpha_clip * f_diff ) )**2


        end subroutine compute_functions_A_B

!===============================================================================

        end module semiimplicit_clip
