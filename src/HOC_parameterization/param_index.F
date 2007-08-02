!-----------------------------------------------------------------------
! $Id: param_index.F,v 1.3 2007-08-02 22:51:40 griffinb Exp $
        module param_index

!       Description:
!       Since f90/95 lacks enumeration, we're stuck numbering each
!       parameter by hand like this.

!       Adding new parameters is relatively simple.  First, the
!       parameter should be added in the common block of the parameters
!       module so it can be used in other parts of the code. Each
!       variable needs a unique number in this module, and nparams must
!       be incremented for the new variable.  Next, the params_list
!       variable in module parameters should have new variable added to
!       it.  The subroutines pack_parameters and uppack_parameters will
!       need to have the variable added to their list, but the order
!       doesn't actually matter, since the i variables in here determine
!       where in the params vector the number is placed.
!       Finally, the namelists initvars and initspread will need to
!       have the parameter added to them.
!-----------------------------------------------------------------------

        implicit none

        integer, parameter :: 
     .  nparams = 45 ! Total tunable parameters

        integer, parameter ::
     .  iC1             =  1, 
     .  iC1b            =  2,
     .  iC1c            =  3, 
     .  iC2             =  4, 
     .  iC2b            =  5,
     .  iC2c            =  6,
     .  iC2rt           =  7,
     .  iC2thl          =  8,
     .  iC2rtthl        =  9,
     .  iC4             = 10, 
     .  iC5             = 11, 
     .  iC6rt           = 12,
     .  iC6rtb          = 13,
     .  iC6rtc          = 14,  
     .  iC6thl          = 15, 
     .  iC6thlb         = 16,
     .  iC6thlc         = 17,
     .  iC7             = 18,
     .  iC7b            = 19,
     .  iC7c            = 20,
     .  iC8             = 21,
     .  iC8b            = 22,
     .  iC10            = 23,
     .  iC11            = 24,
     .  iC11b           = 25,
     .  iC11c           = 26,
     .  iC12            = 27,
     .  iC13            = 28,
     .  iC14            = 44

        integer, parameter ::
     .  inu1            = 29,
     .  inu2            = 30,
     .  inu6            = 31,
     .  inu8            = 32,
     .  inu9            = 45,
     .  inu_r           = 33

        integer, parameter ::
     .  igamma_coef     = 34, 
     .  igamma_coefb    = 35,
     .  igamma_coefc    = 36,
     .  imu             = 37,
     .  ic_K            = 38,
     .  ic_Krrm         = 39,
     .  ibeta           = 40,
     .  ilmin_coef      = 41,
     .  itaumin         = 42,
     .  itaumax         = 43

        end module param_index
!-----------------------------------------------------------------------
