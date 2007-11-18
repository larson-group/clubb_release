!-----------------------------------------------------------------------
! $Id: param_index.F,v 1.4 2007-11-18 17:17:27 griffinb Exp $
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
     .  nparams = 49 ! Total tunable parameters

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
     .  iC14            = 29

        integer, parameter ::
     .  ic_K            = 30,
     .  ic_K1           = 31,
     .  inu1            = 32,
     .  ic_K2           = 33,
     .  inu2            = 34,
     .  ic_K6           = 35,
     .  inu6            = 36,
     .  ic_K8           = 37,
     .  inu8            = 38,
     .  inu9            = 39,
     .  ic_Krrm         = 40,
     .  inu_r           = 41

        integer, parameter ::
     .  igamma_coef     = 42, 
     .  igamma_coefb    = 43,
     .  igamma_coefc    = 44,
     .  imu             = 45,
     .  ibeta           = 46,
     .  ilmin_coef      = 47,
     .  itaumin         = 48,
     .  itaumax         = 49

        end module param_index
!-----------------------------------------------------------------------
