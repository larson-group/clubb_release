!-----------------------------------------------------------------------
! $Id$
module parameter_indices

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

  private ! Default Scope

  integer, parameter, public ::  & 
    nparams = 53 ! Total tunable parameters

!***************************************************************
!                    ***** IMPORTANT *****
! If you change the order of these parameters, you will need to
! change the order of params_list as well or the tuner will
! break!
!                    ***** IMPORTANT *****
!***************************************************************

  integer, parameter, public :: & 
    iC1             =  1, & 
    iC1b            =  2, & 
    iC1c            =  3, & 
    iC2             =  4, & 
    iC2b            =  5, & 
    iC2c            =  6, & 
    iC2rt           =  7, & 
    iC2thl          =  8, & 
    iC2rtthl        =  9, & 
    iC4             = 10, & 
    iC5             = 11, & 
    iC6rt           = 12, & 
    iC6rtb          = 13, & 
    iC6rtc          = 14, & 
    iC6thl          = 15, & 
    iC6thlb         = 16, & 
    iC6thlc         = 17, & 
    iC7             = 18, & 
    iC7b            = 19, & 
    iC7c            = 20, & 
    iC8             = 21, & 
    iC8b            = 22, & 
    iC10            = 23, & 
    iC11            = 24, & 
    iC11b           = 25, & 
    iC11c           = 26, & 
    iC12            = 27, & 
    iC13            = 28, & 
    iC14            = 29, &
    iC15            = 30

  integer, parameter, public :: & 
    ic_K            = 31, & 
    ic_K1           = 32, & 
    inu1            = 33, & 
    ic_K2           = 34, & 
    inu2            = 35, & 
    ic_K6           = 36, & 
    inu6            = 37, & 
    ic_K8           = 38, & 
    inu8            = 39, & 
    ic_K9           = 40, & 
    inu9            = 41, & 
    ic_Krrainm      = 42, & 
    inu_r           = 43, & 
    ic_Ksqd         = 44, &
    inu_hd          = 45

  integer, parameter, public :: & 
    igamma_coef     = 46, & 
    igamma_coefb    = 47, & 
    igamma_coefc    = 48, & 
    imu             = 49, & 
    ibeta           = 50, & 
    ilmin_coef      = 51, & 
    itaumin         = 52, & 
    itaumax         = 53

end module parameter_indices
!-----------------------------------------------------------------------
