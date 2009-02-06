!-----------------------------------------------------------------------
! $Id$
module parameter_indices

!       Description:
!       Since f90/95 lacks enumeration, we're stuck numbering each
!       parameter by hand like this.

!       Adding new parameters is relatively simple.  First, the
!       parameter should be added in the common block of the parametersle
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
nparams = 52 ! Total tunable parameters

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
iC14            = 29

integer, parameter, public :: & 
ic_K            = 30, & 
ic_K1           = 31, & 
inu1            = 32, & 
ic_K2           = 33, & 
inu2            = 34, & 
ic_K6           = 35, & 
inu6            = 36, & 
ic_K8           = 37, & 
inu8            = 38, & 
ic_K9           = 39, & 
inu9            = 40, & 
ic_Krrainm      = 41, & 
inu_r           = 42, & 
ic_Ksqd         = 43, &
inu_hd          = 44

integer, parameter, public :: & 
igamma_coef     = 45, & 
igamma_coefb    = 46, & 
igamma_coefc    = 47, & 
imu             = 48, & 
ibeta           = 49, & 
ilmin_coef      = 50, & 
itaumin         = 51, & 
itaumax         = 52

end module parameter_indices
!-----------------------------------------------------------------------
