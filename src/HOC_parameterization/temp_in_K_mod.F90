! $Id: temp_in_K_mod.F90,v 1.1 2008-07-22 16:04:30 faschinj Exp $ 

         module T_in_K_mod
         
         implicit none
         
         private ! Default scope

         public :: thlm2T_in_K
         
         contains

!-----------------------------------------------------------------------
         elemental function thlm2T_in_K( thlm, exner, rcm )  & 
         result( T_in_K )

!        Description:
!        Calculates absolute temperature from liquid water potential
!        temperature.  (Does not include ice.)

!        References:  Cotton and Anthes (1989), "Storm and Cloud Dynamics",
!                        Eqn. (2.51). 
!-----------------------------------------------------------------------
         use constants, only: & 
         ! Variable(s) 
         Cp,  & ! Dry air specific heat at constant p [J/kg/K]
         Lv  ! Latent heat of vaporization         [J/kg]

         implicit none

         ! Input 
         real, intent(in) :: & 
         thlm,   & ! Liquid potential temperature  [K]
         exner,  & ! Exner function                [-]
         rcm    ! Liquid water mixing ratio     [kg/kg]

         real :: & 
         T_in_K ! Result temperature [K]

         T_in_K = thlm * exner + Lv * rcm / Cp

         return
         end function thlm2T_in_K
!-----------------------------------------------------------------------

         end module T_in_K_mod
