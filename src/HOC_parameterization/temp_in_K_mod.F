         module temp_in_K_mod
         
         implicit none
         
         private

         public :: thlm2temp_in_K
         
         contains

!-----------------------------------------------------------------------
         elemental function thlm2temp_in_K( thlm, exner, rcm ) 
     .   result( temp_in_K )

!        Description:
!        Calculates absolute temperature from liquid potential
!        temperature.

!        References:
!-----------------------------------------------------------------------
         use constants, only: 
     .   Cp, ! Dry air specific heat at constant p [J/kg/K]
     .   Lv  ! Latent heat of vaporization         [J/kg]

         implicit none

         ! Input 
         real, intent(in) ::
     .   thlm,  ! Liquid potential temperature  [K]
     .   exner, ! Exner function                [-]
     .   rcm    ! Liquid water mixing ratio     [kg/kg]

         real ::
     .   temp_in_K ! Result temperature [K]

         temp_in_K = thlm * exner + Lv * rcm / Cp

         return
         end function thlm2temp_in_K
!-----------------------------------------------------------------------

         end module temp_in_K_mod
