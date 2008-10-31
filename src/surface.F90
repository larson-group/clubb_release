!$Id$
module surface

implicit none

public :: prognose_soil_T_in_K, initialize_surface

real, private :: initial_veg_t_sfc, initial_t_sfc, initial_deep_t_sfc

private


contains
!----------------------------------------------------------------------
SUBROUTINE prognose_soil_T_in_K( BTIME, TIME, DT, ITF, ITL, R0, FS, FSIN, FLDO, FLUP, WTQ, &
                   veg_T_in_K, sfc_soil_T_in_K, deep_soil_T_in_K )
!
!
!     IN THIS SUBROUTINE TQ AND QW AT THE SURFACE ARE CALCULATED AS A
!     FUNCTION OF THE SURFACE TEMPERATURE (sfc_soil_T_in_K). sfc_soil_T_in_K IS CALCULATED  FROM
!     THE SURFACE ENERGY BUDGET.
!
! ************************************** 2
!               *                   *
!               *                   *
!             --*         +         *--  1
!               *                   *
! *** SURFACE ************+************* 0
!                        sfc_soil_T_in_K              L (EVEL)
!
!
!     HEAT CONDUCTION (IN A HOMOGENEOUS MEDIUM) CAN BE DESCRIBED BY
!     THE EQUATION:
!
!              D(sfc_soil_T_in_K)/DT=KS*D(D(sfc_soil_T_in_K)/DZ)/DZ
!
!     IN WHICH KS IS THE SOIL THERMAL DIFFUSITY.
!     WE CONSIDER A SEMI HALF INFINITE MEDIUM, INITIALLY AT THE
!     CONSTANT TEMPERATURE deep_soil_T_in_K. IF WE VARY THE SURFACE TEMPERATURE
!     sfc_soil_T_in_K SINUSSODALLY IN TIME WE CAN DEDUCE A RELATION BETWEEN THE
!     SURFACE TEMPERATURE, THE SOIL HEAT FLUX (SHF)  AND deep_soil_T_in_K:
!
!              D(sfc_soil_T_in_K)/DT=C1*SHF - C2*(sfc_soil_T_in_K-deep_soil_T_in_K)
!
!     HOWEVER IN REALITY THE TEMPERATURE deep_soil_T_in_K ALSO VARIES IN TIME, IT
!     MAY BE CALCULATED FROM:
!
!              D(deep_soil_T_in_K)/DT= C2*SHF
!
!     THE EQUATIONS GIVEN ABOVE ARE ANALOGOUS TO THOSE USED BY
!     DEARDORFF  (1978).
!
!-----------------------------------------------------------------------
  use stats_precision, only: time_precision
  use constants, only: pi, stefan_boltzmann

  implicit none

  integer, parameter :: LEAF = 2
  ! NEEDED CONSTANsfc_soil_T_in_K
  real, parameter :: CPD=1005.E0 



  ! Input Variables
  real(kind=time_precision), intent(in) :: BTIME, TIME
  real, intent(in) :: DT
  integer, intent(in) :: ITF, ITL
  !real, intent(in) :: R0, FS, FSIN, FLDO, FLUP, WTQ
  real, intent(in) :: R0, FS, FSIN, FLDO, WTQ
  

  ! Input/Output Variables
  real, dimension(LEAF), intent(inout):: deep_soil_T_in_K,sfc_soil_T_in_K,veg_T_in_K

          REAL CS,KS,RS,C1,C2,C3,D1,SHF,SHFS,FLUP

  integer :: IT
  !FUNCTIONALITY OF THE CODE:
  !UPDATE SURFACE AND SOIL TEMP, SOIL HEAT FLUX, WHILE ASSUMING THAT NET RADIATION
  !AND TURBULENT HEAT FLUXES ARE ALREADY AVAILABLE FROM ANOTHER SUBROUTINE.
  !
  !
!     MEANING OF VARIABLES
!	FLDO= downwelling longwaverad
!	FLUP= upwelling shortwaverad
!	SHF= soil heatflux
!	R0(0)= air density at surface
!	KAPPA= Von Karman Constant
!	CPD= air heatcapacity at constant pressure
!	FSIN= net solar radiation
!	FS= net solar radiation array
!	WTQ= wet bulb equavelnt heatflux (is sum of sensible and latent surfaceheat flux)
!	LV = latent heat of vaporization
!	IX= index for horizontal dimension (always 0 here)
!	IZ= index for vertical height
!	TIME = time
!	BTIME= starttime
!	DT= MODEL TIMESTEP (SHOULD BE SMALLER THAN 60 SEC)
!	C1= COEFFICIENT IN FORCE RESTORE 1
!	C2= COEFFICIENT IN FORCE RESTORE 2
!	C3= COEFFICIENT IN FORCE RESTORE 3
!	CS= SOIL HEAT CAPACITY
!	RS= SOIL DENSITY
!	KS= SOIL HEAT DIFFUSIVITY
!	PI = 3.1415927.....
!	ITF= CURRENT TIME STEP
!	ITL= NEXT TIME STEP
!	HVEL= vector wind speed at first model level
!	PES0TE,T0TET,ATET,BTET= Coefficiensfc_soil_T_in_K in Tetens formula
!	G = Gravity acceleration
!	RD = gas constant dry air
!	RV = gas constant water vapor
!	CPL = heat capacity liquid water
!	QWA= intermediate variable for water vapor
!	QW = water vapor array
!	QVS= Saturated water vapor pressure at surface
!	DEQW = difference between first model level and surface
!	FACT= Resistance factor between first level and surface humidity (See Duynkerke 1991)
!	RCAN= canopy resistance (RCAN about 100 sm-1)
!	KGZ= Surface layer exchange coefficient.


!
!       SOIL PARAMETERS
!
  CS=2.00E3
  RS=1.00E3
  KS=2.00E-7
  D1=SQRT(KS*3600.E0*24.E0)
  C1=2.E0*SQRT(PI)/(RS*CS*D1)
  C2=2.E0*PI/(3600.E0*24.E0)
  C3=SQRT(PI*2.E0)/(EXP(PI/4.E0)*RS*CS*SQRT(KS*3600.E0*24.E0* &
                   365.E0))
  
  !INITIALIZE SURFACE VEGETATION TEMP (veg_T_in_K), SOIL SURF. TEMP (sfc_soil_T_in_K), DEEP SOIL TEMP (deep_soil_T_in_K)

  IF (TIME .EQ. BTIME)  THEN
    DO IT=1,LEAF
      veg_T_in_K(IT) = initial_veg_t_sfc
      sfc_soil_T_in_K(IT) = initial_t_sfc
      deep_soil_T_in_K(IT) = initial_deep_t_sfc
      SHF=0.
    END DO
  ELSE
    FLUP = stefan_boltzmann * (veg_T_in_K(ITF)**4)
    !CALCULATE NET RADIATION MINUS TURBULENT HEAT FLUX
    SHFS=FLDO-FLUP-WTQ*R0*CPD+FS
    !CALCULATE SOIL HEAT FLUX
    SHF=3.0*(veg_T_in_K(ITF)-sfc_soil_T_in_K(ITF))+0.05*FSIN
    !SHF = 1.0
    !UPDATE SURF VEG TEMP
    veg_T_in_K(ITL)=veg_T_in_K(ITF)+DT*5.E-5*(SHFS-SHF)
    !UPDATE SOIL TEMP
    sfc_soil_T_in_K(ITL)=sfc_soil_T_in_K(ITF)+DT*( C1*SHF-C2*( sfc_soil_T_in_K(ITF)-deep_soil_T_in_K(ITF) ) )
    !UPDATE DEEP SOIL TEMP
    deep_soil_T_in_K(ITL)=deep_soil_T_in_K(ITF)+DT*C3*SHF
  ENDIF

  print *, "veg_T_in_K = ", veg_T_in_K(ITF)
  print *, "FLDO = ", FLDO
  print *, "FLUP = ", FLUP
  print *, "WTQ*R0*CPD = ", WTQ*R0*CPD
  print *, "FS = ", FS
  print *, "C1 = ", C1
  print *, "SHF= ", SHF

  RETURN
  end subroutine prognose_soil_T_in_K
!------------------------------------------------------------------------------------------
  subroutine initialize_surface( initial_veg_t_sfc_in, initial_t_sfc_in, initial_deep_t_sfc_in )
    implicit none
    
    ! Input Variables
    real, intent(in) :: &
    initial_veg_t_sfc_in, &
    initial_t_sfc_in, &
    initial_deep_t_sfc_in
    
    initial_veg_t_sfc = initial_veg_t_sfc_in
    initial_t_sfc = initial_t_sfc_in
    initial_deep_t_sfc = initial_deep_t_sfc_in

  end subroutine initialize_surface

end module surface
