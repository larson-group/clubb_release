!$Id$
module surface

implicit none

public :: compute_surface_temp, initialize_surface

real, private :: initial_veg_t_sfc, initial_t_sfc, initial_deep_t_sfc

private


contains
!----------------------------------------------------------------------
SUBROUTINE compute_surface_temp( BTIME, TIME, DT, ITF, ITL, R0, FS, FSIN, FLDO, FLUP, WTQ, &
                   TVE, TS, TD )
!
!
!     IN THIS SUBROUTINE TQ AND QW AT THE SURFACE ARE CALCULATED AS A
!     FUNCTION OF THE SURFACE TEMPERATURE (TS). TS IS CALCULATED  FROM
!     THE SURFACE ENERGY BUDGET.
!
! ************************************** 2
!               *                   *
!               *                   *
!             --*         +         *--  1
!               *                   *
! *** SURFACE ************+************* 0
!                        TS              L (EVEL)
!
!
!     HEAT CONDUCTION (IN A HOMOGENEOUS MEDIUM) CAN BE DESCRIBED BY
!     THE EQUATION:
!
!              D(TS)/DT=KS*D(D(TS)/DZ)/DZ
!
!     IN WHICH KS IS THE SOIL THERMAL DIFFUSITY.
!     WE CONSIDER A SEMI HALF INFINITE MEDIUM, INITIALLY AT THE
!     CONSTANT TEMPERATURE TD. IF WE VARY THE SURFACE TEMPERATURE
!     TS SINUSSODALLY IN TIME WE CAN DEDUCE A RELATION BETWEEN THE
!     SURFACE TEMPERATURE, THE SOIL HEAT FLUX (SHF)  AND TD:
!
!              D(TS)/DT=C1*SHF - C2*(TS-TD)
!
!     HOWEVER IN REALITY THE TEMPERATURE TD ALSO VARIES IN TIME, IT
!     MAY BE CALCULATED FROM:
!
!              D(TD)/DT= C2*SHF
!
!     THE EQUATIONS GIVEN ABOVE ARE ANALOGOUS TO THOSE USED BY
!     DEARDORFF  (1978).
!
!-----------------------------------------------------------------------
  use stats_precision, only: time_precision
  use constants, only: pi

  implicit none

  integer, parameter :: LEAF = 2
  ! NEEDED CONSTANTS
  real, parameter :: CPD=1005.E0 



  ! Input Variables
  real(kind=time_precision), intent(in) :: BTIME, TIME
  real, intent(in) :: DT
  integer, intent(in) :: ITF, ITL
  real, intent(in) :: R0, FS, FSIN, FLDO, FLUP, WTQ
  

  ! Input/Output Variables
  real, dimension(LEAF), intent(inout):: TD,TS,TVE

          REAL CS,KS,RS,C1,C2,C3,D1,SHF,SHFS

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
!	PES0TE,T0TET,ATET,BTET= Coefficients in Tetens formula
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
  
  !INITIALIZE SURFACE VEGETATION TEMP (TVE), SOIL SURF. TEMP (TS), DEEP SOIL TEMP (TD)

  IF (TIME .EQ. BTIME)  THEN
    DO IT=1,LEAF
      TVE(IT) = initial_veg_t_sfc
      TS(IT) = initial_t_sfc
      TD(IT) = initial_deep_t_sfc
      SHF=0.
    END DO
  ELSE
    !CALCULATE NET RADIATION MINUS TURBULENT HEAT FLUX
    SHFS=FLDO-FLUP-WTQ*R0*CPD+FS
    !CALCULATE SOIL HEAT FLUX
    SHF=3.0*(TVE(ITF)-TS(ITF))+0.05*FSIN
    !SHF = 1.0
    !UPDATE SURF VEG TEMP
    TVE(ITL)=TVE(ITF)+DT*5.E-5*(SHFS-SHF)
    !UPDATE SOIL TEMP
    TS(ITL)=TS(ITF)+DT*( C1*SHF-C2*( TS(ITF)-TD(ITF) ) )
    !UPDATE DEEP SOIL TEMP
    TD(ITL)=TD(ITF)+DT*C3*SHF
  ENDIF

  !print *, "TD = ", TD(ITL)
  !print *, "C2 = ", C2
  !print *, "C3 = ", C3

  RETURN
  end subroutine compute_surface_temp
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
