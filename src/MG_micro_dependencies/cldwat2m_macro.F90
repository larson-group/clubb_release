
  module cldwat2m_macro

  !--------------------------------------------------- !
  ! Purpose     : CAM Interface for Cloud Macrophysics !
  ! Author      : Sungsu Park                          !
  ! Description : Park et al. 2010.                    !
  ! For questions, contact Sungsu Park                 !
  !                        e-mail : sungsup@ucar.edu   !
  !                        phone  : 303-497-1375       !
  !--------------------------------------------------- !

   use shr_kind_mod,     only: r8=>shr_kind_r8
   use spmd_utils,       only: masterproc
   use ppgrid,           only: pcols, pver, pverp
   use abortutils,       only: endrun
   use wv_saturation,    only: estblf, hlatv, tmin, hlatf, rgasv, pcf, cp, epsqs, ttrice, & 
                               vqsatd2_water, vqsatd2_water_single, polysvp
   use cam_history,      only: addfld, add_default, phys_decomp, outfld 
   use cam_logfile,      only: iulog

   implicit none
   private
   save

   public :: aist_vector

   ! -------------- !
   ! Set Parameters !
   ! -------------- !

   ! -------------------------- !
   ! Parameters for Ice Stratus !
   ! -------------------------- !

   integer,  parameter          :: iceopt       = 5             ! Option for ice cloud closure. 5 : Modified Slingo Formula. For other options, modify each 2 subroutine.
   real(r8), public,  parameter :: rhmini       = 0.80_r8       ! Minimum rh for ice cloud fraction > 0.
   real(r8), public,  parameter :: rhmaxi       = 1.1_r8        ! rhi at which ice cloud fraction = 1.
   real(r8), parameter          :: qist_min     = 1.e-7_r8      ! Minimum in-stratus ice IWC constraint [ kg/kg ]
   real(r8), parameter          :: qist_max     = 5.e-3_r8      ! Maximum in-stratus ice IWC constraint [ kg/kg ]

   ! ----------------------------- !
   ! Parameters for Liquid Stratus !
   ! ----------------------------- !

   logical,  parameter          :: CAMstfrac    = .false.       ! If .true. (.false.), use Slingo (triangular PDF-based) liquid stratus fraction
   logical,  parameter          :: freeze_dry   = .false.       ! If .true., use 'freeze dry' in liquid stratus fraction formula
   real(r8), parameter          :: qlst_min     = 2.e-5_r8      ! Minimum in-stratus LWC constraint [ kg/kg ]
   real(r8), parameter          :: qlst_max     = 3.e-3_r8      ! Maximum in-stratus LWC constraint [ kg/kg ]
   real(r8), parameter          :: cc           = 0.1_r8        ! For newly formed/dissipated in-stratus CWC ( 0 <= cc <= 1 )
   real(r8), parameter          :: premib       = 700.e2_r8     ! Bottom height for mid-level liquid stratus fraction
   integer,  parameter          :: niter        = 2             ! For iterative computation of QQ with 'ramda' below.
   real(r8), parameter          :: ramda        = 0.5_r8        ! Explicit : ramda = 0, Implicit : ramda = 1 ( 0<= ramda <= 1 )
   real(r8), private            :: rhminl                       ! Critical RH for low-level  liquid stratus clouds
   real(r8), private            :: rhminh                       ! Critical RH for high-level liquid stratus clouds
   real(r8), private            :: premit                       ! Top    height for mid-level liquid stratus fraction

   contains

   subroutine aist_vector( qv_in, T_in, p_in, qi_in, landfrac_in, snowh_in, aist_out, ncol )

   ! --------------------------------------------------------- !
   ! Compute non-physical ice stratus fraction                 ! 
   ! --------------------------------------------------------- !

   use physconst,     only: rair
   use wv_saturation, only: vqsatd2_water

   implicit none
  
   integer,  intent(in)  :: ncol 
   real(r8), intent(in)  :: qv_in(pcols)       ! Grid-mean water vapor[kg/kg]
   real(r8), intent(in)  :: T_in(pcols)        ! Temperature
   real(r8), intent(in)  :: p_in(pcols)        ! Pressure [Pa]
   real(r8), intent(in)  :: qi_in(pcols)       ! Grid-mean ice water content [kg/kg]
   real(r8), intent(in)  :: landfrac_in(pcols) ! Land fraction
   real(r8), intent(in)  :: snowh_in(pcols)    ! Snow depth (liquid water equivalent)
   real(r8), intent(out) :: aist_out(pcols)    ! Non-physical ice stratus fraction ( 0<= aist <= 1 )

   ! Local variables

   real(r8) qv                              ! Grid-mean water vapor[kg/kg]
   real(r8) T                               ! Temperature
   real(r8) p                               ! Pressure [Pa]
   real(r8) qi                              ! Grid-mean ice water content [kg/kg]
   real(r8) landfrac                        ! Land fraction
   real(r8) snowh                           ! Snow depth (liquid water equivalent)
   real(r8) aist                            ! Non-physical ice stratus fraction ( 0<= aist <= 1 )

   real(r8) rhmin                           ! Critical RH
   real(r8) rhwght

   real(r8) a,b,c,as,bs,cs                  ! Fit parameters
   real(r8) Kc                              ! Constant for ice cloud calc (wood & field)
   real(r8) ttmp                            ! Limited temperature
   real(r8) icicval                         ! Empirical IWC value [ kg/kg ]
   real(r8) rho                             ! Local air density
   real(r8) esl                             ! Liq sat vapor pressure
   real(r8) esi                             ! Ice sat vapor pressure
   real(r8) ncf,phi                         ! Wilson and Ballard parameters
   real(r8) esat, qsat, dqsdT
   real(r8) esat_in(pcols)
   real(r8) qsat_in(pcols)
   real(r8) dqsdT_in(pcols)

   real(r8) rhi                             ! grid box averaged relative humidity over ice
   real(r8) minice                          ! minimum grid box avg ice for having a 'cloud'
   real(r8) mincld                          ! minimum ice cloud fraction threshold
   real(r8) icimr                           ! in cloud ice mixing ratio
 ! real(r8) qist_min                        ! minimum in cloud ice mixing ratio
 ! real(r8) qist_max                        ! maximum in cloud ice mixing ratio                
   real(r8) rhdif                           ! working variable for slingo scheme

   integer i
 ! integer iceopt                           ! option for ice cloud closure 
                                            ! 1=wang & sassen 2=schiller (iciwc)  
                                            ! 3=wood & field, 4=Wilson (based on smith)
                                            ! 5=modified slingo (ssat & empyt cloud)

   real(r8) icecrit                         ! Critical RH for ice clouds in Wilson & Ballard closure ( smaller = more ice clouds )

   ! Statement functions
   logical land
   land(i) = nint(landfrac_in(i)) == 1

   ! --------- !
   ! Constants !
   ! --------- !

   ! iceopt = 5

   ! Wang and Sassen IWC paramters ( Option.1 )
     a = 26.87_r8
     b = 0.569_r8
     c = 0.002892_r8
   ! Schiller parameters ( Option.2 )
     as = -68.4202_r8
     bs = 0.983917_r8
     cs = 2.81795_r8
   ! Wood and Field parameters ( Option.3 )
     Kc = 75._r8
   ! Wilson & Ballard closure ( Option.4. smaller = more ice clouds)
     icecrit = 0.93_r8
   ! Slingo modified (option 5)
     minice = 1.e-12_r8
     mincld = 1.e-4_r8
   ! qist_min = 1.e-7_r8
   ! qist_max = 5.e-3_r8

   ! ---------------- !
   ! Main computation !
   ! ---------------- !

     aist_out(:) = 0._r8
     esat_in(:)  = 0._r8
     qsat_in(:)  = 0._r8
     dqsdT_in(:) = 0._r8

     call vqsatd2_water(T_in(1:ncol),p_in(1:ncol),esat_in(1:ncol),qsat_in(1:ncol),dqsdT_in(1:ncol),ncol)
     
     do i = 1, ncol

     landfrac = landfrac_in(i)     
     snowh = snowh_in(i)   
     T = T_in(i)
     qv = qv_in(i)
     p = p_in(i)
     qi = qi_in(i)
     qsat = qsat_in(i)
     esl = polysvp(T,0)
     esi = polysvp(T,1)
          
     if( iceopt.lt.3 ) then
         if( iceopt.eq.1 ) then
             ttmp = max(195._r8,min(T,253._r8)) - 273.16_r8
             icicval = a + b * ttmp + c * ttmp**2._r8
             rho = p/(rair*T)
             icicval = icicval * 1.e-6_r8 / rho 
         else
             ttmp = max(190._r8,min(T,273.16_r8))
             icicval = 10._r8 **(as * bs**ttmp + cs)
             icicval = icicval * 1.e-6_r8 * 18._r8 / 28.97_r8
         endif
         aist =  max(0._r8,min(qi/icicval,1._r8)) 
     elseif( iceopt.eq.3 ) then
         aist = 1._r8 - exp(-Kc*qi/(qsat*(esi/esl)))
         aist = max(0._r8,min(aist,1._r8))
     elseif( iceopt.eq.4) then
         if( p .ge. premib ) then
             if( land(i) .and. (snowh.le.0.000001_r8) ) then
                 rhmin = rhminl - 0.10_r8
             else
                 rhmin = rhminl
             endif
         elseif( p .lt. premit ) then
             rhmin = rhminh
         else
             rhwght = (premib-(max(p,premit)))/(premib-premit)
           ! if( land(i) .and. (snowh.le.0.000001_r8) ) then
           !     rhmin = rhminh*rhwght + (rhminl - 0.10_r8)*(1.0_r8-rhwght)
           ! else
                 rhmin = rhminh*rhwght + rhminl*(1.0_r8-rhwght)
           ! endif
         endif
         ncf = qi/((1._r8 - icecrit)*qsat)
         if( ncf.le.0._r8 ) then 
             aist = 0._r8
         elseif( ncf.gt.0._r8 .and. ncf.le.1._r8/6._r8 ) then 
             aist = 0.5_r8*(6._r8 * ncf)**(2._r8/3._r8)
         elseif( ncf.gt.1._r8/6._r8 .and. ncf.lt.1._r8 ) then
             phi = (acos(3._r8*(1._r8-ncf)/2._r8**(3._r8/2._r8))+4._r8*3.1415927_r8)/3._r8
             aist = (1._r8 - 4._r8 * cos(phi) * cos(phi))
         else
             aist = 1._r8
         endif
             aist = max(0._r8,min(aist,1._r8))
     elseif (iceopt.eq.5) then 
! set rh ice cloud fraction
             rhi= (qv+qi)/qsat * (esl/esi)
             rhdif= (rhi-rhmini) / (rhmaxi - rhmini)
             aist = min(1.0_r8, max(rhdif,0._r8)**2)

! limiter to remove empty cloud and ice with no cloud
! and set icecld fraction to mincld if ice exists

             if (qi.lt.minice) then
                aist=0._r8
             else
                aist=max(mincld,aist)
             endif

! enforce limits on icimr
             if (qi.ge.minice) then
                icimr=qi/aist

!minimum
                if (icimr.lt.qist_min) then
                   aist = max(0._r8,min(1._r8,qi/qist_min))
                endif
!maximum
                if (icimr.gt.qist_max) then
                   aist = max(0._r8,min(1._r8,qi/qist_max))
                endif

             endif
     endif 

   ! 0.999_r8 is added to prevent infinite 'ql_st' at the end of instratus_condensate
   ! computed after updating 'qi_st'.  

     aist = max(0._r8,min(aist,0.999_r8))

     aist_out(i) = aist

     enddo

   return
   end subroutine aist_vector

   ! ----------------- !
   ! End of subroutine !
   ! ----------------- !

end module cldwat2m_macro

