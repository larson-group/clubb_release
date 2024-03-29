! CVS:  $Id: bugs_swr.F,v 1.2 2005-11-30 22:24:38 dschanen Exp $
! CVS:  $Name: not supported by cvs2svn $

!-----------------------------------------------------------------------
#ifdef USE_BUGSrad_ocast_random

      subroutine bugs_swr
     +                (  ncol ,     nlm ,       pp ,      ppl
     +,                    dp ,      tt ,     rmix ,    cwrho
     +,                 cirho ,   o3mix ,       ts ,     amu0
     +,                   slr ,   alvdf ,    alndf ,    alvdr
     +,                 alndr ,  cldamt ,   cldmax ,    umco2
     +,                 umch4 ,   umn2o ,       b1 ,       b2
     +,                    b3 ,      b4 ,     fdsw ,   fdswcl
     +,                  fusw ,  fuswcl ,   radvbc , radvbccl
     +,                radvdc ,radvdccl ,   radnbc , radnbccl
     +,                radndc ,radndccl ,sel_rules,  c_maximal
     +,                cf_max , cf_rndm
     +                )

#else

      subroutine bugs_swr
     +                (  ncol ,     nlm ,       pp ,      ppl
     +,                    dp ,      tt ,     rmix ,    cwrho
     +,                 cirho ,   o3mix ,       ts ,     amu0
     +,                   slr ,   alvdf ,    alndf ,    alvdr
     +,                 alndr ,  cldamt ,   cldmax ,    umco2
     +,                 umch4 ,   umn2o ,       b1 ,       b2
     +,                    b3 ,      b4 ,     fdsw ,   fdswcl
     +,                  fusw ,  fuswcl ,   radvbc , radvbccl
     +,                radvdc ,radvdccl ,   radnbc , radnbccl
     +,                radndc ,radndccl ,sel_rules
     +                )

#endif

      use kinds
      use gases_ckd, only: gases, stanps, pscale
!      use gases_ckd, only:  stanps, pscale

      implicit none

!-----------------------------------------------------------------------
! REFERENCES:
! bugs_swr replaces crcswr written by G. Stephens. BUGSswr computes the
! downward and upward SW radiative fluxes, and SW heating rates.
! Laura D. Fowler (slikrock/08-20-97).

! send comments to laura@slikrock.atmos.colostate.edu and
! partain@atmos.colostate.edu.

! MODIFICATIONS:
! * moved the computation of the all-sky and clear-sky radiative heating
!   rates to bugs_rad.
!   Laura D. Fowler and Phil Partain(slikrock/01-27-99).

! * added effective radii of cloud droplets and ice crystals that are
!   dependent on the cloud water and cloud ice contents.
!   Laura D. Fowler/slikrock (06-08-00).

! * cleaned up the argument list to remove variables related to short
!   wave radiative transfer.
!   Laura D. Fowler/slikrock (02-01-00).

! * changed declarations to adapt the code from BUGS4 to BUGS5.
!   Laura D. Fowler/slikrock (02-01-00).

! SUBROUTINES CALLED:

!     pscale         : pressure scaling.
!     cloudg          : computes cloud optical properties of water/ice
!                       clouds.
!     rayle           : Computes Rayleigh scattering properties.
!     comscp1         : combines optical properties for gray absorption
!                       (clouds and water vapor continuum).
!     comscp2         : combines optical properties for non-gray gaseous
!                       absorption.
!     gases           : computes gases absorption.
!     two_rt_sw       : two-stream parameterization.

! FUNCTIONS CALLED:
!     none. 

! INCLUDED COMMONS:
!     none.

! ARGUMENT LIST VARIABLES:
! All arrays indexed as nlm correspond to variables defined in the
! middle of layers. All arrays indexed as nlm+1 correspond to variables
! defined at levels at the top and bottom of layers.

!     INPUT ARGUMENTS:
!     ----------------
      logical (kind=log_kind), intent(in)::
     &  sel_rules

      integer (kind=int_kind), intent(in)::
     &  ncol   !Length of sub-domain.
     &, nlm    !Number of layers.

      real (kind=dbl_kind), intent(in)::
     &  umco2  !Concentration of CO2                              (ppm).
     &, umch4  !Concentration of CH4                              (???).
     &, umn2o  !Concentration of N2o                              (???).

      real (kind=dbl_kind), intent(in), dimension(ncol)::
     &  ts!Surface temperature                                      (K).
     &, amu0   !Cosine of solar zenith angle                        (-).
     &, slr    !Fraction of daylight                                (-).
     &, alvdr  !Visible direct surface albedo                       (-).
     &, alndr  !Near-IR direct surface albedo                       (-).
     &, alvdf  !Visible diffuse surface albedo                      (-).
     &, alndf  !Near-IR diffuse surface albedo                      (-).
     &, cldmax !Maximum cloud fraction                              (-).

      real (kind=dbl_kind), intent(in), dimension(ncol,nlm)::
     &  ppl     !Layer pressure                                   (hPa).
     &, dp      !Layer thickness                                  (hPa).
     &, tt      !Temperature                                        (K).
     &, rmix    !Water vapor mixing ratio                       (kg/kg).
     &, cwrho   !Cloud water water content                     (g/m^-3).
     &, cirho   !Cloud ice content                             (g/m^-3).
     &, o3mix   !Ozone mixing ratio                             (kg/kg).
     &, cldamt  !Cloud fraction                                     (-).
     &, b1      !Cloud overlap parameter                            (-).
     &, b2      !Cloud overlap parameter                            (-).
     &, b3      !Cloud overlap parameter                            (-).
     &, b4      !Cloud overlap parameter                            (-).

      real (kind=dbl_kind), intent(in), dimension(ncol,nlm+1)::
     &  pp      !Level pressure                                   (hPa).

#ifdef USE_BUGSrad_ocast_random
      real (kind=dbl_kind), intent(in), dimension(ncol) ::
     & c_maximal !                                                  (-).

      real (kind=dbl_kind), intent(in), dimension(ncol,nlm) ::
     & cf_max    !                                                  (-).
     &,cf_rndm   !                                                  (-).
#endif


!     OUTPUT ARGUMENTS:
!     -----------------
      real (kind=dbl_kind), intent(out), dimension(ncol)::
     &  radvbc  !SFC all-sky visible direct net SW radiation   (W/m^-2).
     &, radvbccl!SFC clear-sky visible direct net SW radiation (W/m^-2).
     &, radvdc  !SFC all-sky visible direct net SW radiation   (W/m^-2).
     &, radvdccl!SFC clear-sky visible direct net SW radiation (W/m^-2).
     &, radnbc  !SFC all-sky near-ir direct net SW radiation   (W/m^-2).
     &, radnbccl!SFC clear-sky near-ir direct net SW radiation (W/m^-2).
     &, radndc  !SFC all-sky near-ir direct net SW radiation   (W/m^-2).
     &, radndccl!SFC clear-sky near-ir direct net SW radiation (W/m^-2).

      real (kind=dbl_kind), intent(out), dimension(ncol,nlm+1)::
     &  fdsw    !Downward SW flux                              (W/m^-2).
     &, fdswcl  !Downward clear-ksy SW flux                    (W/m^-2).
     &, fusw    !Upward SW flux                                (W/m^-2).
     &, fuswcl  !Upward clear-sky SW flux                      (W/m^-2).
     
! LOCAL VARIABLES:  

      integer (kind=int_kind)::
     &  mb      !Total number of spectral intervals.
     &, mbs     !Number of shortwave (SW) spectral intervals.
     &, mbir    !Number of shortwave (LW) spectral intervals.
      parameter(mb=18,mbs=6,mbir=12)        

      integer (kind=int_kind) ::    
     &  i       !Horizontal index.
     &, l       !Vertical index.
     &, ib      !Index of spectral interval.
     &, ig      !Index of k-distribution.

      integer (kind=int_kind), dimension(ncol,nlm)::
     &  ip1     !Used in conjunction with pressure weigthing.
     &, ip2     !Used in conjunction with pressure weigthing.

      real (kind=dbl_kind)
     &  hk      !Weighted spectral solar constant              (W/m^-2).       .
     &, tmax    !Temperature threshold                              (K).
     &, eps     !Threshold for cloud optical properties                .
     &, pdist
      data eps,tmax,pdist /1.e-05,340.,2./     
      
      real (kind=dbl_kind), dimension(mbs):: 
     &  kg      !Nb of k-distributions per spectral intervals.  
      data kg /10,8,12,7,12,5/  

      real (kind=dbl_kind), dimension(mbs)::
     &  asym_wat!Spectral asymmetry factor of water clouds.
     &, asym_ice!Spectral asymmetry factor of ice clouds.
     &, ri      !Coefficients related to Rayleigh absorption.     
      data ri / 0.9022e-5, 0.5282e-6, 0.5722e-7
     &,         0.1433e-7, 0.4526e-8, 0.1529e-8 /

      real (kind=dbl_kind), dimension(mb)::
     &  cnrw    !Real part of refractive index (Water clouds).
     &, cniw    !Imaginary part of refractive index (Water clouds).
     &, cnri    !Real part of refractive index (Ice clouds).
     &, cnii    !Imaginary part of refractive indec (Ice clouds).
     &, xlam    !Center of spectral band.

      real (kind=dbl_kind), dimension(ncol,mbs)::
     &  asdir   !Spectral direct surface albedo                     (-).
     &, asdif   !Spectral diffuse surface albedo                    (-).

      real (kind=dbl_kind), dimension(ncol,nlm)::
     &  rew     !Effective radius for cloud water                  (mu).
     &, rei     !Effective radius for cloud ice                    (mu).
     &, ttem    !Local temperature                                  (K).
     &, pkd     !
     &, tau1    !All-sky optical depth                              (-).
     &, tauclr1 !Clear-sky optical depth                            (-).
     &, tau     !All-sky optical depth                              (-).
     &, tauclr  !Clear-sky optical depth                            (-).
     &, taer    !Aerosol optical depth                              (-).
     &, tray    !Rayley optical depth                               (-).
     &, tg      !Gases optical depth                                (-).
     &, tgm     !WV continuum optical depth                         (-).
     &, tcldi   !Ice cloud optical depth                            (-).
     &, tcldw   !Water cloud optical depth                          (-).
     &, wc      !All-sky single scattering albedo                   (-).
     &, wcclr   !Clear-sky single scattering albedo                 (-).
     &, waer    !Aerosol single scattering albedo                   (-).
     &, wray    !Rayley single scattering albedo                    (-).
     &, wcldi   !Ice cloud single scattering albedo                 (-).
     &, wcldw   !Water cloud single scattering albedo               (-).
     &, asym    !All-sky asymmetry factor                           (-).
     &, asyclr  !Clear-sky asymmetry factor                         (-).
     &, asyaer  !Aerosol asymmetry factor                           (-).
     &, asycldi !Ice cloud asymmetry factor                         (-).
     &, asycldw !Water cloud asymmetry factor                       (-).
     &, fwclr   !
     &, fwcld   !

      real (kind=dbl_kind), dimension(ncol,nlm+1)::
     &  fdgdir!Spectral direct downward flux                    (W/m^2).
     &, fdgcldir!Spectral direct clear-sky downward flux        (W/m^2).
     &, fdgdif!Spectral diffuse downward flux                   (W/m^2).
     &, fdgcldif!Spectral diffuse clear-sky downward flux       (W/m^2).
     &, fugdif!Spectral diffuse upward flux                     (W/m^2).
     &, fugcldif!Spectral diffuse clear-sky upward flux         (W/m^2).

!     shortwave asymmetry parameters:
!     (assumes: re=10 for water; re=30 for ice)
      data asym_wat / 0.8625, 0.8469, 0.8287, 0.8182, 0.9472, 0.7630 /
      data asym_ice / 0.8678, 0.8640, 0.8653, 0.8615, 0.9526, 0.8293 /

!---  cnrw and cniw (water clouds):
      data cnrw/1.3422,1.3281,1.3174,1.2901,1.3348,1.3700,1.3191,1.2821
     &,         1.3160,1.3030,1.2739,1.2319,1.1526,1.1981,1.3542,1.4917
     &,         1.5463,1.8718/
      data cniw/6.4790e-9,1.3417e-06,1.2521e-4,7.1533e-4,4.2669e-2
     &,         4.3785e-3,1.3239e-2 ,1.5536e-2,5.3894e-2,3.4346e-2
     &,         3.7490e-2,4.7442e-2 ,1.2059e-1,3.3546e-1,4.1698e-1
     &,         4.0674e-1,3.6362e-1 ,5.2930e-1/

!--- cnri and cnii (ice clouds):
      data cnri/1.3266,1.2986,1.2826,1.2556,1.2963,1.3956
     &,         1.3324,1.2960,1.3121,1.3126,1.2903,1.2295
     &,         1.1803,1.5224,1.5572,1.5198,1.4993,1.7026/
      data cnii/7.0696e-9,9.1220e-7,1.2189e-4,5.7648e-4,4.3144e-2
     &,         8.2935e-3,1.5540e-2,2.5594e-2,5.9424e-2,5.1511e-2
     &,         4.0325e-2,4.7994e-2,2.3834e-1,3.0697e-1,1.1852e-1
     &,         4.3048e-2,6.3218e-2,1.5843e-1/

!---- spectral band center:
      data xlam/0.45  ,1.0   ,1.6  ,2.2  ,3.0   ,3.75  ,4.878 ,5.556
     &,         6.452 ,7.547 ,8.511,9.615,11.236,13.605,16.529,21.277 
     &,         29.412,71.403/

!-----------------------------------------------------------------------

!---- 0. initialize local and output arrays:

      radvbc(:)   = 0.
      radvbccl(:) = 0.
      radvdc(:)   = 0.
      radvdccl(:) = 0.
      radnbc(:)   = 0.
      radnbccl(:) = 0.
      radndc(:)   = 0.
      radndccl(:) = 0.

      fdsw(:,:)   = 0.
      fdswcl(:,:) = 0.
      fusw(:,:)   = 0.
      fuswcl(:,:) = 0.

      rew(:,:)  = 10.
      rei(:,:)  = 30.

      fdgdir(:,:) = 0.0
      fdgcldir(:,:) = 0.0
      fdgdif(:,:) = 0.0
      fdgcldif(:,:) = 0.0
      fugdif(:,:) = 0.0
      fugcldif(:,:) = 0.0

      do l = 1, nlm
         do i = 1, ncol
            ttem(i,l) = min(tmax,tt(i,l))
          enddo
      enddo

      do i = 1, ncol
        asdir(i,1)   = alvdr(i)
        asdir(i,2:6) = alndr(i)
        asdif(i,1)   = alvdf(i)
        asdif(i,2:6) = alndf(i)
      enddo

!--   pressure scaling:

       call pscale(ncol,nlm,ppl,stanps,pkd,ip1,ip2)

!---- 1. loop over the mbs spectral intervals starts here:

      do ib = 1, mbs         
      
         tgm(:,:)      = 0.
         taer(:,:)     = 0.
         waer(:,:)     = 0.
         asyaer(:,:)   = 1.
!        fdswband(:,:) = 0.
!        fuswband(:,:) = 0.

!---- 1.1 rayleigh absorption:

         call rayle
     +            (ncol ,  nlm , mbs ,  ib
     +,            amu0 ,   ri ,  pp , ppl
     +,            tray , wray
     +            )

!---- 1.2 optical properties of water and ice clouds (as in crcswr for
!        now):

         call cloudg
     +           (   ncol ,     nlm ,    mb ,    ib
     +,                pp ,      tt , cwrho ,   rew
     +,             pdist ,    cnrw ,  cniw ,  cnri
     +,              cnii ,    xlam , tcldw , wcldw
     +,           asycldw , .false.
     +           )

         call cloudg
     +           (   ncol ,   nlm   ,    mb ,    ib
     +,                pp ,    tt   , cirho ,   rei
     +,             pdist ,  cnrw   ,  cniw ,  cnri
     +,              cnii ,  xlam   , tcldi , wcldi
     +,           asycldi , .true.
     +           )

!     the asymmetry factor for water and ice clouds are fixed as
!     functions of the spectral intervals.

         do l = 1, nlm
            do i = 1, ncol
               if(cwrho(i,l).ge.eps) asycldw(i,l) = asym_wat(ib)
               if(cirho(i,l).ge.eps) asycldi(i,l) = asym_ice(ib)
            enddo
         enddo

!---- 1.3 combines single-scattering properties for gray absorption:

        call comscp1
     +          (   ncol ,     nlm ,  taer ,   tcldi
     +,            tcldw ,     tgm ,  tray ,    waer
     +,            wcldi ,   wcldw ,  wray ,  asyaer
     +,          asycldi , asycldw ,  tau1 , tauclr1
     +,             asym ,  asyclr , fwcld ,   fwclr
     +          )

!---- loop over the k-probability distributions starts here:

         do ig = 1, int( kg(ib) ) ! D.S. for g95 

!---- 1.4 non-gray gaseous absorption:         

            call gases
     +              ( ncol ,   nlm ,    ib ,    ig
     +,                 pp ,    dp ,  ttem ,  rmix
     +,              o3mix , umco2 , umch4 , umn2o
     +,                 hk ,    tg ,   pkd ,   ip1
     +,                ip2
     +              )

!---- 1.5 combines single-scattering properties:

            call comscp2
     +              (  ncol ,  nlm ,      tg , fwcld
     +,               fwclr , tau1 , tauclr1 ,   tau
     +,              tauclr ,   wc ,   wcclr
     +              )

!---- 1.6 two-stream approximation:
#ifdef nooverlap
! No overlap
            call two_rt_sw
     +              (  ncol ,    nlm ,       mbs ,     ib
     +,                 slr ,   amu0 ,        wc ,   asym
     +,                 tau ,  asdir ,     asdif , fugdif
     +,              fdgdir , fdgdif , sel_rules
     +              )

            call two_rt_sw
     +              (    ncol ,     nlm ,       mbs ,       ib
     +,                   slr ,    amu0 ,     wcclr ,   asyclr
     +,                tauclr ,   asdir ,     asdif , fugcldif
     +,              fdgcldir ,fdgcldif , sel_rules
     +              )
#else
#ifdef USE_BUGSrad_ocast_random
!New overcast-random method
            call two_rt_sw_gsolap
     +              (
     +                   ncol ,     nlm ,      mbs ,        ib
     +,                   slr ,    amu0 ,    wcclr ,    asyclr
     +,                tauclr ,      wc ,     asym ,       tau
     +,                 asdir ,   asdif , fdgcldir ,  fdgcldif
     +,              fugcldif ,  fdgdir ,   fdgdif ,    fugdif
     +,                cldamt , c_maximal,  cf_max ,   cf_rndm
     +              )
#else
!Standard max/random
            call two_rt_sw
     +              (
     +                   ncol ,     nlm ,       mbs ,       ib
     +,                   slr ,    amu0 ,     wcclr ,   asyclr
     +,                tauclr ,   asdir ,     asdif , fugcldif
     +,              fdgcldir ,fdgcldif , sel_rules
     +              )

            call  two_rt_sw_bs
     +              (
     +                   ncol ,      nlm ,       mbs ,       ib
     +,                   slr ,     amu0 ,        wc ,    wcclr
     +,                  asym ,   asyclr ,       tau ,   tauclr
     +,                 asdir ,    asdif ,    fugdif ,   fdgdir
     +,                fdgdif ,sel_rules ,        b1 ,       b2
     +,                    b3 ,       b4
     +              )
#endif
#endif

            fdsw(:,:)     = fdsw(:,:)   
     +                    + (fdgdir(:,:)+fdgdif(:,:))*hk
            fusw(:,:)     = fusw(:,:)   + fugdif(:,:)*hk
            fdswcl(:,:)   = fdswcl(:,:) 
     +                    + (fdgcldir(:,:)+fdgcldif(:,:)) * hk
            fuswcl(:,:)   = fuswcl(:,:) + fugcldif(:,:)*hk

!---- 1.7 computes the surface visible and near infrared net radiation.
            select case (ib)

               case(1)
                  radvbc(:)   = radvbc(:) + fdgdir(:,nlm+1)*hk
                  radvbccl(:) = radvbccl(:) + fdgcldir(:,nlm+1)*hk
                  radvdc(:)   = radvdc(:) + fdgdif(:,nlm+1)*hk
                  radvdccl(:) = radvdccl(:) + fdgcldif(:,nlm+1)*hk

               case(2:6)
                  radnbc(:) = radnbc(:) + fdgdir(:,nlm+1)*hk
                  radnbccl(:) = radnbccl(:) + fdgcldir(:,nlm+1)*hk
                  radndc(:) = radndc(:) + fdgdif(:,nlm+1)*hk
                  radndccl(:) = radndccl(:) + fdgcldif(:,nlm+1)*hk

            end select 

         enddo ! end k-distribution

!           print *, 'ib= ', ib, 'fdsw= ', fdsw
!       print *, 'ib= ', ib, 'wc_swr= ', wc
!       print *, 'ib= ', ib, 'asym_swr= ', asym
      enddo ! end spectral interval

      return

      end subroutine bugs_swr

!-----------------------------------------------------------------------

