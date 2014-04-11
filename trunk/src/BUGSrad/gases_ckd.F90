! CVS:  $Id$
! CVS:  $Name: not supported by cvs2svn $

module gases_ckd
!Correlated-k coefficients and routines for calculating gaseous absorption per Fu, 1991
!The hk* variables are the weights for the k-distribution intervals for each band
!The variable fk1o3 is the ozone absorption coefficient in (atm-cm)^-1
!The various c2h2o, c3h2o, c10ch4, c10n2o etc. variables are coefficients
!required to calculate the absorption coefficient in (atm-cm)^-1.

   use kinds, only: int_kind, dbl_kind
   implicit none
   private

   public :: &
    gases          & ! Subroutine for calculation gas transmission
   ,pscale

!This include file contains all the correlated-k coefficient
!values and associated parameters, (including some public ones)
#include "gases_ckd_data.h"

contains
   
subroutine gases  ( ncol ,   nlm ,    ib ,    ig   &
     ,                pp ,    dp ,    tt ,  rmix   &
     ,             o3mix , umco2 , umch4 , umn2o   &
     ,                hk ,    tg ,   pkd ,   ip1   &
     ,               ip2                           &
                  )

!-----------------------------------------------------------------------
! Calculates the non-gray spectral absorption due to H2O, O3, CO2,
! CHA4, and N2O at long and short wavelengths.

   ! INPUT ARGUMENTS:
   integer (kind=int_kind), intent(in)::   &
    ncol        &! Length of sub-domain.
   ,nlm         &! Number of layers.
   ,ib          &! Index of spectral interval.
   ,ig           ! Index of k-distribution.

   integer (kind=int_kind), intent(in), dimension(ncol,nlm):: &
    ip1         &!Used in conjunction with pressure weighting.
   ,ip2          !Used in conjunction with pressure weighting.

   real (kind=dbl_kind), intent(in):: &
    umco2       &!Concentration of CO2                          (ppm).
   ,umch4       &!Concentration of CH4                          (???).
   ,umn2o        !Concentration of N2O                          (???).
  
   real (kind=dbl_kind), intent(in), dimension(ncol,nlm):: &
    dp          &!Layer Thickness                               (hPa).
   ,tt          &!Layer Temperature                               (K).
   ,rmix        &!Layer Water Vapor Mixing Ratio              (kg/kg).
   ,o3mix       &!Layer Ozone Mixing Ratio                    (kg/kg).    
   ,pkd          !Used in conjunction with pressure weighting     (-).
  
   real (kind=dbl_kind), intent(in), dimension(ncol,nlm+1):: &
    pp           !Level Pressure                                (hPa).   
  
   ! OUTPUT ARGUMENTS:
   real (kind=dbl_kind), intent(out):: &
    hk           !Weight                                          (-).
     
   real (kind=dbl_kind), intent(out), dimension(ncol,nlm):: &
    tg           !Transmission                                    (-).
    
   ! LOCAL VARIABLES:

   real (kind=dbl_kind), dimension(ncol,nlm):: &
    fkg , fkga , fkgb , pq &
   ,tg1 ,  tg2 ,  tg3

   integer :: i ! Added by dschanen UWM for bug fix.
     
!-----------------------------------------------------------------------

   select case (ib)
   case (1)
      ! (50000-14500 cm^-1 ):  nongray absorption by O3.
      ! 619.618 is the solar energy contained in that band (Wm^-2).
      call qopo3s(fk1o3(ig),dp,o3mix,tg)
      hk = 619.618*hk1(ig)

   case (2)
      ! (14500-7700 cm^-1):  nongray absorption by H2O.
      ! 484.295 is the solar energy contained in the band (Wm^-2).
      call qk(ncol,nlm,c2h2o(:,:,ig),tt,fkg,pkd,ip1,ip2)
      call qoph2o(fkg,dp,rmix,tg)
      hk = 484.295*hk2(ig)

   case (3)
      ! (7700-5250 cm^-1):  nongray absorption by H2O.
      ! 149.845 is the solar energy contained in the band (Wm^-2).
      call qk(ncol,nlm,c3h2o(:,:,ig),tt,fkg,pkd,ip1,ip2)
      call qoph2o(fkg,dp,rmix,tg)
      hk = 149.845*hk3(ig)

   case (4)
      ! (5250-4000 cm^-1):  nongray absorption by H2O.
      ! 48.7302 is the solar energy contained in the band (Wm^-2).	
      call qk(ncol,nlm,c4h2o(:,:,ig),tt,fkg,pkd,ip1,ip2)
      call qoph2o(fkg,dp,rmix,tg)      
      hk = 48.7302*hk4(ig)

   case (5)
      ! (4000-2850 cm^-1):  nongray absorption by H2O.
      ! 31.6576 is the solar energy contained in the band (Wm^-2).
      call qk(ncol,nlm,c5h2o(:,:,ig),tt,fkg,pkd,ip1,ip2)
      call qoph2o(fkg,dp,rmix,tg)      
      hk = 31.6576*hk5(ig)

   case (6)
      ! (2850-2500 cm^-1):  nongray absorption by H2O.
      ! 5.79927 is the solar energy contained in the band (Wm^-2).	
      call qk(ncol,nlm,c6h2o(:,:,ig),tt,fkg,pkd,ip1,ip2)
      call qoph2o(fkg,dp,rmix,tg)      
      hk = 5.79927*hk6(ig)

   case (7)
      ! (2200-1900 cm^-1):  nongray absorption by H2O.
      call qk(ncol,nlm,c7h2o(:,:,ig),tt,fkg,pkd,ip1,ip2)
      call qoph2o(fkg,dp,rmix,tg)
      hk = hk7(ig)

   case (8)      
      ! (1900-1700 cm^-1):  nongray absorption by H2O.
      call qk(ncol,nlm,c8h2o(:,:,ig),tt,fkg,pkd,ip1,ip2)
      call qoph2o(fkg,dp,rmix,tg)      
      hk = hk8(ig)

   case (9)
      ! (1700-1400 cm^-1):  nongray absorption by H2O.
      call qk(ncol,nlm,c9h2o(:,:,ig),tt,fkg,pkd,ip1,ip2)
      call qoph2o(fkg,dp,rmix,tg)      
      hk = hk9(ig)     

   case (10)
      ! (1400-1250 cm^-1):  overlapping absorption by H2O, CH4, and N2O
      ! using approach one of Fu(1991).	
      call qk(ncol,nlm,c10h2o(:,:,ig),tt,fkg,pkd,ip1,ip2)
      call qoph2o(fkg,dp,rmix,tg1)      
      call qk(ncol,nlm,c10ch4,tt,fkg,pkd,ip1,ip2)
      call qopch4(fkg,dp,tg2)
      call qk(ncol,nlm,c10n2o,tt,fkg,pkd,ip1,ip2)
      call qopn2o(fkg,dp,tg3)

      tg = tg1 + tg2/1.6*umch4 + tg3/0.28*umn2o
      hk = hk10(ig)     

   case (11)
      ! (1250-1100 cm^-1):  overlapping absorption by H2O, CH4, and N2O
      ! using approach one of Fu(1991).      
      call qk(ncol,nlm,c11h2o(:,:,ig),tt,fkg,pkd,ip1,ip2) 
      call qoph2o(fkg,dp,rmix,tg1)        
      call qk(ncol,nlm,c11ch4,tt,fkg,pkd,ip1,ip2)
      call qopch4(fkg,dp,tg2)      
      call qk(ncol,nlm,c11n2o,tt,fkg,pkd,ip1,ip2)
      call qopn2o(fkg,dp,tg3)
      
      tg = tg1 + tg2/1.6*umch4 + tg3/0.28*umn2o
      hk = hk11(ig)

   case(12)
      ! (1100-980 cm^-1):  overlapping absorption by H2O and O3
      ! using approach one of Fu(1991).	
      call qkio3(ncol,nlm,c12o3(:,:,ig),tt,fkg,pkd,ip1,ip2)
      call qopo3i(fkg,dp,o3mix,tg1)
      call qk(ncol,nlm,c12h2o,tt,fkg,pkd,ip1,ip2)
      call qoph2o(fkg,dp,rmix,tg2)      
      
      tg = tg1 + tg2
      hk = hk12(ig)

   case (13)
      ! (980-800 cm^-1 ):  nongray absorption by H2O.
      call qk(ncol,nlm,c13h2o(:,:,ig),tt,fkg,pkd,ip1,ip2)
      call qoph2o(fkg,dp,rmix,tg)       
      hk = hk13(ig)

   case (14)
      ! (800-670 cm^-1):  overlapping absorption by H2O and CO2
      ! using approach two of Fu(1991).

      forall (i=1:ncol) ! Forall added by Dave Schanen UWM for a bug fix
         where (pp(i,:) .ge. 63.1)
            pq(i,:) = rmix(i,:)
         else where
            pq(i,:) = 0._dbl_kind
         end where
      end forall

      call qk(ncol,nlm,c14hca(:,:,ig),tt,fkga,pkd,ip1,ip2)
      call qk(ncol,nlm,c14hcb(:,:,ig),tt,fkgb,pkd,ip1,ip2)
      
      fkg = fkga/330.0*umco2+pq*fkgb
      call qophc(fkg,dp,tg)
      hk = hk14(ig)

   case (15)
      ! (670-540 cm^-1):  overlapping absorption by H2O and CO2
      ! using approach two of Fu(1991).

      forall (i=1:ncol) ! Forall added by Dave Schanen UWM for a bug fix
         where (pp(i,:) .ge. 63.1)
            pq(i,:) = rmix(i,:)
         else where
            pq(i,:) = 0._dbl_kind
         end where
      end forall

      call qk(ncol,nlm,c15hca(:,:,ig),tt,fkga,pkd,ip1,ip2)
      call qk(ncol,nlm,c15hcb(:,:,ig),tt,fkgb,pkd,ip1,ip2)
      
      fkg = fkga/330.0*umco2+pq*fkgb
      call qophc(fkg,dp,tg)
      hk = hk15(ig)

   case (16)
      ! (540-400 cm^-1 ):  nongray absorption by H2O.
      call qk(ncol,nlm,c16h2o(:,:,ig),tt,fkg,pkd,ip1,ip2)
      call qoph2o(fkg,dp,rmix,tg)      
      hk = hk16(ig)

   case (17)
      ! (400-280 cm^-1 ):  nongray absorption by H2O.
      call qk(ncol,nlm,c17h2o(:,:,ig),tt,fkg,pkd,ip1,ip2)
      call qoph2o(fkg,dp,rmix,tg)         
      hk = hk17(ig)

   case (18)
      ! (280-000 cm^-1 ):  nongray absorption by H2O.
      call qk(ncol,nlm,c18h2o(:,:,ig),tt,fkg,pkd,ip1,ip2)
      call qoph2o(fkg,dp,rmix,tg)       
      hk = hk18(ig)

   case default
      stop
   end select

   return
end subroutine gases


subroutine qopch4 (fkg,   dp,   tg)
!-----------------------------------------------------------------------
! Absorption by CH4 in the 1400-1250cm^-1 and 1250-1100cm^-1
! spectral intervals.
      
   use bugsrad_physconst, only:  molar_volume,gravity,MW_dry_air
   implicit none

! ARGUMENT LIST VARIABLES:
!  INPUT ARGUMENTS:
!  ----------------
   real (kind=dbl_kind), intent(in), dimension(:,:)::  &
     dp      &!Layer Thickness                       (hPa).
    ,fkg      !Gaseous absorption coefficients (cm-atm)^-1.
     
!  OUTPUT ARGUMENTS:
!  -----------------
   real (kind=dbl_kind), intent(out), dimension(:,:)::  &
     tg       !Transmission (-).
     
   real (kind=dbl_kind) :: &
     ch4_conc                   !ppv

   ch4_conc = 1.6e-6_dbl_kind
   tg = fkg*dp*molar_volume*10.0/gravity*ch4_conc/MW_dry_air

   return
end subroutine qopch4


subroutine qopn2o (fkg , dp,  tg)
!-----------------------------------------------------------------------
! Absorption by N2O in the 1400-1250cm^-1 and 1250-1100cm^-1
! spectral intervals.

   use bugsrad_physconst, only:  molar_volume, gravity, MW_dry_air
   implicit none

   real (kind=dbl_kind), intent(in), dimension(:,:)::  &
     dp      & !Layer Thickness                       (hPa).
    ,fkg       !Gaseous absorption coefficients (cm-atm)^-1.
     
   !  OUTPUT ARGUMENTS:
   real (kind=dbl_kind), intent(out), dimension(:,:)::  &
     tg   !Transmission                             (-).
     
   real (kind = dbl_kind)::  &
     n2o_conc              !ppv

   n2o_conc = 0.28e-6_dbl_kind

   tg = fkg*dp*molar_volume*10.0/gravity*n2o_conc/MW_dry_air

   return
end subroutine qopn2o


subroutine qoph2o (fkg , dp ,  rmix ,  tg)
!-----------------------------------------------------------------------
! Absorption by H2O in each SW and LW spectral interval.

   use bugsrad_physconst, only:  molar_volume, gravity, MW_h2o
   implicit none

   real (kind=dbl_kind), intent(in), dimension(:,:)::  &
     dp          & !Layer Thickness                       (hPa).
    ,rmix        & !Water Vapor Mixing Ratio            (kg/kg).
    ,fkg          !Gaseous absorption coefficients (cm-atm)^-1.
     
   !  OUTPUT ARGUMENTS:
   real (kind=dbl_kind), intent(out), dimension(:,:)::  &
     tg   !Transmisssion                           (-).
     
   tg = fkg*rmix*dp*molar_volume/MW_h2o*10.0/gravity

   return
end subroutine qoph2o


subroutine qopo3i (fkg , dp , o3mix ,  tg )
!-----------------------------------------------------------------------
! Absorption by O3 in the infrared band (1100-980 cm^-1).
      
   use bugsrad_physconst, only:  molar_volume, gravity, MW_o3
   implicit none
   
   real (kind=dbl_kind), intent(in), dimension(:,:)::  &
      dp        &  !Layer Thickness                       (hPa).
     ,o3mix     &  !Ozone Mixing Ratio                  (kg/kg).
     ,fkg          !Gaseous absorption coefficients (cm-atm)^-1.
  
   !  OUTPUT ARGUMENTS:
   real (kind=dbl_kind), intent(out), dimension(:,:)::  &
      tg    !Transmission                            (-).
  
   tg = fkg*o3mix*dp*molar_volume/MW_o3*10.0/gravity

   return
end subroutine qopo3i


subroutine qopo3s (fk , dp , o3mix ,  tg )
!-----------------------------------------------------------------------      
! Computes non-gray absorption of O3.

   use bugsrad_physconst, only:  molar_volume, gravity, MW_o3
   implicit none

   real (kind=dbl_kind), intent(in):: &
      fk            !
  
   real (kind=dbl_kind), intent(in), dimension(:,:)::  &
     dp           & !Layer Thickness      (hPa).
    ,o3mix          !Ozone Mixing Ratio (kg/kg).

   !  OUTPUT ARGUMENTS:
   real (kind=dbl_kind), intent(out), dimension(:,:)::  &
     tg             !Transmission           (-).
  
   tg = o3mix*dp*fk*molar_volume/MW_o3*10.0/gravity

   return
end subroutine qopo3s


subroutine qophc(fkg,dp,tg)
!-----------------------------------------------------------------------
! Overlapping of H2O and CO2 bands.
      
   implicit none

   ! INPUT ARGUMENTS
   real (kind=dbl_kind), intent(in), dimension(:,:)::  &
     dp       & !Layer Thickness                       (hPa).
    ,fkg        !Gaseous absorption coefficients (cm-atm)^-1.
     
   !  OUTPUT ARGUMENTS:
   real (kind=dbl_kind), intent(out), dimension(:,:)::  &
     tg   !Transmission                            (-).
     
   tg = fkg*dp
   return
end subroutine qophc

subroutine qk(ncol,nlm,coefk,tt,fkg,pkd,ip1,ip2)
!-----------------------------------------------------------------------
! Interpolates the correlated k coeffients to the temperatures tt
! using the form
!   ln k = a + b * ( t - 245 ) + c * ( t - 245 ) ** 2
! and linearly interpolates in pressure if needed.

   implicit none
      
   ! INPUT ARGUMENTS:
   integer (kind=int_kind), intent(in)::  &
    ncol         & !Length of sub-domain.
   ,nlm            !Number of layers.

   integer (kind=int_kind), intent(in), dimension(:,:)::  &
    ip1          & !Used in conjunction with pressure weigthing.
   ,ip2            !Used in conjunction with pressure weigthing.

   real (kind=dbl_kind), intent(in), dimension(:,:)::  &
    coefk          !Pre-computed coefficients.
     
   real (kind=dbl_kind), intent(in), dimension(:,:)::  &
    tt           & !Layer temperature (K).
   ,pkd            !Scaled pressure (hPa).
     
   ! OUTPUT ARGUMENTS:
   real (kind=dbl_kind), intent(out), dimension(:,:)::  &
    fkg    !Gaseous absorption coefficients (cm-atm)^-1.  

   ! LOCAL VARIABLES:
   integer (kind=int_kind)::  &
    i, l, i1, i2

   real (kind=dbl_kind)::  &
    x1,x2,y1
     

   do l = 1, nlm
      do i = 1, ncol
         i1=ip1(i,l)
         i2=ip2(i,l)
         y1 = tt(i,l) - 245.0
         if (y1.lt.-65.) y1 = -65.
         if (y1.gt.75.) y1 = 75.
         x1 = exp(coefk(1,i1)+coefk(2,i1)*y1+coefk(3,i1)*y1**2)
         if(i1.eq.i2) then
            fkg(i,l) = x1 * pkd(i,l)
         else
            x2 = exp(coefk(1,i2)+coefk(2,i2)*y1+coefk(3,i2)*y1**2)
            fkg(i,l) = x1 + ( x2 - x1 ) * pkd(i,l)
         endif
      enddo
   enddo

   return
end subroutine qk

subroutine qkio3 (ncol , nlm , coefk ,  tt , fkg , pkd ,    ip1 , ip2)
! Computes the gaseous absorption coefficients in units of (cm-atm)^-1
! for a given cumulative probability in nlm layers. coefk are
! the coefficients used to calculate the absorption coefficient at the
! temperature t for the 19 pressures by
!   ln k = a + b * ( t - 250 ) + c * ( t - 250 ) ** 2
! and the absorption coefficient at conditions other than those nineteen
! pressures is interpolated linearly with pressure (Fu, 1991).
      
   implicit none

   ! INPUT ARGUMENTS:
   integer (kind=int_kind), intent(in)::  &
    ncol         & !Length of sub-domain.
   ,nlm            !Number of layers.

   integer (kind=int_kind), intent(in), dimension(:,:)::  &
    ip1          & !Used in conjunction with pressure weigthing.
   ,ip2            !Used in conjunction with pressure weigthing.

   real (kind=dbl_kind), intent(in), dimension(:,:)::  &
    coefk          !Pre-computed coefficients.
  
   real (kind=dbl_kind), intent(in), dimension(:,:)::  &
    tt           & !Layer temperature  (K).
   ,pkd            !Scaled pressure  (hPa).
 
 
   ! OUTPUT ARGUMENTS:
   real (kind=dbl_kind), intent(out), dimension(:,:)::  &
    fkg    !Gaseous absorption coefficients (cm-atm)^-1.  


   ! LOCAL VARIABLES:
   integer (kind=int_kind)::  &
    i, l, i1, i2     

   real (kind=dbl_kind)::  &
    x1,x2,y1
      
   do l = 1, nlm
      do i = 1, ncol
         i1 = ip1(i,l)
         i2 =ip2(i,l)
         y1 = tt(i,l) - 250.0
         if (y1.lt.-70.) y1 = -70.
         if (y1.gt.75.) y1 = 70.
         x1 = exp(coefk(1,i1)+coefk(2,i1)*y1+coefk(3,i1)*y1**2)
         if(i1.eq.i2) then
            fkg(i,l) = x1 * pkd(i,l)
         else
            x2 = exp(coefk(1,i2)+coefk(2,i2)*y1+coefk(3,i2)*y1**2)
            fkg(i,l) = x1 + ( x2 - x1 ) * pkd(i,l)
         endif
      enddo
   enddo

   return
end subroutine qkio3

subroutine pscale (ncol , nlm , ppl, stanp, pkd , ip1 , ip2 )
!-----------------------------------------------------------------------
! pscale is used to scale pressure for use in the k-distribu
! tion routine.
 
   implicit none

   ! INPUT ARGUMENTS:
   integer (kind=int_kind), intent(in)::  &
    ncol         & !Length of sub-domain.
   ,nlm            !Number of layers.
     
   real (kind=dbl_kind), dimension(:,:)::  &
    ppl            !Pressure          (hPa).

   real (kind=dbl_kind), dimension(:):: &
    stanp

   ! OUTPUT ARGUMENTS:
   integer (kind=int_kind), dimension(:,:)::  &
    ip1            & !pointer.
   ,ip2              !pointer.
     
   real (kind=dbl_kind), dimension(:,:)::  &
    pkd              !Weighted pressure (hPa).

   ! LOCAL VARIABLES:
   integer (kind=int_kind)::  &
    n_stanp

   integer (kind=int_kind)::  &
    i1, icol, l

   n_stanp = size(stanp)

   do icol = 1, ncol
      i1 = 1
      do l = 1, nlm
         if( ppl(icol,l) .lt. stanp(1) ) then
            pkd(icol,l) = ppl(icol,l) / stanp(1)
            ip1(icol,l) = 1
            ip2(icol,l) = 1
         elseif ( ppl(icol,l) .ge. stanp(n_stanp) ) then
            pkd(icol,l) = (ppl(icol,l) - stanp(n_stanp-1))  &
                        / (stanp(n_stanp) - stanp(n_stanp-1))
            ip1(icol,l) = n_stanp - 1
            ip2(icol,l) = n_stanp
         else
30          continue
            if ( ppl(icol,l) .ge. stanp(i1) ) goto 20
            pkd(icol,l) = (ppl(icol,l)-stanp(i1-1))  &
                        / (stanp(i1)-stanp(i1-1))
            ip1(icol,l) = i1-1
            ip2(icol,l) = i1
            goto 5
20          i1 = i1 + 1
            goto 30
         endif
5     enddo
   enddo
   return
end subroutine pscale


!-----------------------------------------------------------------------
end module gases_ckd
