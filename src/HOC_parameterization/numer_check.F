!------------------------------------------------------------------------
! $Id: numer_check.F,v 1.8 2008-02-01 21:31:55 faschinj Exp $
        module numerical_check

        implicit none

!       Made isnan2d public so it may be used
!       for finding code that cause NaNs
!       Joshua Fasching November 2007

!       *_check subroutines were added to ensure that the
!       subroutines they are checking perform correctly
!       Joshua Fasching February 2008

        public :: invalid_model_arrays, isnan2d, 
     .            rad_clipping, parameterization_check,
     .            sfc_var_check, pdf_closure_new_check
        
        private :: check_negative, nan_report, isNaN


!       Abstraction of nan_report
        interface nan_report
           module procedure nan_report_sclr, nan_report_2d
        end interface
        
        contains
!--------------------------------------------------------------------------------- 
         subroutine length_check( Lscale, lup, ldown, isValid )
!
!        Description: This subroutine determines if any of the output
!        variables for the length_new subroutine carry values that
!        are NaNs.
!
!        Joshua Fasching February 2008
!
!---------------------------------------------------------------------------------
         use grid_class
         implicit none

         ! Input Variables
         real, dimension(gr%nnzp), intent(in) :: 
     .   Lscale, ! Mixing length                 [m]
     .   lup,    ! Upward mixing length          [m]
     .   ldown   ! Downward mixing length        [m]
         
         ! Output Variable
         logical, intent(out) ::
     .   isValid

         ! Local Variables
         character(*), parameter :: proc_name = "length"
!-----------------------------------------------------------------------------
         call nan_report( Lscale, "Lscale", proc_name, isValid )
         call nan_report( lup, "lup", proc_name, isValid )
         call nan_report( ldown, "ldown", proc_name, isValid )

         end subroutine length_check
!---------------------------------------------------------------------------
         subroutine pdf_closure_new_check( wp4, wprtp2, wp2rtp, wpthlp2,
     .                           wp2thlp, cf, rcm, wpthvp, wp2thvp,
     .                           rtpthvp, thlpthvp, wprcp, wp2rcp,
     .                           rtprcp, thlprcp, rcp2, wprtpthlp,
     .                           crt1, crt2, cthl1, cthl2, pdf_parms,
     .                           isValid, sclrpthvp, sclrprcp, wpsclrp2,
     .                           wpsclrprtp, wpsclrpthlp, wp2sclrp )
!        Description: This subroutine determines if any of the output
!        variables for the pdf_closure_new subroutine carry values that
!        are NaNs.
!
!        Joshua Fasching February 2008
!---------------------------------------------------------------------------
         use constants, only:sclr_dim
         implicit none

         ! Input Variables
         real, intent(in) ::
     .   wp4,            ! w'^4                  [m^4/s^4]
     .   wprtp2,         ! w' r_t'               [(m kg)/(s kg)]
     .   wp2rtp,         ! w'^2 r_t'             [(m^2 kg)/(s^2 kg)]
     .   wpthlp2,        ! w' th_l'^2            [(m K^2)/s]
     .   wp2thlp,        ! w'^2 th_l'            [(m^2 K)/s^2]
     .   cf,             ! Cloud fraction        [%]
     .   rcm,            ! Mean liquid water     [kg/kg]
     .   wpthvp,         ! Buoyancy flux         [(K m)/s] 
     .   wp2thvp,        ! w'^2 th_v'            [(m^2 K)/s^2]
     .   rtpthvp,        ! r_t' th_v'            [(kg K)/kg]
     .   thlpthvp,       ! th_l' th_v'           [K^2]
     .   wprcp,          ! w' r_c'               [(m kg)/(s kg)]
     .   wp2rcp,         ! w'^2 r_c'             [(m^2 kg)/(s^2 kg)]
     .   rtprcp,         ! r_t' r_c'             [(kg^2)/(kg^2)]
     .   thlprcp,        ! th_l' r_c'            [(K kg)/kg]
     .   rcp2,           ! r_c'^2                [(kg^2)/(kg^2)]
     .   wprtpthlp,      ! w' r_t' th_l'         [(m kg K)/(s kg)]
     .   crt1, crt2, 
     .   cthl1, cthl2

        real, intent(in), dimension(26):: 
     .  pdf_parms        ! pdf paramters         [units vary]
   
        ! Input (Optional passive scalar variables)  
        real, optional, dimension(sclr_dim),intent(in) :: 
     .  sclrpthvp, 
     .  sclrprcp, 
     .  wpsclrp2,
     .  wpsclrprtp,
     .  wpsclrpthlp,
     .  wp2sclrp
        
        ! Output Variable        
        logical, intent(out) :: 
     .  isValid          ! Are the outputs usable numbers?
        
        ! Local Variable
        character(*), parameter :: proc_name = "pdf_closure_new"
!-------------------------------------------------------------------------------
        call nan_report( wp4,"wp4", proc_name, isValid )
        call nan_report( wprtp2,"wprtp2", proc_name, isValid )
        call nan_report( wp2rtp,"wp2rtp", proc_name, isValid )
        call nan_report( wpthlp2,"wpthlp2", proc_name, isValid )
        call nan_report( wp2thlp,"wp2thlp", proc_name, isValid )
        call nan_report( cf,"cf", proc_name, isValid )
        call nan_report( rcm,"rcm", proc_name, isValid ) 
        call nan_report( wpthvp, "wpthvp", proc_name, isValid )
        call nan_report( wp2thvp, "wp2thvp", proc_name, isValid )
        call nan_report( rtpthvp, "rtpthvp", proc_name, isValid )
        call nan_report( thlpthvp, "thlpthvp", proc_name, isValid )
        call nan_report( wprcp, "wprcp", proc_name, isValid )
        call nan_report( wp2rcp, "wp2rcp", proc_name, isValid )
        call nan_report( rtprcp, "rtprcp", proc_name, isValid )
        call nan_report( thlprcp, "thlprcp", proc_name, isValid )
        call nan_report( rcp2, "rcp2", proc_name, isValid)
        call nan_report( wprtpthlp, "wprtpthlp", proc_name, isValid )
        call nan_report( crt1, "crt1", proc_name, isValid )
        call nan_report( crt2, "crt2", proc_name, isValid )
        call nan_report( cthl1, "cthl1", proc_name, isValid )
        call nan_report( cthl2, "cthl2", proc_name, isValid ) 
        call nan_report( pdf_parms, "pdf_parms", proc_name,isValid )
        
        if( present(sclrpthvp) ) then
            call nan_report(sclrpthvp,"sclrpthvp",
     .                      proc_name, isValid)
        endif
        if( present(sclrprcp) ) then
             call nan_report(sclrprcp,"sclrprcp",
     .                    proc_name,isValid)
        endif
        if( present(wpsclrprtp) ) then
             call nan_report(wpsclrprtp,"wpsclrprtp", 
     .                    proc_name,isValid)
        endif   
        if( present(wpsclrp2) ) then
             call nan_report(wpsclrp2,"wpsclrp2", 
     .                    proc_name,isValid)
        endif
        if( present(wpsclrpthlp) ) then
             call nan_report(wpsclrpthlp,"wpsclrtlp", 
     .                    proc_name,isValid)
        endif
        if( present(wp2sclrp) ) then
             call nan_report(wp2sclrp,"wps2clrp", 
     .                    proc_name,isValid)
        endif
        end subroutine pdf_closure_new_check
!-------------------------------------------------------------------------------------
        subroutine parameterization_check( thlm_forcing, rtm_forcing,
     .                           rtm_mc, thlm_mc, wmm, wmt, p, rhom, 
     .                           rhot, exner, wpthlp_sfc, wprtp_sfc, 
     .                           upwp_sfc, vpwp_sfc, um, upwp, vm, vpwp,
     .                           up2, vp2, rtm, wprtp, thlm, 
     .                           wpthlp, wp2, wp3, Scm, rtp2, thlp2,
     .                           rtpthlp, taum, rcm, Ncm, Ncnm, Nim,
     .                           hydromet, cf 
#ifdef SCALARS
     /                           ,wpsclrp_sfc, wpedsclrp_sfc,
                                 sclrm, sclrm_forcing, edsclrm
#endif /*SCALARS*/ 

     .  )
!
!       Description: This subroutine determines what input variables may have NaN values.
!       In addition it checks to see if rhom, rhot, exner, up2, vp2, rtm, thlm, wp2, Scm,
!       rtp2, thlp2, taum, rcm, Ncm, Ncnm, Nim, hydromet, or cf have negative values.
!
!
!----------------------------------------------------------------------------------------
        use grid_class
        use constants   
        implicit none

        real, intent(in), dimension(gr%nnzp) :: 
     .  thlm_forcing,  ! theta_l forcing.        [K/s]
     .  rtm_forcing,   ! r_t forcing.            [(kg/kg)/s] 
     .  rtm_mc,        ! r_t microphysical       [(kg/kg)/s] 
     .  thlm_mc,       ! thlm microphysical      [K/s] 
     .  wmm,           ! wm on moment. grid.     [m/s]
     .  wmt,           ! wm on thermo. grid.     [m/s]
     .  p,             ! Pressure.               [Pa] 
     .  rhom,          ! Density on moment. grid [kg/m^3]
     .  rhot,          ! Density on thermo. grid [kg/m^3] 
     .  exner          ! Exner function.         [-]

        real, intent(in) :: 
     .  wpthlp_sfc,  ! w' theta_l' at surface.   [(m K)/s]
     .  wprtp_sfc,   ! w' r_t' at surface.       [(kg m)/( kg s)]
     .  upwp_sfc,    ! u'w' at surface.          [m^2/s^2]
     .  vpwp_sfc     ! v'w' at surface.          [m^2/s^2]

        ! These are prognostic or are planned to be in the future
        real, intent(in), dimension(gr%nnzp) :: 
     .  um,      ! u wind.                       [m/s]
     .  upwp,    ! u'w'.                         [m^2/s^2]
     .  vm,      ! v wind.                       [m/s]
     .  vpwp,    ! u'w'.                         [m^2/s^2]
     .  up2,     ! u'^2                          [m^2/s^2]
     .  vp2,     ! v'^2                          [m^2/s^2]
     .  rtm,     ! r_t Total water mixing ratio. [kg/kg]
     .  wprtp,   ! w' r_t'.                      [(m kg)/(s kg)]
     .  thlm,    ! th_l Liquid potential temp.   [K]
     .  wpthlp,  ! w' th_l'.                     [(m K)/s]
     .  wp2,     ! w'^2.                         [m^2/s^2]
     .  wp3,     ! w'^3.                         [m^3/s^3]
     .  Scm,     ! Sc on moment. grid.           [-]
     .  rtp2,    ! r_t'^2.                       [(kg/kg)^2]
     .  thlp2,   ! th_l'^2.                      [K^2]
     .  rtpthlp, ! r_t' th_l'.                   [(kg K)/kg]
     .  taum,    ! Tau on moment. grid.          [s]
     .  rcm,     ! Liquid water mixing ratio.    [kg/kg]
     .  Ncm,     ! Cloud droplet number conc.    [num/kg]
     .  Ncnm,    ! Cloud nuclei number conc.     [num/m^3]
     .  Nim      ! Ice crystal number conc.      [num/m^3]

        real, intent(in), dimension(gr%nnzp,hydromet_dim) ::
     .  hydromet  ! Hydrometer fields                     [units vary]
        ! (:,1) rrm      Rain water mixing ratio          [kg/kg]
        ! (:,2) Nrm      Rain droplet number conc.        [num/kg]
        ! (:,3) rsnow    Snow water mixing ratio          [kg/kg]
        ! (:,4) rice     Ice water mixing ratio           [kg/kg]
        ! (:,5) rgraupel Graupel water mixing ratio       [kg/kg]

        real, intent(in), dimension(gr%nnzp) :: 
     .  cf ! Cloud fraction.     [%]

#ifdef SCALARS
        ! Optional Input Variables
         real, intent(in),  dimension(sclr_dim) ::
     .   wpsclrp_sfc,  ! Scalar flux at surface [units m/s]
     .   wpedsclrp_sfc ! Eddy-Scalar flux at surface      [units m/s]

        ! Optional Input/Output Variables
        real, intent(in),dimension(gr%nnzp,sclr_dim) ::
     .  sclrm,          ! Passive scalar mean.
     .  sclrm_forcing,  ! Passive scalar forcing.
     .  edsclrm         !Eddy passive scalar mean.

#endif /*SCALARS*/
        ! Local Variables
                                                    
        ! Name of the procedure using parameterization_check
        character(*), parameter :: 
     .  proc_name = "parameterization timestep"
        
!-------- Input Nan Check ----------------------------------------------
        call nan_report( thlm_forcing, "thlm_forcing", proc_name )
        call nan_report( rtm_forcing,"rtm_forcing", proc_name )
        call nan_report( rtm_mc, "rtm_mc", proc_name )
        call nan_report( thlm_mc, "thlm_mc", proc_name )
        call nan_report( wmm,"wmm", proc_name )
        call nan_report( wmt,"wmt", proc_name )
        call nan_report( p,"p", proc_name )
        call nan_report( rhom,"rhom", proc_name )
        call nan_report( rhot,"rhot", proc_name )
        call nan_report( exner,"exner", proc_name )
        call nan_report( um,"um", proc_name )
        call nan_report( upwp,"upwp", proc_name )
        call nan_report( vm,"vm", proc_name )
        call nan_report( vpwp,"vpwp", proc_name )
        call nan_report( up2,"up2", proc_name )
        call nan_report( vp2,"vp2", proc_name )
        call nan_report( rtm,"rtm", proc_name )
        call nan_report( wprtp,"wprtp", proc_name )
        call nan_report( thlm,"thlm", proc_name )
        call nan_report( wpthlp,"wpthlp", proc_name )
        call nan_report( wp2,"wp2", proc_name )
        call nan_report( wp3,"wp3", proc_name )
        call nan_report( Scm,"Scm", proc_name )
        call nan_report( rtp2,"rtp2", proc_name )
        call nan_report( thlp2,"thlp2", proc_name )
        call nan_report( rtpthlp, "rtpthlp", proc_name )
        call nan_report( taum,"taum", proc_name )
        call nan_report( rcm,"rcm", proc_name )
        call nan_report( Ncm,"Ncm", proc_name )
        call nan_report( Ncnm,"Ncnm", proc_name )
        call nan_report( Nim,"Nim", proc_name )
        call nan_report( hydromet(:,1),"rrm", proc_name )
        call nan_report( hydromet(:,2),"Nrm", proc_name )
        call nan_report( hydromet(:,3),"rsnowm", proc_name )
        call nan_report( hydromet(:,4),"ricem", proc_name )
        call nan_report( hydromet(:,5),"rgraupel", proc_name )
        call nan_report( cf,"cf", proc_name )
        
#ifdef SCALARS
        call nan_report( wpsclrp_sfc,"wpsclrp_sfc", proc_name )
        call nan_report( wpedsclrp_sfc,"wpedsclrp_sfc", proc_name )
        call nan_report( sclrm,"sclrm", proc_name )
        call nan_report( sclrm_forcing,"sclrm_forcing", proc_name ) 
        call nan_report( edsclrm,"edsclrm", proc_name )
#endif /*SCALARS*/ 
!---------------------------------------------------------------------
        call check_negative( p,"p", proc_name )
        call check_negative( rhot,"rhot", proc_name )
        call check_negative( rhom,"rhom", proc_name )
        call check_negative( exner,"exner", proc_name )
        call check_negative( up2,"up2", proc_name )
        call check_negative( vp2,"vp2", proc_name )
        call check_negative( wp2,"wp2", proc_name )
        call check_negative( rtm,"rtm", proc_name )
        call check_negative( thlm,"thlm", proc_name )
        call check_negative( Scm,"Scm", proc_name )
        call check_negative( rtp2,"rtp2", proc_name )
        call check_negative( thlp2,"thlp2", proc_name )
        call check_negative( taum,"taum", proc_name )
        call check_negative( rcm,"rcm", proc_name )
        call check_negative( Ncm,"Ncm", proc_name )
        call check_negative( Ncnm,"Ncnm", proc_name )
        call check_negative( Nim,"Nim", proc_name )
        call check_negative( hydromet(:,1),"rrm", proc_name )
        call check_negative( hydromet(:,2),"Nrm", proc_name )
        call check_negative( hydromet(:,3),"rsnowm", proc_name )
        call check_negative( hydromet(:,4),"ricem", proc_name )
        call check_negative( hydromet(:,5),"rgraupelm", proc_name )
        call check_negative( cf,"cf", proc_name )
       
!-------------------------------------------------------------
        end subroutine parameterization_check
!-----------------------------------------------------------------------
        subroutine sfc_var_check( wp2_sfc, up2_sfc, vp2_sfc, thlp2_sfc,
     .           rtp2_sfc, rtpthlp_sfc, isValid,
     .           sclrp2_sfc, sclrprtp_sfc, sclrpthlp_sfc )
!
!       Description:This subroutine determines if any of the output
!       variables for the sfc_var subroutine carry values that
!       are NaNs.
!
!       Joshua Fasching February 2008
!
!
!-----------------------------------------------------------------------
        use constants
        implicit none
        
        ! Input Variables
        real,intent(in) :: 
     .  wp2_sfc,    ! Vertical velocity variance        [m^2/s^2]
     .  up2_sfc,    ! u'^2                              [m^2/s^2]
     .  vp2_sfc,    ! u'^2                              [m^2/s^2]
     .  thlp2_sfc,  ! thetal variance                   [K^2]
     .  rtp2_sfc,   ! rt variance                       [kg/kg]
     .  rtpthlp_sfc ! thetal rt covariance              [kg K]
  

        real, optional, dimension(sclr_dim), intent(in) ::
     .  sclrp2_sfc,   ! Passive scalar variance                 [units^2]
     .  sclrprtp_sfc, ! Passive scalar r_t covariance           [units kg/kg]
     .  sclrpthlp_sfc ! Passive scalar theta_l covariance       [units K]
     
        ! Output Variable
        logical, intent(out) :: isValid    ! Are these outputs valid?

        ! Local Variables

        ! Name of the subroutine calling the check
        character(*),parameter :: proc_name = "sfc_var"

        integer :: k

        call nan_report( wp2_sfc, "wp2_sfc", proc_name, isValid)
        call nan_report( up2_sfc, "up2_sfc", proc_name, isValid)
        call nan_report( vp2_sfc, "vp2_sfc", proc_name, isValid)
        call nan_report( thlp2_sfc, "thlp2_sfc", proc_name, isValid)
        call nan_report( rtp2_sfc, "rtp2_sfc", proc_name, isValid)
        call nan_report( rtpthlp_sfc, "rtpthlp_sfc", 
     .                  proc_name, isValid)
        
        if( present(sclrp2_sfc) ) then
                call nan_report(sclrp2_sfc,"sclrp2_sfc",
     .                      proc_name, isValid)
        endif
        if( present(sclrprtp_sfc) ) then
             call nan_report(sclrprtp_sfc,"sclrprtp_sfc",
     .                    proc_name,isValid)
             endif
        if( present(sclrp2_sfc) ) then
             call nan_report(sclrpthlp_sfc,"sclrpthlp_sfc", 
     .                    proc_name,isValid)
        endif
 
        end subroutine sfc_var_check 
!-----------------------------------------------------------------------
        subroutine mixing_clipping(rtm,thlm)
!
!
!       Description: This subroutine clips negative values caused by the mixing
!       subroutine. The clipped values are set to zero.
!
!------------------------------------------------------------------------
        use constants, only: fstderr
        use grid_class
        implicit none
        
        ! Input/Output variables
        real, dimension(gr%nnzp), intent(inout) ::
     .  rtm,    ! Total Water Mixing Ration            [Kg/Kg]
     .  thlm    ! Liquid Water Potential Temperature   [K/s]

        ! Local variables
        integer :: k 
!-------------------------------------------------------------------------
        do k = 1, gr%nnzp, 1
                if ( thlm(k) < 0.0 ) then
                        write(fstderr,*) "Negative value in thlm at ",k,
     .                  " caused in mixing"              
                        thlm(k) = 0.0
                endif
                if( rtm(k) < 0.0 ) then
                        write(fstderr,*) "Negative value in rtm at ", k,
     .                  " caused in mixing" 
                        rtm(k) = 0
                endif

        end do        
        return
        end subroutine mixing_clipping
!-----------------------------------------------------------------------
        subroutine rad_clipping(nz,thlm,rcm,rtm,
     .     ricem, 
     .     cf,p,exner,rhom)
!       Description:
!       Clips radiation input variables if they are < 0 and reports
!       its action to the console
!------------------------------------------------------------------------
        use constants, only: fstderr
        implicit none
        
        ! Input variable
        integer, intent(in) :: nz
        !Input/Output variables
        real, dimension(nz), intent(inout) ::
     .  thlm,          ! Liquid Water Potential Temperature   [K/s]
     .  rcm,           ! Liquid Water Mixing Ratio            [Kg/Kg]
     .  rtm,           ! Total Water Mixing Ratio             [Kg/Kg]
     .  ricem,         ! Ice Water Mixing Ratio               [Kg/Kg]
     .  cf,            ! Cloud Fraction                       [%]
     .  p,             ! Pressure                             [Pa]
     .  exner,         ! Exner Function                       [-]
     .  rhom           ! Density                              [-]

       ! Local variables
        real,dimension(nz) :: rvm
        integer :: k 
!-------------------------------------------------------------------------
        rvm = rtm - rcm
        do k = 1, nz, 1
                if ( thlm(k) < 0.0 ) then
                        write(fstderr,*) "Negative value in thlm at ",k,
     .                  " before BUGSrad." 
                        thlm(k) = 0.0
                endif
                if( rcm(k) < 0.0 ) then
                        write(fstderr,*) "Negative value in rcm at ", k,
     .                  " before BUGSrad." 
                        rcm(k) = 0.0
                endif
                if( rtm(k) < 0.0 ) then
                        write(fstderr,*) "Negative value in rtm at ", k,
     .                  " before BUGSrad." 
                        rtm(k) = 0
                endif
                if( rvm(k) < 0.0 ) then
                        write(fstderr,*) "Negative value in rvm at ", k,
     .                  " before BUGSrad." 
                        rcm(k) = rtm(k)
                endif
                if ( ricem(k) < 0.0 ) then
                        write(fstderr,*) "Negative value",
     .                                           " in ricem at ", k,
     .                  " before BUGSrad." 
                        ricem(k) = 0.0
                endif
                if ( cf(k) < 0.0 ) then
                        write(fstderr,*) "Negative value in cf at ", k,
     .                  " before BUGSrad." 
                        cf(k) = 0.0
                endif
                if (p(k) < 0.0 ) then
                        write(fstderr,*) "Negative value in p at ", k,
     .                  " before BUGSrad." 
                        p(k) = 0.0
                endif
                if (exner(k) < 0.0) then
                        write(fstderr,*) "Negative value", 
     .                                           " in exner at ", k,
     .                  " before BUGSrad." 
                        exner(k) = 0.0
                endif
                if (rhom(k) < 0.0) then
                        write(fstderr,*) "Negative value in rhom at ",k,
     .                  " before BUGSrad." 
                        rhom(k) = 0.0
                endif
        end do        
        return
        end subroutine rad_clipping

!-----------------------------------------------------------------------
        logical function invalid_model_arrays( ) 

!       Description:
!       Checks for invalid floating point values in select model arrays.

!       References:
!       None
!------------------------------------------------------------------------

        use diagnostic_variables
        use prognostic_variables
        use model_flags, only: lcoamps_micro, kk_rain
        use constants, only: sclr_dim, hydromet_dim, fstderr

        implicit none

        ! Local Variables
        integer :: i

        invalid_model_arrays = .false.
           
        ! Check whether any variable array contains a NaN for
        ! um, vm, thlm, rtm, rtp2, thlp2, wprtp, wpthlp, rtpthlp, 
        ! wp2, & wp3.        
        if ( isnan2d( um ) ) then
          write(fstderr,*) "NaN in um model array" 
!         write(fstderr,*) "um= ", um
          invalid_model_arrays = .true.
!         return
        end if

        if ( isnan2d( vm ) ) then
          write(fstderr,*) "NaN in vm model array" 
!         write(fstderr,*) "vm= ", vm
          invalid_model_arrays = .true.
!         return
        end if

        if ( isnan2d( wp2 ) ) then
          write(fstderr,*) "NaN in wp2 model array" 
!         write(fstderr,*) "wp2= ", wp2
          invalid_model_arrays = .true.
!         return
        end if

        if ( isnan2d( wp3 ) ) then
          write(fstderr,*) "NaN in wp3 model array" 
!         write(fstderr,*) "wp3= ", wp3
          invalid_model_arrays = .true.
!         return
        end if

        if ( isnan2d( rtm ) ) then
          write(fstderr,*) "NaN in rtm model array" 
!         write(fstderr,*) "rtm= ", rtm
          invalid_model_arrays = .true.
!         return
        end if

        if ( isnan2d( thlm ) ) then
          write(fstderr,*) "NaN in thlm model array" 
!         write(fstderr,*) "thlm= ", thlm
          invalid_model_arrays = .true.
!         return
        end if

        if ( isnan2d( rtp2 ) ) then
          write(fstderr,*) "NaN in rtp2 model array" 
!         write(fstderr,*) "rtp2= ", rtp2
          invalid_model_arrays = .true.
!         return
        end if

        if ( isnan2d( thlp2 ) ) then
          write(fstderr,*) "NaN in thlp2 model array" 
!         write(fstderr,*) "thlp2= ", thlp2
          invalid_model_arrays = .true.
!         return
        end if

        if ( isnan2d( wprtp ) ) then
          write(fstderr,*) "NaN in wprtp model array" 
!         write(fstderr,*) "wprtp= ", wprtp
          invalid_model_arrays = .true.
!         return
        end if

        if ( isnan2d( wpthlp ) ) then
          write(fstderr,*) "NaN in wpthlp model array" 
!         write(fstderr,*) "wpthlp= ", wpthlp
          invalid_model_arrays = .true.
!         return
        end if

        if ( isnan2d( rtpthlp ) ) then
          write(fstderr,*) "NaN in rtpthlp model array" 
!         write(fstderr,*) "rtpthlp= ", rtpthlp
          invalid_model_arrays = .true.
!         return
        end if

        if ( kk_rain .or. lcoamps_micro ) then
          do i = 1, hydromet_dim, 1
            if ( isnan2d( hydromet(:,i) ) ) then
              write(fstderr,*) "NaN in a hydrometeor model array" 
!             write(fstderr,*) "hydromet= ", hydromet
              invalid_model_arrays = .true.
!             return
            end if
          end do
        end if

!       if ( isnan2d( wmt ) ) then
!         write(fstderr,*) "NaN in wmt model array" 
!         write(fstderr,*) "wmt= ", wmt
!         invalid_model_arrays = .true.
!         return
!       end if

        if ( isnan2d( wp2thvp ) ) then
          write(fstderr,*) "NaN in wp2thvp model array" 
!         write(fstderr,*) "wp2thvp = ", wp2thvp
          invalid_model_arrays = .true.
!         return
        end if

        if ( isnan2d( rtpthvp ) ) then
          write(fstderr,*) "NaN in rtpthvp model array" 
!         write(fstderr,*) "rtpthvp = ", rtpthvp
          invalid_model_arrays = .true.
        end if

        if ( isnan2d( thlpthvp ) ) then
          write(fstderr,*) "NaN in thlpthvp model array" 
!         write(fstderr,*) "thlpthvp = ", thlpthvp
          invalid_model_arrays = .true.
        end if

#ifdef SCALARS
        do i = 1, sclr_dim, 1
          if ( isnan2d( sclrm(:,i) ) ) then
            write(fstderr,*) "NaN in sclrm", i, "model array"
!           write(fstderr,'(a6,i2,a1)') "sclrm(", i, ")"
!           write(fstderr,*) sclrm(:,i)
            invalid_model_arrays = .true.
          end if
          if ( isnan2d( edsclrm(:,i) ) ) then
            write(fstderr,*) "NaN in `edsclrm", i, "model array"
!           write(fstderr,'(a8,i2,a1)') "edsclrm(", i, ")"
!           write(fstderr,*) edsclrm(:,i)
            invalid_model_arrays = .true.
          end if
        end do
#endif /*SCALARS*/

        return
        end function invalid_model_arrays
!------------------------------------------------------------------------
        logical function isNaN( x2d )

!       Description:
!       Checks if a given real is a NaN, +inf or -inf.

!       Notes: 
!       Got a tip from Andy Vaught to use transfer( ) for portability.
!       Ideally this would be done with intrinsic module ieee in 
!       Fortran 2003 or with the isnan( ) extension provided by some 
!       (not pgroup) compilers.

!       Certain compiler optimizations may cause variables with invalid
!       results to flush to zero.  Avoid these!
!------------------------------------------------------------------------
        
        implicit none

        ! External
        intrinsic :: transfer
    
        ! Input Variables
        real, intent(in) :: x2d

        ! Local Variables
        integer(kind=4) :: nanbits

        real(kind=4) :: PosInf, NegInf

        data nanbits /Z"7F800000"/

        PosInf = transfer( nanbits, PosInf )
        NegInf = -( transfer( nanbits, NegInf ) )

        ! This works on compilers with standardized floating point, 
        ! because the IEEE 754 spec defines that subnormals and NaNs 
        ! should not equal themselves.
        ! However, all compilers do not seem to follow this. 
        if (x2d /= x2d ) then
          isNaN = .true.

        ! This a second check, assuming the above does not work as 
        ! expected.
        else if ( x2d == PosInf ) then
          isNaN = .true.

        ! This may never be needed, it's here just in case.
        else if ( x2d == NegInf ) then
          isNaN = .true.

        else 
          isNaN = .false. ! Yippee, real numbers.

        end if

        return
        end function isNaN
!------------------------------------------------------------------------

!------------------------------------------------------------------------
        logical function isnan2d( x2d )

!       Description:
!       Checks if a given real is a NaN, +inf or -inf.

!       Notes: 
!       Got a tip from Andy Vaught to use transfer( ) for portability.
!       Ideally this would be done with intrinsic module ieee in 
!       Fortran 2003 or with the isnan( ) extension provided by some 
!       (not pgroup) compilers.

!       Certain compiler optimizations may cause variables with invalid
!       results to flush to zero.  Avoid these!
!------------------------------------------------------------------------
        
        implicit none
        intrinsic :: any,transfer
        ! Input Variables
        real, dimension(:), intent(in) :: x2d
        ! Local Variables
        integer(kind=4) :: nanbits

        real(kind=4) :: PosInf, NegInf

        data nanbits /Z"7F800000"/

        PosInf = transfer( nanbits, PosInf )
        NegInf = -( transfer( nanbits, NegInf ) )

        ! This works on compilers with standardized floating point, 
        ! because the IEEE 754 spec defines that subnormals and NaNs 
        ! should not equal themselves.
        ! However, all compilers do not seem to follow this. 
        if ( any( x2d /= x2d ) ) then
          isnan2d = .true.

        ! This a second check, assuming the above does not work as 
        ! expected.
        else if ( any( x2d == PosInf ) ) then
          isnan2d = .true.

        ! This may never be needed, it's here just in case.
        else if ( any( x2d == NegInf ) ) then
          isnan2d = .true.

        else 
          isnan2d = .false. ! Yippee, real numbers.

        end if

        return
        end function isnan2d
!------------------------------------------------------------------------
        subroutine check_negative( var, varname, operation, isValid )
!
!       Description: Checks for negative values in the var array and reports
!       them.
!
!-----------------------------------------------------------------------
        use constants, only:fstderr
        use grid_class
        
        implicit none

        ! Input Variables
        real, intent(in), dimension(:) :: var
        
        character(len=*), intent(in):: 
     .      varname,                   ! Varible being examined
     .      operation                  ! Procedure calling check_zero

        ! Output Variable
        logical, optional, intent(inout) :: isValid



        if(any(var < 0.0)) then
           write(fstderr,*) varname," is less than zero in ",operation
           if(present(isValid)) isValid = .false.
        endif
         
        end subroutine check_negative
!------------------------------------------------------------------------        
        subroutine nan_report_2d( var, varname, operation, isValid )
!
!       Description: Checks for a NaN in the var array and reports it.
!
!
!------------------------------------------------------------------------
        use constants, only:fstderr
        use grid_class
        
        implicit none
        ! Input variables
        real, intent(in), dimension(:) :: var ! Variable being examined

        character(len=*), intent(in):: 
     .      varname,                   ! Name of variable
     .      operation                  ! Procedure calling nan_report
        
        logical,optional,intent(out) :: isValid
        
        if (isnan2d(var)) then
           write(fstderr,*) varname, " is NaN in ",operation
           if(present(isValid)) isValid = .false.
        endif
        end subroutine nan_report_2d
!-----------------------------------------------------------------------        
        subroutine nan_report_sclr( var, varname, operation, isValid )
!
!       Description: Checks for a NaN in the scalar var then reports it.
!
!-----------------------------------------------------------------------        
        use constants, only:fstderr
        use grid_class
        
        implicit none
        ! Input Variables
        real, intent(in) :: var        ! Variable being examined

        character(len=*), intent(in):: 
     .      varname,                   ! Name of variable being examined
     .      operation                  ! Procedure calling nan_report
        
        ! Output variable
        logical, optional, intent(out) :: isValid 
!--------------------------------------------------------------------
        if (isNaN(var)) then
           write(fstderr,*) varname, " is NaN in ",operation
           if (present(isValid)) isValid = .false.
        endif
        end subroutine nan_report_sclr
!-------------------------------------------------------------------------       
        end module numerical_check
