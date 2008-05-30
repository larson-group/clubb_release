! $Id: mixing_length.F,v 1.12 2008-05-30 16:24:24 faschinj Exp $
!-----------------------------------------------------------------------------
        module mixing_length
        
        implicit none 
        
        private ! Default Scope
        
        public :: compute_length
        
        contains

!---------------------------------------------------------------------------        
        subroutine compute_length( thvm, thlm, rtm, rcm, em, p, 
     .      exner, Lscale, err_code ) 
!       Description:
!       Larson's 5th moist, nonlocal length scale

!       References:
!       Section 3b ( /Eddy length formulation/ ) of
!       ``A PDF-Based Model for Boundary Layer Clouds. Part I:
!         Method and Model Description'' Golaz, et al. (2002)
!       JAS, Vol. 59, pp. 3540--3551.

!-----------------------------------------------------------------------------

!       mu = (1/M)dM/dz > 0.  mu=0 for no entrainment.
!       Siebesma recommends mu=2e-3, although most schemes use mu=1e-4
!       When mu was fixed, we used the value mu = 6.e-4

        use constants, only: 
            ! Variable(s)
     .      Cp,   ! Dry air specific heat at constant p              [J/kg/K]
     .      Rd,   ! Dry air gas constant                             [J/kg/K]
     .      ep,   ! Rd / Rv                                          [-]
     .      ep1,  ! (1-ep)/ep                                        [-]
     .      ep2,  ! 1/ep                                             [-]
     .      Lv,   ! Latent heat of vaporiztion                       [J/kg/K]
     .      grav, ! Gravitational acceleration                       [m/s^2]
     .      Lscale_max, ! Maximum value for Lscale                   [m]
     .      fstderr

        use parameters, only:
        ! Variable(s)
     .      mu,  ! Fractional entrainment rate per unit altitude    [1/m]
     .      T0,  ! Reference temperature                            [K]
     .      lmin ! Minimum value for Lscale                         [m]

        use grid_class, only:
     .      gr, ! Variable(s)
     .      zm2zt ! Procedure(s)

        use numerical_check, only: 
     .      length_check ! Procedure(s)

        use saturation, only: 
     .      rsat ! Procedure(s)

        use diagnostic_variables, only: 
     .      Lup, ! Variable(s)
     .      Ldown

        use error_code, only: 
     .      clubb_var_equals_NaN, ! Variable(s)
     .      clubb_at_debug_level ! Procedure(s)

        implicit none
        
        ! External
        intrinsic :: max, sqrt

        ! Constant Parameters
        real, parameter :: 
     .  zlmin = 0.1,
     .  zeps  = 1.e-10

        ! Input Variables
        real, dimension(gr%nnzp), intent(in) :: 
     .  thvm,   ! Virtual potential temp. on themodynamic level   [K]
     .  thlm,   ! Liquid potential temp. on themodynamic level    [K]
     .  rtm,    ! Total water mixing ratio on themodynamic level  [kg/kg]
     .  rcm,    ! Cloud water mixing ration on themodynamic level [kg/kg]
     .  em,     ! em = 3/2 * w'^2; on momentum level              [m^2/s^2]
     .  exner,  ! Exner function on thermodynamic level           [-]
     .  p       ! Pressure on thermodynamic level                 [Pa]

        ! Output Variables
        real, dimension(gr%nnzp), intent(out) :: 
     .  Lscale  ! Mixing length                 [m]
        
        integer, intent(inout) ::
     .  err_code

        ! Local Variables
        integer :: i, j

        real :: tke_i, CAPE_incr, dthv

        real :: zincr

        ! Minimum value for Lscale that will taper off with height
        real :: lminh

        ! Parcel quantities at grid level j
        real :: thl_par_j, rt_par_j, rc_par_j, thv_par_j

        ! Used in latent heating calculation
        real :: tl_par_j, rsl_par_j, beta_par_j,
     .    s_par_j, rsat_j

        ! Parcel quantities at grid level j-1
        real :: thl_par_j_minus_1, rt_par_j_minus_1, rc_par_j_minus_1 

        ! Parcel quantities at grid level j+1
        real :: thl_par_j_plus_1, rt_par_j_plus_1, rc_par_j_plus_1 

        ! Variables to make L nonlocal
        real :: Lup_max_alt, Ldown_min_alt

!---------- Mixing length computation ----------------------------------

        ! Avoid uninitialized memory (these values are not used in Lscale) 
        ! -dschanen 12 March 2008
        Lup(1)   = 0.0
        Ldown(1) = 0.0

        ! Upwards loop

        Lup_max_alt = 0.
        do i=2,gr%nnzp

          tke_i = zm2zt( em, i )                   ! tke on thermo level

          Lup(i) = zlmin
          j = i + 1

          thl_par_j_minus_1 = thlm(i)
          rt_par_j_minus_1  = rtm(i)
          rc_par_j_minus_1  = rcm(i)

          do while ((tke_i > 0.) .and. (j < gr%nnzp))

            ! thl, rt of parcel are conserved except for entrainment

            thl_par_j = ( 1 - mu/gr%dzm(j-1) ) * thl_par_j_minus_1
     .                + ( mu/gr%dzm(j-1) ) * thlm(j)

            rt_par_j = ( 1 - mu/gr%dzm(j-1) ) * rt_par_j_minus_1
     .               + ( mu/gr%dzm(j-1) ) * rtm(j)

!           Include effects of latent heating on Lup 6/12/00
!           Use thermodynamic formula of Bougeault 1981 JAS Vol. 38, 2416
!           Probably should use properties of bump 1 in Gaussian, not mean!!!

            tl_par_j = thl_par_j*exner(j)
            rsl_par_j = rsat(p(j),tl_par_j)
            ! SD's beta (eqn. 8)
            beta_par_j = ep*(Lv/(Rd*tl_par_j))*(Lv/(cp*tl_par_j))
            ! s from Lewellen and Yoh 1993 (LY) eqn. 1
            s_par_j = (rt_par_j-rsl_par_j)/(1+beta_par_j*rsl_par_j)
            rc_par_j = max( s_par_j, 0. )

            ! theta_v of entraining parcel
            thv_par_j = thl_par_j + ep1 * T0 * rt_par_j
     .       + ( Lv / (exner(j)*cp) - ep2 * T0 ) * rc_par_j

            CAPE_incr = ( ( grav/thvm(j) ) / gr%dzm(j-1) ) 
     .                  * ( thv_par_j - thvm(j) ) 

            if (tke_i+CAPE_incr > 0.) then
              Lup(i) = Lup(i) + gr%zt(j) - gr%zt(j-1)
            else
              Lup(i) = Lup(i) + ( gr%zt(j) - gr%zt(j-1) ) 
     .                          * tke_i/max( zeps, -CAPE_incr )
            end if

            thl_par_j_minus_1 = thl_par_j
            rt_par_j_minus_1 = rt_par_j
            rc_par_j_minus_1 = rc_par_j

            tke_i = tke_i + CAPE_incr
            j = j + 1

          end do

          ! Make Lup nonlocal

          Lup_max_alt = max( Lup_max_alt, Lup(i)+gr%zt(i) )
          if ( ( gr%zt(i) + Lup(i) ) < Lup_max_alt ) then
            Lup(i) = Lup_max_alt - gr%zt(i)
          end if

        end do

        ! Do it again for downwards particle motion.
        ! For now, do not include latent heat 

        ! Chris Golaz modification to include effects on latent heating
        ! on Ldown

        Ldown_min_alt = gr%zt(gr%nnzp)
        do i=gr%nnzp,2,-1

          tke_i = zm2zt( em, i )  ! tke on thermo level

          Ldown(i) = zlmin
          j = i - 1

          thl_par_j_plus_1 = thlm(i)
          rt_par_j_plus_1 = rtm(i)
          rc_par_j_plus_1 = rcm(i)

          do while ( (tke_i > 0.) .and. (j >= 2) )

            ! thl, rt of parcel are conserved except for entrainment

            thl_par_j = ( 1 - mu/gr%dzm(j) ) * thl_par_j_plus_1
     .                +  ( mu/gr%dzm(j) ) * thlm(j)

            rt_par_j = ( 1 - mu/gr%dzm(j) ) * rt_par_j_plus_1
     .               +  ( mu/gr%dzm(j) ) * rtm(j)

           ! Include effects of latent heating on Ldown
           ! Use thermodynamic formula of Bougeault 1981 JAS Vol. 38, 2416
           ! Probably should use properties of bump 1 in Gaussian, not mean!!!

            tl_par_j = thl_par_j*exner(j)
            rsl_par_j = rsat(p(j),tl_par_j)
            ! SD's beta (eqn. 8)
            beta_par_j = ep*(Lv/(Rd*tl_par_j))*(Lv/(cp*tl_par_j))
            ! s from Lewellen and Yoh 1993 (LY) eqn. 1
            s_par_j = (rt_par_j-rsl_par_j)/(1+beta_par_j*rsl_par_j)
            rc_par_j = max( s_par_j, 0. )

            ! theta_v of entraining parcel
            thv_par_j = thl_par_j + ep1 * T0 * rt_par_j
     .       + ( Lv / (exner(j)*cp) - ep2 * T0 ) * rc_par_j

           ! New code: CAPE_incr including moisture effects

           CAPE_incr = -( ( grav/thvm(j) ) / gr%dzm(j) ) 
     .          * ( thv_par_j - thvm(j) ) 

           ! Old code: CAPE_incr without including moisture effects

!            CAPE_incr = - grav/thvm(j) * (thvm(i)-thvm(j)) / gr%dzm(j)

            if (tke_i+CAPE_incr > 0.) then
              Ldown(i) = Ldown(i) + gr%zt(j+1) - gr%zt(j)
            else
              Ldown(i) = Ldown(i) + ( gr%zt(j+1) - gr%zt(j) ) 
     .                              * tke_i/max( zeps, -CAPE_incr )
            end if

            ! Bug fix from Brian 1/25/04: missing update

            thl_par_j_plus_1 = thl_par_j
            rt_par_j_plus_1  = rt_par_j
            rc_par_j_plus_1  = rc_par_j

            tke_i = tke_i + CAPE_incr
            j = j - 1

          end do

          ! Make Ldown nonlocal
!         Ldown_min_alt = max( Ldown_min_alt, gr%zt(i)-Ldown(i) )
          Ldown_min_alt = min( Ldown_min_alt, gr%zt(i)-Ldown(i) ) ! %% test
          if ( (gr%zt(i)-Ldown(i)) > Ldown_min_alt ) then
            Ldown(i) = gr%zt(i) - Ldown_min_alt
          end if

        end do

        do i=2,gr%nnzp
          ! Make lminh a linear function starting at value lmin at the
          ! bottom and going to zero at 500 meters in altitude.
          ! -dschanen 27 April 2007
          lminh = max( 0., 500. - gr%zt(i) ) * ( lmin / 500. )

          Lup(i)    = max( lminh, Lup(i) )
          Ldown(i)  = max( lminh, Ldown(i) )

          Lscale(i) = sqrt( Lup(i)*Ldown(i) )

        end do

        Lscale(1) = Lscale(2)
        Lscale(gr%nnzp) = Lscale(gr%nnzp-1)

        ! Vince Larson limited Lscale to allow host model to take over deep convection.  13 Feb 2008.
        !Lscale = min( Lscale, 1e5 )
        Lscale = min( Lscale, Lscale_max )
        
        if( clubb_at_debug_level( 2 ) ) then
                
                ! Ensure that the output from this subroutine is valid.
                call length_check( Lscale, Lup, Ldown, err_code )
                ! Joshua Fasching January 2008

        !       Error Reporting
        !       Joshua Fasching February 2008
                
        !       isValid replaced with err_code
        !       Joshua Fasching March 2008
                if ( err_code == clubb_var_equals_NaN ) then
                        
                   write(fstderr,*) "Errors in length subroutine"
                   
                   write(fstderr,*) "Intent(in)"
                   
                   write(fstderr,*) "thvm = ", thvm
                   write(fstderr,*) "thlm = ", thlm
                   write(fstderr,*) "rtm = ", rtm
                   write(fstderr,*) "rcm = ", rcm
                   write(fstderr,*) "em = ", em
                   write(fstderr,*) "exner = ", exner
                   write(fstderr,*) "p = ", p
                   
                   write(fstderr,*) "Intent(out)"

                   write(fstderr,*) "Lscale = ", Lscale
                   write(fstderr,*) "Lup = ", Lup
                   
                endif ! err_code == clubb_var_equals_NaN

        endif ! clubb_debug_level

        return

        end subroutine compute_length

        end module mixing_length
