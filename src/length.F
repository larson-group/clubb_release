!-----------------------------------------------------------------------------
! Subroutine length
!
! Larson's 5th moist, nonlocal length scale
!

        subroutine length(thvm, thlm, rtm, rlm, em, p, exner,
     .                    Lscale, lup, ldown, mu)

        use constants
        use grid_class

        implicit none

! External

        real, external :: rsat

! Input

        real, dimension(1:gr%nnzp), intent(in) :: thvm    ! themodynamic level
        real, dimension(1:gr%nnzp), intent(in) :: thlm    ! themodynamic level
        real, dimension(1:gr%nnzp), intent(in) :: rtm     ! themodynamic level
        real, dimension(1:gr%nnzp), intent(in) :: rlm     ! themodynamic level
        real, dimension(1:gr%nnzp), intent(in) :: em      ! momentum level
        real, dimension(1:gr%nnzp), intent(in) :: exner
        real, dimension(1:gr%nnzp), intent(in) :: p
!       mu = fractional entrainment rate per unit altitude [1/m].
!       mu = (1/M)dM/dz > 0.  mu=0 for no entrainment.
!       Siebesma recommends mu=2e-3, although most schemes use mu=1e-4
!       When mu was fixed, we used the value mu = 6.e-4
        real, intent(in) :: mu

! Output

        real, dimension(1:gr%nnzp), intent(out) :: Lscale
        real, dimension(1:gr%nnzp), intent(out) :: lup
        real, dimension(1:gr%nnzp), intent(out) :: ldown

! Constants

        real, parameter :: zlmin = 0.1
        real, parameter :: zeps = 1.e-10

! Internal

        integer i,j

        real tke_i,CAPE_incr,dthv

        real zincr

!       parcel quantities at grid level j
        real thl_par_j, rt_par_j, rl_par_j, thv_par_j
!       used in latent heating calculation
        real tl_par_j, rsl_par_j, beta_par_j
        real s_par_j, rsat_j
!       parcel quantities at grid level j-1
        real thl_par_j_minus_1, rt_par_j_minus_1, rl_par_j_minus_1 
!       parcel quantities at grid level j+1
        real thl_par_j_plus_1, rt_par_j_plus_1, rl_par_j_plus_1 

!       variables to make L nonlocal
        real lup_max_alt, ldown_min_alt



! Mixing length computation

!       upwards loop

        lup_max_alt = 0.
        do i=2,gr%nnzp

          tke_i = 0.5 * ( em(i-1) + em(i) )        ! tke on thermo level

          lup(i) = zlmin
          j = i + 1

          thl_par_j_minus_1 = thlm(i)
          rt_par_j_minus_1 = rtm(i)
          rl_par_j_minus_1 = rlm(i)

          do while ((tke_i.gt.0.).and.(j.lt.gr%nnzp))

!           thl, rt of parcel are conserved except for entrainment

            thl_par_j = ( 1 - mu/gr%dzm(j-1) ) * thl_par_j_minus_1
     .                + ( mu/gr%dzm(j-1) ) * thlm(j)

            rt_par_j = ( 1 - mu/gr%dzm(j-1) ) * rt_par_j_minus_1
     .               + ( mu/gr%dzm(j-1) ) * rtm(j)

!           Include effects of latent heating on lup 6/12/00
!           Use thermodynamic formula of Bougeault 1981 JAS Vol. 38, 2416
!           Probably should use properties of bump 1 in Gaussian, not mean!!!

            tl_par_j = thl_par_j*exner(j)
            rsl_par_j = rsat(p(j),tl_par_j)
!           SD's beta (eqn. 8)
            beta_par_j = ep*(Lv/(R*tl_par_j))*(Lv/(cp*tl_par_j))
!           s from Lewellen and Yoh 1993 (LY) eqn. 1
            s_par_j = (rt_par_j-rsl_par_j)/(1+beta_par_j*rsl_par_j)
            rl_par_j = max(s_par_j,0.)

!           theta_v of entraining parcel
            thv_par_j = thl_par_j + ep1 * T0 * rt_par_j
     .       + ( Lv / (exner(j)*cp) - ep2 * T0 ) * rl_par_j

            CAPE_incr = ( ( g/thvm(j) ) / gr%dzm(j-1) ) 
     .                  * ( thv_par_j - thvm(j) ) 

            if (tke_i+CAPE_incr.gt.0.) then
              lup(i) = lup(i) + gr%zt(j) - gr%zt(j-1)
            else
              lup(i) = lup(i) + ( gr%zt(j) - gr%zt(j-1) ) 
     .                          * tke_i/max(zeps,-CAPE_incr)
            end if

            thl_par_j_minus_1 = thl_par_j
            rt_par_j_minus_1 = rt_par_j
            rl_par_j_minus_1 = rl_par_j

            tke_i = tke_i + CAPE_incr
            j = j + 1

          end do

!         make lup nonlocal

          lup_max_alt = max(lup_max_alt,lup(i)+gr%zt(i))
          if ( (gr%zt(i)+lup(i)).lt.lup_max_alt ) then
            lup(i) = lup_max_alt - gr%zt(i)
          end if

        end do

!       do it again for downwards particle motion
!       for now, do not include latent heat 

!       Chris Golaz modification to include effects on latent heating
!       on ldown

        ldown_min_alt = gr%zt(gr%nnzp)
        do i=gr%nnzp,2,-1

          tke_i = 0.5 * ( em(i-1) + em(i) )        ! tke on thermo level

          ldown(i) = zlmin
          j = i - 1

          thl_par_j_plus_1 = thlm(i)
          rt_par_j_plus_1 = rtm(i)
          rl_par_j_plus_1 = rlm(i)

          do while ((tke_i.gt.0.).and.(j.ge.2))

!           thl, rt of parcel are conserved except for entrainment

            thl_par_j = ( 1 - mu/gr%dzm(j) ) * thl_par_j_plus_1
     .                +  ( mu/gr%dzm(j) ) * thlm(j)

            rt_par_j = ( 1 - mu/gr%dzm(j) ) * rt_par_j_plus_1
     .               +  ( mu/gr%dzm(j) ) * rtm(j)

!           Include effects of latent heating on ldown
!           Use thermodynamic formula of Bougeault 1981 JAS Vol. 38, 2416
!           Probably should use properties of bump 1 in Gaussian, not mean!!!

            tl_par_j = thl_par_j*exner(j)
            rsl_par_j = rsat(p(j),tl_par_j)
!           SD's beta (eqn. 8)
            beta_par_j = ep*(Lv/(R*tl_par_j))*(Lv/(cp*tl_par_j))
!           s from Lewellen and Yoh 1993 (LY) eqn. 1
            s_par_j = (rt_par_j-rsl_par_j)/(1+beta_par_j*rsl_par_j)
            rl_par_j = max(s_par_j,0.)

!           theta_v of entraining parcel
            thv_par_j = thl_par_j + ep1 * T0 * rt_par_j
     .       + ( Lv / (exner(j)*cp) - ep2 * T0 ) * rl_par_j

!           New code: CAPE_incr including moisture effects

           CAPE_incr = -( ( g/thvm(j) ) / gr%dzm(j) ) 
     .          * ( thv_par_j - thvm(j) ) 


!           Old code: CAPE_incr without including moisture effects

!            CAPE_incr = - g/thvm(j) * (thvm(i)-thvm(j)) / gr%dzm(j)

            if (tke_i+CAPE_incr.gt.0.) then
              ldown(i) = ldown(i) + gr%zt(j+1) - gr%zt(j)
            else
              ldown(i) = ldown(i) + ( gr%zt(j+1) - gr%zt(j) ) 
     .                              * tke_i/max(zeps,-CAPE_incr)
            end if

!           Bug fix from Brian 1/25/04: missing update

            thl_par_j_plus_1 = thl_par_j
            rt_par_j_plus_1 = rt_par_j
            rl_par_j_plus_1 = rl_par_j

            tke_i = tke_i + CAPE_incr
            j = j - 1

          end do

!         make ldown nonlocal
          ldown_min_alt = max(ldown_min_alt,gr%zt(i)-ldown(i))
          if ( (gr%zt(i)-ldown(i)).gt.ldown_min_alt ) then
            ldown(i) = gr%zt(i) - ldown_min_alt
          end if

        end do

        do i=2,gr%nnzp

          lup(i)=max(lmin,lup(i))
          ldown(i)=max(lmin,ldown(i))
          Lscale(i) = sqrt(lup(i)*ldown(i))

        end do

        Lscale(1) = Lscale(2)
        Lscale(gr%nnzp) = Lscale(gr%nnzp-1)

        return
        end subroutine length
