!------------------------------------------------------------------------
! $Id$
!===============================================================================
module advance_windm_edsclrm_module

  implicit none

  private ! Set Default Scope

  public :: advance_windm_edsclrm

  private :: windm_edsclrm_solve, &
    compute_uv_tndcy,  &
    windm_edsclrm_lhs, &
    windm_edsclrm_rhs, &
    xpwp_fnc

  contains

  !=============================================================================
  subroutine advance_windm_edsclrm &
             ( dt, wm_zt, Kh_zm, ug, vg, um_ref, vm_ref, &
               wp2, up2, vp2, um_forcing, vm_forcing, &
               edsclrm_forcing, upwp_sfc, vpwp_sfc, &
               wpedsclrp_sfc, fcor, l_implemented, &
               um, vm, edsclrm, &
               upwp, vpwp, wpedsclrp, err_code )
    ! Description:
    ! Solves for both mean horizontal wind components, um and vm, and for the
    ! eddy-scalars (passive scalars that don't use the high-order closure).

    ! Uses the LAPACK tridiagonal solver subroutine with 2 + # of scalar(s)
    ! back substitutions (since the left hand side matrix is the same for all
    ! input variables).

    ! References:
    ! Eqn. 8 & 9 on p. 3545 of
    ! ``A PDF-Based Model for Boundary Layer Clouds. Part I:
    ! Method and Model Description'' Golaz, et al. (2002)
    ! JAS, Vol. 59, pp. 3540--3551.
    !-----------------------------------------------------------------------

    use grid_class, only:  &
      gr  ! Variables(s)

    use parameters_model, only:  &
      ts_nudge,  & ! Variable(s)
      edsclr_dim

    use model_flags, only:  &
      l_uv_nudge,  & ! Variable(s)
      l_tke_aniso

    use stats_precision, only:  &
      time_precision  ! Variable(s)

    use stats_type, only: &
      stat_begin_update, & ! Subroutines
      stat_end_update, &
      stat_update_var

    use stats_variables, only: &
      ium_ref, &
      ivm_ref, &
      ivm_bt, & ! Variables
      ium_bt, &
      ium_sdmp, &
      ivm_sdmp, &
      iwindm_matrix_condt_num, &
      zt,     &
      l_stats_samp

    use clip_explicit, only:  &
      clip_covariance  ! Procedure(s)

    use error_code, only:  & 
      lapack_error,  & ! Procedure(s)
      clubb_at_least_debug_level

    use constants, only:  & 
        fstderr, &  ! Constant
        eps

    use sponge_layer_damping, only: &
      uv_sponge_damp_settings, &  
      uv_sponge_damp_profile, &  
      sponge_damp_xm ! Procedure(s)

    implicit none

    ! Input Variables
    real(kind=time_precision), intent(in) ::  &
      dt             ! Model timestep                                [s]

    real, dimension(gr%nnzp), intent(in) ::  &
      wm_zt,       & ! w wind component on thermodynamic levels      [m/s]
      Kh_zm,       & ! Eddy diffusivity on momentum levels           [m^2/s]
      ug,          & ! u (west-to-east) geostrophic wind component   [m/s]
      vg,          & ! v (south-to-north) geostrophic wind component [m/s]
      um_ref,      & ! Reference u wind component for nudging        [m/s]
      vm_ref,      & ! Reference v wind component for nudging        [m/s]
      wp2,         & ! w'^2 (momentum levels)                        [m^2/s^2]
      up2,         & ! u'^2 (momentum levels)                        [m^2/s^2]
      vp2,         & ! v'^2 (momentum levels)                        [m^2/s^2]
      um_forcing,  & ! u forcing                                     [m/s/s]
      vm_forcing     ! v forcing                                     [m/s/s]

    real, dimension(gr%nnzp,edsclr_dim), intent(in) ::  &
      edsclrm_forcing  ! Eddy scalar large-scale forcing             [{units vary}/s]

    real, intent(in) ::  &
      upwp_sfc,    & ! u'w' at the surface (momentum level 1)        [m^2/s^2]
      vpwp_sfc,    & ! v'w' at the surface (momentum level 1)        [m^2/s^2]
      fcor           ! Coriolis parameter                            [s^-1]
    real, dimension(edsclr_dim), intent(in) :: &
      wpedsclrp_sfc  ! w'edsclr' at the surface (momentum level 1)    [units vary]

    logical, intent(in) ::  &
      l_implemented  ! Flag for CLUBB being implemented in a larger model.

    ! Input/Output Variables
    real, dimension(gr%nnzp), intent(inout) ::  &
      um,          & ! Mean u (west-to-east) wind component          [m/s]
      vm             ! Mean v (south-to-north) wind component        [m/s]

    ! Input/Output Variable for eddy-scalars
    real, dimension(gr%nnzp,edsclr_dim), intent(inout) ::  &
      edsclrm        ! Mean eddy scalar quantity                     [units vary]

    ! Output Variables
    real, dimension(gr%nnzp), intent(out) ::  &
      upwp,        & ! u'w' (momentum levels)                        [m^2/s^2]
      vpwp           ! v'w' (momentum levels)                        [m^2/s^2]

    ! Output Variable for eddy-scalars
    real, dimension(gr%nnzp,edsclr_dim), intent(out) ::  &
      wpedsclrp      ! w'edsclr' (momentum levels)                   [units vary]

    integer, intent(out) :: &
      err_code       ! clubb_singular_matrix when matrix is singular

    ! Local Variables
    real, dimension(gr%nnzp) ::  &
      um_tndcy,    & ! u wind component tendency                     [m/s^2]
      vm_tndcy       ! v wind component tendency                     [m/s^2]

    real, dimension(gr%nnzp) ::  &
      upwp_chnge,  & ! Net change of u'w' due to clipping            [m^2/s^2]
      vpwp_chnge     ! Net change of v'w' due to clipping            [m^2/s^2]

    real, dimension(3,gr%nnzp) :: &
      lhs ! The implicit part of the tridiagonal matrix              [units vary]

    real, dimension(gr%nnzp,max(2,edsclr_dim)) :: &
      rhs,     &! The explicit part of the tridiagonal matrix        [units vary]
      solution  ! The solution to the tridiagonal matrix             [units vary]

    real, dimension(gr%nnzp) :: &
      wind_speed  ! wind speed; sqrt(u^2 + v^2)                      [m/s]

    real :: &
      u_star_sqd  ! Surface friction velocity, u_star, squared       [m/s]

    logical :: &
      l_imp_sfc_momentum_flux  ! Flag for implicit momentum surface fluxes.

    integer :: i     ! Array index

    !--------------------------- Begin Code ------------------------------------

    if ( l_stats_samp ) then

      ! xm total time tendency (1st calculation)
      call stat_begin_update( ium_bt, real( um / dt ), zt )

      call stat_begin_update( ivm_bt, real( vm / dt ), zt )

    endif

    !----------------------------------------------------------------
    ! Prepare tridiagonal system for horizontal winds, um and vm
    !----------------------------------------------------------------

    ! Compute Coriolis, geostrophic, and other prescribed wind forcings for um.
    call compute_uv_tndcy( "um", fcor, vm, vg, um_forcing, l_implemented, & ! in
                           um_tndcy )                                       ! out

    ! Compute Coriolis, geostrophic, and other prescribed wind forcings for vm.
    call compute_uv_tndcy( "vm", fcor, um, ug, vm_forcing, l_implemented, & ! in
                           vm_tndcy )                                       ! out

    ! Momentum surface fluxes, u'w'|_sfc and v'w'|_sfc, are applied to through
    ! an implicit method, such that:
    !    x'w'|_sfc = - ( u_star(t)^2 / wind_speed(t) ) * xm(t+1).
    l_imp_sfc_momentum_flux = .true.
    ! Compute wind speed (use threshold "eps" to prevent divide-by-zero error).
    wind_speed = max( sqrt( um**2 + vm**2 ), eps )
    ! Compute u_star_sqd according to the definition of u_star.
    u_star_sqd = sqrt( upwp_sfc**2 + vpwp_sfc**2 )

    ! Compute the explicit portion of the um equation.
    ! Build the right-hand side vector.
    rhs(1:gr%nnzp,1) = windm_edsclrm_rhs( "um", dt, Kh_zm, um, um_tndcy,  &   ! in
                                          l_imp_sfc_momentum_flux, upwp_sfc ) ! in

    ! Compute the explicit portion of the vm equation.
    ! Build the right-hand side vector.
    rhs(1:gr%nnzp,2) = windm_edsclrm_rhs( "vm", dt, Kh_zm, vm, vm_tndcy,  &   ! in
                                          l_imp_sfc_momentum_flux, vpwp_sfc ) ! in


    ! Store momentum flux (explicit component)

    ! The surface flux, x'w'(1) = x'w'|_sfc, is set elsewhere in the model.
    upwp(1) = upwp_sfc
    vpwp(1) = vpwp_sfc

    ! Solve for x'w' at all intermediate model levels.
    ! A Crank-Nicholson timestep is used.

    upwp(2:gr%nnzp-1) = - 0.5 * xpwp_fnc( Kh_zm(2:gr%nnzp-1), um(2:gr%nnzp-1), & ! in
                                          um(3:gr%nnzp), gr%dzm(2:gr%nnzp-1) )   ! in

    vpwp(2:gr%nnzp-1) = - 0.5 * xpwp_fnc( Kh_zm(2:gr%nnzp-1), vm(2:gr%nnzp-1), & ! in
                                          vm(3:gr%nnzp), gr%dzm(2:gr%nnzp-1) )   ! in

    ! A zero-flux boundary condition at the top of the model, d(xm)/dz = 0,
    ! means that x'w' at the top model level is 0,
    ! since x'w' = - K_zm * d(xm)/dz.
    upwp(gr%nnzp) = 0.
    vpwp(gr%nnzp) = 0.


    ! Compute the implicit portion of the um and vm equations.
    ! Build the left-hand side matrix.
    call windm_edsclrm_lhs( dt, wm_zt, Kh_zm, wind_speed, u_star_sqd,  & ! in
                            l_implemented, l_imp_sfc_momentum_flux,  &   ! in
                            lhs )                                        ! out

    ! Decompose and back substitute for um and vm
    call windm_edsclrm_solve( 2, iwindm_matrix_condt_num, & ! in
                              lhs, rhs, &                   ! in/out
                              solution, err_code )          ! out

    !----------------------------------------------------------------
    ! Update zonal (west-to-east) component of mean wind, um
    !----------------------------------------------------------------
    um(1:gr%nnzp) = solution(1:gr%nnzp,1)

    !----------------------------------------------------------------
    ! Update meridional (south-to-north) component of mean wind, vm
    !----------------------------------------------------------------
    vm(1:gr%nnzp) = solution(1:gr%nnzp,2)

    ! The values of um(1) and vm(1) are located below the model surface and do
    ! not effect the rest of the model.  The values of um(1) or vm(1) are simply
    ! set to the values of um(2) and vm(2), respectively, after the equation
    ! matrices has been solved.  Even though um and vm would sharply decrease
    ! to a value of 0 at the surface, this is done to avoid confusion on plots
    ! of the vertical profiles of um and vm.
    um(1) = um(2)
    vm(1) = vm(2)


    if ( uv_sponge_damp_settings%l_sponge_damping ) then
      if( l_stats_samp ) then
        call stat_begin_update( ium_sdmp, real( um/dt ), zt )
        call stat_begin_update( ivm_sdmp, real( vm/dt ), zt )
      endif

      um(1:gr%nnzp) = sponge_damp_xm( dt, um_ref(1:gr%nnzp), um(1:gr%nnzp), &
                                      uv_sponge_damp_profile )
      vm(1:gr%nnzp) = sponge_damp_xm( dt, vm_ref(1:gr%nnzp), vm(1:gr%nnzp), &
                                      uv_sponge_damp_profile )
      if( l_stats_samp ) then
        call stat_end_update( ium_sdmp, real( um/dt ), zt )
        call stat_end_update( ivm_sdmp, real( vm/dt ), zt )
      endif

    endif


    if ( l_stats_samp ) then

      ! xm total time tendency (2nd calculation)
      call stat_end_update( ium_bt, real( um / dt ), zt ) ! in

      call stat_end_update( ivm_bt, real( vm / dt ), zt ) ! in

      ! Implicit contributions to um and vm
      call windm_edsclrm_implicit_stats( "um", um ) ! in

      call windm_edsclrm_implicit_stats( "vm", vm ) ! in

    endif ! l_stats_samp


    ! Second part of momentum (implicit component)

    ! Solve for x'w' at all intermediate model levels.
    ! A Crank-Nicholson timestep is used.

    upwp(2:gr%nnzp-1) = upwp(2:gr%nnzp-1) &
      - 0.5 * xpwp_fnc( Kh_zm(2:gr%nnzp-1), um(2:gr%nnzp-1), um(3:gr%nnzp), gr%dzm(2:gr%nnzp-1) )!in

    vpwp(2:gr%nnzp-1) = vpwp(2:gr%nnzp-1) &
      - 0.5 * xpwp_fnc( Kh_zm(2:gr%nnzp-1), vm(2:gr%nnzp-1), vm(3:gr%nnzp), gr%dzm(2:gr%nnzp-1) )!in


    ! Adjust um and vm if nudging is turned on.
    if ( l_uv_nudge ) then
      um(1:gr%nnzp) = real( um(1:gr%nnzp) - ((um(1:gr%nnzp) - um_ref(1:gr%nnzp)) * (dt/ts_nudge)) )
      vm(1:gr%nnzp) = real( vm(1:gr%nnzp) - ((vm(1:gr%nnzp) - vm_ref(1:gr%nnzp)) * (dt/ts_nudge)) )
    endif

    if( l_stats_samp ) then
      call stat_update_var(ium_ref, um_ref, zt)
      call stat_update_var(ivm_ref, vm_ref, zt)
    end if

    if ( l_tke_aniso ) then

      ! Clipping for u'w'
      !
      ! Clipping u'w' at each vertical level, based on the
      ! correlation of u and w at each vertical level, such that:
      ! corr_(u,w) = u'w' / [ sqrt(u'^2) * sqrt(w'^2) ];
      ! -1 <= corr_(u,w) <= 1.
      !
      ! Since u'^2, w'^2, and u'w' are each advanced in different subroutines from
      ! each other in advance_clubb_core, clipping for u'w' has to be done three
      ! times during each timestep (once after each variable has been updated).
      ! This is the third instance of u'w' clipping.
      call clip_covariance( "upwp", .false.,      & ! intent(in)
                            .true., dt, wp2, up2, & ! intent(in)
                            upwp, upwp_chnge )      ! intent(inout)

      ! Clipping for v'w'
      !
      ! Clipping v'w' at each vertical level, based on the
      ! correlation of v and w at each vertical level, such that:
      ! corr_(v,w) = v'w' / [ sqrt(v'^2) * sqrt(w'^2) ];
      ! -1 <= corr_(v,w) <= 1.
      !
      ! Since v'^2, w'^2, and v'w' are each advanced in different subroutines from
      ! each other in advance_clubb_core, clipping for v'w' has to be done three
      ! times during each timestep (once after each variable has been updated).
      ! This is the third instance of v'w' clipping.
      call clip_covariance( "vpwp", .false.,      & ! intent(in)
                            .true., dt, wp2, vp2, & ! intent(in)
                            vpwp, vpwp_chnge )      ! intent(inout)

    else

      ! In this case, it is assumed that
      !   u'^2 == v'^2 == w'^2, and the variables `up2' and `vp2' do not interact with
      !   any other variables.

      call clip_covariance( "upwp", .false.,      & ! intent(in)
                            .true., dt, wp2, wp2, & ! intent(in)
                            upwp, upwp_chnge )      ! intent(inout)

      call clip_covariance( "vpwp", .false.,      & ! intent(in)
                            .true., dt, wp2, wp2, & ! intent(in)
                            vpwp, vpwp_chnge )      ! intent(inout)

    endif ! l_tke_aniso


    !----------------------------------------------------------------
    ! Prepare tridiagonal system for eddy-scalars
    !----------------------------------------------------------------

    if ( edsclr_dim > 0 ) then

      ! Eddy-scalar surface fluxes, x'w'|_sfc, are applied through an explicit
      ! method.
      l_imp_sfc_momentum_flux = .false.

      ! Compute the explicit portion of eddy scalar equation.
      ! Build the right-hand side vector.
      ! Because of statistics, we have to use a DO rather than a FORALL here
      ! -dschanen 7 Oct 2008
!HPF$ INDEPENDENT
      do i = 1, edsclr_dim
        rhs(1:gr%nnzp,i)  &
        = windm_edsclrm_rhs( "scalars", dt, Kh_zm, edsclrm(:,i), edsclrm_forcing,  & ! in
                             l_imp_sfc_momentum_flux, wpedsclrp_sfc(i) )     ! in
      enddo


      ! Store momentum flux (explicit component)

      ! The surface flux, x'w'(1) = x'w'|_sfc, is set elsewhere in the model.
      wpedsclrp(1,1:edsclr_dim) =  wpedsclrp_sfc(1:edsclr_dim)

      ! Solve for x'w' at all intermediate model levels.
      ! A Crank-Nicholson timestep is used.
      ! Here we use a forall and high performance fortran directive to try to
      ! parallelize this computation.  Note that FORALL is more restrictive than DO.
!HPF$ INDEPENDENT, REDUCTION(wpedsclrp)
      forall( i = 1:edsclr_dim )
        wpedsclrp(2:gr%nnzp-1,i) = &
          - 0.5 * xpwp_fnc( Kh_zm(2:gr%nnzp-1), edsclrm(2:gr%nnzp-1,i), & ! in
                            edsclrm(3:gr%nnzp,i), gr%dzm(2:gr%nnzp-1) )   ! in
      end forall

      ! A zero-flux boundary condition at the top of the model, d(xm)/dz = 0,
      ! means that x'w' at the top model level is 0,
      ! since x'w' = - K_zm * d(xm)/dz.
      wpedsclrp(gr%nnzp,1:edsclr_dim) = 0.


      ! Compute the implicit portion of the xm (eddy-scalar) equations.
      ! Build the left-hand side matrix.
      call windm_edsclrm_lhs( dt, wm_zt, Kh_zm, wind_speed, u_star_sqd,  & ! in
                              l_implemented, l_imp_sfc_momentum_flux,  &   ! in
                              lhs )                                        ! out

      ! Decompose and back substitute for all eddy-scalar variables
      call windm_edsclrm_solve( edsclr_dim, 0, &     ! in
                                lhs, rhs, &          ! in/out
                                solution, err_code ) ! out

      !----------------------------------------------------------------
      ! Update Eddy-diff. Passive Scalars
      !----------------------------------------------------------------
      edsclrm(1:gr%nnzp,1:edsclr_dim) = solution(1:gr%nnzp,1:edsclr_dim)

      ! The value of edsclrm(1) is located below the model surface and does not
      ! effect the rest of the model.  The value of edsclrm(1) is simply set to
      ! the value of edsclrm(2) after the equation matrix has been solved.
      forall( i=1:edsclr_dim )
        edsclrm(1,i) = edsclrm(2,i)
      end forall

      ! Second part of momentum (implicit component)

      ! Solve for x'w' at all intermediate model levels.
      ! A Crank-Nicholson timestep is used.
!HPF$ INDEPENDENT, REDUCTION(wpedsclrp)
      forall( i = 1:edsclr_dim )
        wpedsclrp(2:gr%nnzp-1,i) = wpedsclrp(2:gr%nnzp-1,i) &
          - 0.5 * xpwp_fnc( Kh_zm(2:gr%nnzp-1), edsclrm(2:gr%nnzp-1,i), & ! in
                            edsclrm(3:gr%nnzp,i), gr%dzm(2:gr%nnzp-1) )   ! in
      end forall

      ! Note that the w'edsclr' terms are not clipped, since we don't compute the
      ! variance of edsclr anywhere. -dschanen 7 Oct 2008

    endif


    ! Error report
    ! Joshua Fasching February 2008
    if ( lapack_error( err_code ) .and.  &
         clubb_at_least_debug_level( 1 ) ) then

      write(fstderr,*) "Error in advance_windm_edsclrm"

      write(fstderr,*) "Intent(in)"

      write(fstderr,*) "dt = ", dt
      write(fstderr,*) "wm_zt = ", wm_zt
      write(fstderr,*) "Kh_zm = ", Kh_zm
      write(fstderr,*) "ug = ", ug
      write(fstderr,*) "vg = ", vg
      write(fstderr,*) "um_ref = ", um_ref
      write(fstderr,*) "vm_ref = ", vm_ref
      write(fstderr,*) "wp2 = ", wp2
      write(fstderr,*) "up2 = ", up2
      write(fstderr,*) "vp2 = ", vp2
      write(fstderr,*) "upwp_sfc = ", upwp_sfc
      write(fstderr,*) "vpwp_sfc = ", vpwp_sfc
      write(fstderr,*) "fcor = ", fcor
      write(fstderr,*) "l_implemented = ", l_implemented

      write(fstderr,*) "Intent(inout)"

      write(fstderr,*) "um = ", um
      write(fstderr,*) "vm = ", vm
      write(fstderr,*) "edsclrm = ", edsclrm

      write(fstderr,*) "Intent(out)"

      write(fstderr,*) "upwp = ", upwp
      write(fstderr,*) "vpwp = ", vpwp
      write(fstderr,*) "wpedsclrp = ", wpedsclrp

      return

    endif

    return
  end subroutine advance_windm_edsclrm

  !=============================================================================
  subroutine windm_edsclrm_solve( nrhs, ixm_matrix_condt_num, &
                                  lhs, rhs, solution, err_code )

    ! Description:
    ! Solves the horizontal wind or eddy-scalar time-tendency equation, and
    ! diagnoses the turbulent flux.  A Crank-Nicholson time-stepping algorithm
    ! is used in solving the turbulent advection term and in diagnosing the
    ! turbulent flux.
    !
    ! The rate of change of an eddy-scalar quantity, xm, is:
    !
    ! d(xm)/dt = - w * d(xm)/dz - d(x'w')/dz + xm_forcings.
    !
    !
    ! The Turbulent Advection Term
    ! ----------------------------
    !
    ! The above equation contains a turbulent advection term:
    !
    ! - d(x'w')/dz;
    !
    ! where the momentum flux, x'w', is closed using a down gradient approach:
    !
    ! x'w' = - K_zm * d(xm)/dz.
    !
    ! The turbulent advection term becomes:
    !
    ! + d [ K_zm * d(xm)/dz ] / dz;
    !
    ! which is the same as an eddy-diffusion term.  Thus, the turbulent
    ! advection term is treated and solved in the same way that an
    ! eddy-diffusion term would be solved.  The term is discretized as follows:
    !
    ! The values of xm are found on the thermodynamic levels, while the values
    ! of K_zm are found on the momentum levels.  The derivatives (d/dz) of xm
    ! are taken over the intermediate momentum levels.  At the intermediate
    ! momentum levels, d(xm)/dz is multiplied by K_zm.  Then, the derivative of
    ! the whole mathematical expression is taken over the central thermodynamic
    ! level, which yields the desired result.
    !
    ! ---xm(kp1)----------------------------------------------- t(k+1)
    !
    ! =============d(xm)/dz===K_zm(k)========================== m(k)
    !
    ! ---xm(k)---------------------------d[K_zm*d(xm)/dz]/dz--- t(k)
    !
    ! =============d(xm)/dz===K_zm(km1)======================== m(k-1)
    !
    ! ---xm(km1)----------------------------------------------- t(k-1)
    !
    ! The vertical indices t(k+1), m(k), t(k), m(k-1), and t(k-1) correspond
    ! with altitudes zt(k+1), zm(k), zt(k), zm(k-1), and zt(k-1), respectively.
    ! The letter "t" is used for thermodynamic levels and the letter "m" is used
    ! for momentum levels.
    !
    ! dzt(k)   = 1 / ( zm(k) - zm(k-1) )
    ! dzm(k)   = 1 / ( zt(k+1) - zt(k) )
    ! dzm(k-1) = 1 / ( zt(k) - zt(k-1) )
    !
    ! The vertically discretized form of the turbulent advection term (treated
    ! as an eddy diffusion term) is written out as:
    !
    ! + dzt(k) * [   K_zm(k) * dzm(k) * ( xm(k+1) - xm(k) )
    !              - K_zm(k-1) * dzm(k-1) * ( xm(k) - xm(k-1) ) ].
    !
    ! For this equation, a Crank-Nicholson (semi-implicit) diffusion scheme is
    ! used to solve the d [ K_zm * d(xm)/dz ] / dz eddy-diffusion term.  The
    ! discretized implicit portion of the term is written out as:
    !
    ! +(1/2)*dzt(k)*[  K_zm(k) * dzm(k) * ( xm(k+1,<t+1>) - xm(k,<t+1>) )
    !                - K_zm(k-1) * dzm(k-1) * ( xm(k,<t+1>) - xm(k-1,<t+1>) ) ].
    !
    ! Note:  When the implicit term is brought over to the left-hand side,
    !        the sign is reversed and the leading "+" in front of the term
    !        is changed to a "-".
    !
    ! The discretized explicit portion of the term is written out as:
    !
    ! + (1/2)*dzt(k)*[  K_zm(k) * dzm(k) * ( xm(k+1,<t>) - xm(k,<t>) )
    !                 - K_zm(k-1) * dzm(k-1) * ( xm(k,<t>) - xm(k-1,<t>) ) ].
    !
    ! Timestep index (t) stands for the index of the current timestep, while
    ! timestep index (t+1) stands for the index of the next timestep, which is
    ! being advanced to in solving the d(xm)/dt equation.
    !
    !
    ! Boundary Conditions:
    !
    ! An eddy-scalar quantity is not allowed to flux out the upper boundary.
    ! Thus, a zero-flux boundary condition is used for the upper boundary in the
    ! eddy-diffusion equation.
    !
    ! The lower boundary condition is much more complicated.  It is neither a
    ! zero-flux nor a fixed-point boundary condition.  Rather, it is a
    ! fixed-flux boundary condition.  This term is a turbulent advection term,
    ! but with the eddy-scalars, the only value of x'w' relevant in solving the
    ! d(xm)/dt equation is the value of x'w' at the surface (the first momentum
    ! level), which is written as x'w'|_sfc.
    !
    ! 1) x'w' surface flux; generalized explicit form
    !
    !    The x'w' surface flux is applied to the d(xm)/dt equation through the
    !    turbulent advection term, which is:
    !
    !    - d(x'w')/dz.
    !
    !    At most vertical levels, a substitution can be made for x'w', such
    !    that:
    !
    !    x'w' = - K_zm * d(xm)/dz.
    !
    !    However, the same substitution cannot be made at the surface (momentum
    !    level 1), as x'w'|_sfc is a surface flux that is explicitly computed
    !    elsewhere in the model code.
    !
    !    The lower boundary condition, which in this case needs to be applied to
    !    the d(xm)/dt equation at level 2, is discretized as follows:
    !
    !    ---xm(3)------------------------------------------------ t(3)
    !
    !    =============[x'w'(2) = -K_zm(2)*d(xm)/dz]============== m(2)
    !
    !    ---xm(2)--------------------------d[K_zm*d(xm)/dz]/dz--- t(2)
    !
    !    =============[x'w'|_sfc]================================ m(1) (surface)
    !
    !    ---xm(1)-------(below surface; not applicable)---------- t(1)
    !
    !    The vertically discretized form of the turbulent advection term
    !    (treated as an eddy diffusion term), with the explicit surface flux,
    !    x'w'|_sfc, in place, is written out as:
    !
    !    + dzt(2) * [ K_zm(2) * dzm(2) * ( xm(3) - xm(2) ) + x'w'|_sfc ];
    !
    !    which can be re-written as:
    !
    !    + dzt(2) * [ K_zm(2) * dzm(2) * ( xm(3) - xm(2) ) ]
    !       + dzt(2) * x'w'|_sfc.
    !
    !    For this equation, a Crank-Nicholson (semi-implicit) diffusion scheme
    !    is used to solve the d [ K_zm * d(xm)/dz ] / dz eddy-diffusion term.
    !    The discretized implicit portion of the term is written out as:
    !
    !    + (1/2)*dzt(2) * [ K_zm(2) * dzm(2) * ( xm(3,<t+1>) - xm(2,<t+1>) ) ].
    !
    !    Note:  When the implicit term is brought over to the left-hand side,
    !           the sign is reversed and the leading "+" in front of the term
    !           is changed to a "-".
    !
    !    The discretized explicit portion of the term is written out as:
    !
    !    + (1/2)*dzt(2) * [ K_zm(2) * dzm(2) * ( xm(3,<t>) - xm(2,<t>) ) ]
    !       + dzt(2) * x'w'|_sfc.
    !
    !    Note:  The x'w'|_sfc portion of the term written above has been pulled
    !           away from the rest of the explicit form written above because
    !           the (1/2) factor due to Crank-Nicholson time_stepping does not
    !           apply to it, as there isn't an implicit portion for x'w'|_sfc.
    !
    !    Timestep index (t) stands for the index of the current timestep, while
    !    timestep index (t+1) stands for the index of the next timestep, which
    !    is being advanced to in solving the d(xm)/dt equation.
    !
    ! 2) x'w' surface flux; implicit form for momentum fluxes u'w' and v'w'
    !
    !    The x'w' surface flux is applied to the d(xm)/dt equation through the
    !    turbulent advection term, which is:
    !
    !    - d(x'w')/dz.
    !
    !    At most vertical levels, a substitution can be made for x'w', such
    !    that:
    !
    !    x'w' = - K_zm * d(xm)/dz.
    !
    !    However, the same substitution cannot be made at the surface (momentum
    !    level 1), as x'w'|_sfc is a surface momentum flux that is found by the
    !    following equation:
    !
    !    x'w'|_sfc = - [ u_star^2 / sqrt( um^2 + vm^2 ) ] * xm;
    !
    !    where x'w'|_sfc and xm are either u'w'|_sfc and um or v'w'|_sfc and vm
    !    (um and vm are located at the first thermodynamic level above the
    !    surface, which is thermodynamic level 2), sqrt( um^2 + vm^2 ) is the
    !    wind speed (also at thermodynamic level 2), and u_star is defined as:
    !
    !    u_star = ( u'w'|_sfc^2 + v'w'|_sfc^2 )^(1/4);
    !
    !    and thus u_star^2 is defined as:
    !
    !    u_star^2 = sqrt( u'w'|_sfc^2 + v'w'|_sfc^2 ).
    !
    !    The value of u_star is either set to a constant value or computed
    !    (through function diag_ustar) based on the surface wind speed, the
    !    height above surface of the surface wind speed (as compared to the
    !    roughness height), and the buoyancy flux at the surface.  Either way,
    !    u_star is computed elsewhere in the model, and the values of u'w'|_sfc
    !    and v'w'|_sfc are based on it and computed along with it.  The values
    !    of u'w'|_sfc and v'w'|_sfc are then passed into advance_clubb_core,
    !    and are eventually passed into advance_windm_edsclrm.  In subroutine
    !    advance_windm_edsclrm, the value of u_star_sqd is then recomputed
    !    based on u'w'|_sfc and v'w'|_sfc.  The value of sqrt( u_star_sqd ) is
    !    consistent with the value of the original computation of u_star.
    !
    !    The equation listed above is substituted for x'w'|_sfc.  The lower
    !    boundary condition, which in this case needs to be applied to the
    !    d(xm)/dt equation at level 2, is discretized as follows:
    !
    !    ---xm(3)------------------------------------------------ t(3)
    !
    !    =========[x'w'(2) = -K_zm(2)*d(xm)/dz]================== m(2)
    !
    !    ---xm(2)--------------------------d[K_zm*d(xm)/dz]/dz--- t(2)
    !
    !    =========[x'w'|_sfc = -[u_star^2/sqrt(um^2+vm^2)]*xm]=== m(1) (surface)
    !
    !    ---xm(1)-------(below surface; not applicable)---------- t(1)
    !
    !    The vertically discretized form of the turbulent advection term
    !    (treated as an eddy diffusion term), with the implicit surface momentum
    !    flux in place, is written out as:
    !
    !    + dzt(2) * [ K_zm(2) * dzm(2) * ( xm(3) - xm(2) )
    !                  - [ u_star^2 / sqrt( um(2)^2 + vm(2)^2 ) ] * xm(2) ];
    !
    !    which can be re-written as:
    !
    !    + dzt(2) * [ K_zm(2) * dzm(2) * ( xm(3) - xm(2) ) ]
    !    - dzt(2) * [ u_star^2 / sqrt( um(2)^2 + vm(2)^2 ) ] * xm(2).
    !
    !    For this equation, a Crank-Nicholson (semi-implicit) diffusion scheme
    !    is used to solve the d [ K_zm * d(xm)/dz ] / dz eddy-diffusion term.
    !    The discretized implicit portion of the term is written out as:
    !
    !    + (1/2)*dzt(2) * [ K_zm(2) * dzm(2) * ( xm(3,<t+1>) - xm(2,<t+1>) ) ]
    !    - dzt(2) * [u_star^2/sqrt( um(2,<t>)^2 + vm(2,<t>)^2 )] * xm(2,<t+1>).
    !
    !    Note:  When the implicit term is brought over to the left-hand side,
    !           the signs are reversed and the leading "+" in front of the first
    !           part of the term is changed to a "-", while the leading "-" in
    !           front of the second part of the term is changed to a "+".
    !
    !    Note:  The x'w'|_sfc portion of the term written above has been pulled
    !           away from the rest of the implicit form written above because
    !           the (1/2) factor due to Crank-Nicholson time_stepping does not
    !           apply to it.  The x'w'|_sfc portion of the term is treated
    !           completely implicitly in order to enhance numerical stability.
    !
    !    The discretized explicit portion of the term is written out as:
    !
    !    + (1/2)*dzt(2) * [ K_zm(2) * dzm(2) * ( xm(3,<t>) - xm(2,<t>) ) ].
    !
    !    Timestep index (t) stands for the index of the current timestep, while
    !    timestep index (t+1) stands for the index of the next timestep, which
    !    is being advanced to in solving the d(xm)/dt equation.
    !
    !
    ! The lower boundary condition for the implicit and explicit portions of the
    ! turbulent advection term, without the x'w'|_sfc portion of the term, can
    ! easily be invoked by using the zero-flux boundary conditions found in the
    ! generalized diffusion function (function diffusion_zt_lhs), which is used
    ! for many other equations in this model.  Either the generalized explicit
    ! surface flux needs to be added onto the explicit term after the diffusion
    ! function has been called from subroutine windm_edsclrm_rhs, or the
    ! implicit momentum surface flux needs to be added onto the implicit term
    ! after the diffusion function has been called from subroutine
    ! windm_edsclrm_lhs.  However, all other equations in this model that use
    ! zero-flux diffusion have level 1 as the level to which the lower boundary
    ! condition needs to be applied.  Thus, an adjuster will have to be used at
    ! level 2 to call diffusion_zt_lhs with level 1 as the input level (the last
    ! variable being passed in during the function call).  However, the other
    ! variables passed in (K_zm, gr%dzt, and gr%dzm variables) will have to be
    ! passed in as solving for level 2.
    !
    ! The value of xm(1) is located below the model surface and does not effect
    ! the rest of the model.  Since xm can be either a horizontal wind component
    ! or a generic eddy scalar quantity, the value of xm(1) is simply set to the
    ! value of xm(2) after the equation matrix has been solved.
    !
    !
    ! Conservation Properties:
    !
    ! When a fixed-flux lower boundary condition is used (combined with a
    ! zero-flux upper boundary condition), this technique of discretizing the
    ! turbulent advection term (treated as an eddy-diffusion term) leads to
    ! conservative differencing.  When the implicit momentum surface flux is
    ! either zero or not used, the column totals for each column in the
    ! left-hand side matrix (for the turbulent advection term) should be equal
    ! to 0.  Otherwise, the column total for the second column will be equal to
    ! the surface flux.   When the generalized explicit surface flux is either
    ! zero or not used, the column total for the right-hand side vector (for the
    ! turbulent advection term) should be equal to 0.  Otherwise, the column
    ! total for the right-hand side vector (for the turbulent advection term)
    ! should be equal to the surface flux.  This ensures that the total amount
    ! of quantity xm over the entire vertical domain is only changed by the
    ! surface flux (neglecting any forcing terms).  The total amount of change
    ! is equal to the surface flux.
    !
    ! To see that this conservation law is satisfied by the left-hand side
    ! matrix, compute the turbulent advection (treated as eddy diffusion) of xm,
    ! neglecting any implicit momentum surface flux, and integrate vertically.
    ! In discretized matrix notation (where "i" stands for the matrix column
    ! and "j" stands for the matrix row):
    !
    !  0 = Sum_j Sum_i ( 1/dzt )_i ( 0.5 * dzt * (K_zm*dzm) )_ij (xm<t+1>)_j.
    !
    ! The left-hand side matrix, ( 0.5 * dzt * (K_zm*dzm) )_ij, is partially
    ! written below.  The sum over i in the above equation removes dzt
    ! everywhere from the matrix below.  The sum over j leaves the column totals
    ! that are desired, which are 0.
    !
    ! Left-hand side matrix contributions from the turbulent advection term
    ! (treated as an eddy-diffusion term using a Crank-Nicholson timestep);
    ! first five vertical levels:
    !
    !     ------------------------------------------------------------------------------->
    !k=1 |  0             0                        0                          0
    !    |
    !k=2 |  0   +0.5*dzt(k)*           -0.5*dzt(k)*                           0
    !    |        K_zm(k)*dzm(k)         K_zm(k)*dzm(k)
    !    |
    !k=3 |  0   -0.5*dzt(k)*           +0.5*dzt(k)*               -0.5*dzt(k)*
    !    |        K_zm(k-1)*dzm(k-1)     [ K_zm(k)*dzm(k)           K_zm(k)*dzm(k)
    !    |                                +K_zm(k-1)*dzm(k-1) ]
    !    |
    !k=4 |  0             0            -0.5*dzt(k)*               +0.5*dzt(k)*
    !    |                               K_zm(k-1)*dzm(k-1)         [ K_zm(k)*dzm(k)
    !    |                                                           +K_zm(k-1)*dzm(k-1) ]
    !    |
    !k=5 |  0             0                        0              -0.5*dzt(k)*
    !    |                                                          K_zm(k-1)*dzm(k-1)
    !   \ /
    !
    ! Note:  The superdiagonal term from level 4 and both the main diagonal and
    !        superdiagonal terms from level 5 are not shown on this diagram.
    !
    ! Note:  If an implicit momentum surface flux is used, an additional term,
    !        + dzt(2) * [ u_star^2 / sqrt( um(2,<t>)^2 + vm(2,<t>)^2 ) ], is
    !        added to row 2 (k=2), column 2.
    !
    ! To see that the above conservation law is satisfied by the right-hand side
    ! vector, compute the turbulent advection (treated as eddy diffusion) of xm,
    ! neglecting any generalized explicit surface flux, and integrate
    ! vertically.  In discretized matrix notation (where "i" stands for the
    ! matrix column and "j" stands for the matrix row):
    !
    !  0 = Sum_j Sum_i ( 1/dzt )_i ( rhs_vector )_j.
    !
    ! The right-hand side vector, ( rhs_vector )_j, is partially written below.
    ! The sum over i in the above equation removes dzt everywhere from the
    ! vector below.  The sum over j leaves the column total that is desired,
    ! which is 0.
    !
    ! Right-hand side vector contributions from the turbulent advection term
    ! (treated as an eddy-diffusion term using a Crank-Nicholson timestep);
    ! first five vertical levels:
    !
    !     --------------------------------------------
    !k=1 |                      0                     |
    !    |                                            |
    !    |                                            |
    !k=2 | +0.5*dzt(k)*                               |
    !    |       [ K_zm(k)*dzm(k)*                    |
    !    |                  (xm(k+1,<t>)-xm(k,<t>)) ] |
    !    |                                            |
    !k=3 | +0.5*dzt(k)*                               |
    !    |       [ K_zm(k)*dzm(k)*                    |
    !    |                  (xm(k+1,<t>)-xm(k,<t>))   |
    !    |        -K_zm(k-1)*dzm(k-1)*                |
    !    |                  (xm(k,<t>)-xm(k-1,<t>)) ] |
    !    |                                            |
    !k=4 | +0.5*dzt(k)*                               |
    !    |       [ K_zm(k)*dzm(k)*                    |
    !    |                  (xm(k+1,<t>)-xm(k,<t>))   |
    !    |        -K_zm(k-1)*dzm(k-1)*                |
    !    |                  (xm(k,<t>)-xm(k-1,<t>)) ] |
    !    |                                            |
    !k=5 | +0.5*dzt(k)*                               |
    !    |       [ K_zm(k)*dzm(k)*                    |
    !    |                  (xm(k+1,<t>)-xm(k,<t>))   |
    !    |        -K_zm(k-1)*dzm(k-1)*                |
    !    |                  (xm(k,<t>)-xm(k-1,<t>)) ] |
    !   \ /                                          \ /
    !
    ! Note:  If a generalized explicit surface flux is used, an additional term,
    !        + dzt(2) * x'w'|_sfc, is added to row 2 (k=2).
    !
    ! Note:  Only the contributions by the turbulent advection term are shown
    !        for both the left-hand side matrix and the right-hand side vector.
    !        There are more terms in the equation, and thus more factors to be
    !        added to both the left-hand side matrix (such as time tendency and
    !        mean advection) and the right-hand side vector (such as xm
    !        forcings).  The left-hand side matrix is set-up so that a singular
    !        matrix is not encountered.

    ! References:
    ! Eqn. 8 & 9 on p. 3545 of
    ! ``A PDF-Based Model for Boundary Layer Clouds. Part I:
    !   Method and Model Description'' Golaz, et al. (2002)
    ! JAS, Vol. 59, pp. 3540--3551.
    !-----------------------------------------------------------------------

    use grid_class, only: & 
      gr ! Variable(s)

    use lapack_wrap, only:  & 
      tridag_solve, & ! Procedure(s)
      tridag_solvex

    use stats_variables, only: & 
      sfc,     &   ! Variable(s)
      l_stats_samp

    use stats_type, only:  &
      stat_update_var_pt  ! Subroutine

    use constants, only:  & 
      fstderr ! Variable(s)

    implicit none

    ! Constant parameters

    integer, parameter :: &
      kp1_tdiag = 1,    & ! Thermodynamic superdiagonal index.
      k_tdiag   = 2,    & ! Thermodynamic main diagonal index.
      km1_tdiag = 3       ! Thermodynamic subdiagonal index.

    ! Input Variables

    integer, intent(in) :: &
      nrhs   ! Number of right-hand side (explicit) vectors.
    ! Number of solution vectors.

    integer, intent(in) :: &
      ixm_matrix_condt_num  ! Stats index of the condition numbers

    real, dimension(3,gr%nnzp), intent(inout) :: &
      lhs    ! Implicit contributions to um, vm, and eddy scalars  [units vary]

    real, dimension(gr%nnzp,nrhs), intent(inout) :: &
      rhs    ! Right-hand side (explicit) contributions.

    real, dimension(gr%nnzp,nrhs), intent(out) :: &
      solution ! Solution to the system of equations    [units vary]

    integer, intent(out) :: & 
      err_code ! clubb_singular_matrix when matrix is singular

    ! Local variables
    real :: rcond ! Estimate of the reciprocal of the condition number on the LHS matrix

    ! Solve tridiagonal system for xm.
    if ( l_stats_samp .and. ixm_matrix_condt_num > 0 ) then
      call tridag_solvex & 
           ( "windm_edsclrm", gr%nnzp, nrhs, &                          ! Intent(in) 
             lhs(kp1_tdiag,:), lhs(k_tdiag,:), lhs(km1_tdiag,:), rhs, & ! Intent(inout)
             solution, rcond, err_code )                                ! Intent(out)

      ! Est. of the condition number of the variance LHS matrix
      call stat_update_var_pt( ixm_matrix_condt_num, 1, 1.0 / rcond, &  ! Intent(in)
                               sfc )                                       ! Intent(inout)
    else

      call tridag_solve( "windm_edsclrm", gr%nnzp, nrhs, &                           ! In
                         lhs(kp1_tdiag,:),  lhs(k_tdiag,:), lhs(km1_tdiag,:), rhs, & ! Inout
                         solution, err_code )                                        ! Out
    end if

    return
  end subroutine windm_edsclrm_solve

  !=============================================================================
  subroutine windm_edsclrm_implicit_stats( solve_type, xm )

    ! Description:
    ! Compute implicit contributions to um and vm

    ! References:
    ! None
    !-----------------------------------------------------------------------

    use stats_variables, only: & 
      ium_ma,  & ! Variables
      ium_ta,  & 
      ivm_ma,  &
      ivm_ta,  & 
      ztscr01, & 
      ztscr02, & 
      ztscr03, & 
      ztscr04, & 
      ztscr05, & 
      ztscr06, & 
      zt

    use stats_type, only:  &
      stat_end_update_pt,  & ! Subroutines
      stat_update_var_pt

    use constants, only:  & 
      fstderr ! Variable(s)

    use stats_precision, only:  & 
      time_precision ! Variable(s)

    use grid_class, only: &
      gr ! Derived type variable

    implicit none

    ! Input variables
    character(len=*), intent(in) :: & 
      solve_type     ! Desc. of what is being solved for

    real, dimension(gr%nnzp), intent(in) :: &
      xm !  Computed value um or vm at <t+1>    [m/s]

    ! Local variables
    integer :: k, kp1, km1 ! Array indices

    ! Budget indices
    integer :: ixm_ma, ixm_ta

    select case ( trim( solve_type ) )
    case ( "um" )
      ixm_ma = ium_ma
      ixm_ta = ium_ta

    case ( "vm" )
      ixm_ma = ivm_ma
      ixm_ta = ivm_ta

    case default
      ixm_ma = 0
      ixm_ta = 0

    end select


    ! Finalize implicit contributions for xm

    do k = 2, gr%nnzp-1, 1

      km1 = max( k-1, 1 )
      kp1 = min( k+1, gr%nnzp )

      ! xm mean advection
      ! xm term ma is completely implicit; call stat_update_var_pt.
      call stat_update_var_pt( ixm_ma, k, &
             ztscr01(k) * xm(km1) &
           + ztscr02(k) * xm(k) &
           + ztscr03(k) * xm(kp1), zt )

      ! xm turbulent transport (implicit component)
      ! xm term ta has both implicit and explicit components;
      ! call stat_end_update_pt.
      call stat_end_update_pt( ixm_ta, k, &
             ztscr04(k) * xm(km1) &
           + ztscr05(k) * xm(k) &
           + ztscr06(k) * xm(kp1), zt )

    enddo

    ! Upper boundary condition:
    ! xm turbulent transport (implicit component)
    ! xm term ta has both implicit and explicit components;
    ! call stat_end_update_pt.
    k   = gr%nnzp
    km1 = max( k-1, 1 )
    call stat_end_update_pt( ixm_ta, k, &
           ztscr04(k) * xm(km1) &
         + ztscr05(k) * xm(k), zt )

    return
  end subroutine windm_edsclrm_implicit_stats

  !=============================================================================
  subroutine compute_uv_tndcy( solve_type, fcor, perp_wind_m, perp_wind_g, xm_forcing, &
                               l_implemented, xm_tndcy )

    ! Description:
    ! Computes the explicit tendency for the um and vm wind components.
    !
    ! The only explicit tendency that is involved in the d(um)/dt or d(vm)/dt
    ! equations is the Coriolis tendency.
    !
    ! The d(um)/dt equation contains the term:
    !
    ! - f * ( v_g - vm );
    !
    ! where f is the Coriolis parameter and v_g is the v component of the
    ! geostrophic wind.
    !
    ! Likewise, the d(vm)/dt equation contains the term:
    !
    ! + f * ( u_g - um );
    !
    ! where u_g is the u component of the geostrophic wind.
    !
    ! This term is treated completely explicitly.  The values of um, vm, u_g,
    ! and v_g are all found on the thermodynamic levels.
    !
    ! Wind forcing from the GCSS cases is also added here.
    !
    ! References:
    !-----------------------------------------------------------------------

    use grid_class, only: & 
        gr

    use stats_type, only: & 
        stat_update_var

    use stats_variables, only:      &
        ium_gf, & 
        ium_cf, & 
        ivm_gf, & 
        ivm_cf, & 
        ium_f,  &
        ivm_f,  &
        zt, & 
        l_stats_samp

    implicit none

    ! Input Variables
    character(len=*), intent(in) ::  &
      solve_type      ! Description of what is being solved for

    real, intent(in) ::  & 
      fcor            ! Coriolis parameter     [s^-1]

    real, dimension(gr%nnzp), intent(in) :: & 
      perp_wind_m,  & ! Perpendicular component of the mean wind (e.g. v, for the u-eqn) [m/s]
      perp_wind_g,  & ! Perpendicular component of the geostropic wind (e.g. vg)         [m/s]
      xm_forcing      ! Prescribed wind forcing                                          [m/s/s]

    logical, intent(in) :: & 
      l_implemented   ! Flag for CLUBB being implemented in a larger model.

    ! Output Variables
    real, dimension(gr%nnzp), intent(out) ::  &
      xm_tndcy        ! xm tendency            [m/s^2]

    ! Local Variables
    integer :: & 
      ixm_gf, & 
      ixm_cf, &
      ixm_f

    real, dimension(gr%nnzp) :: & 
      xm_gf, & 
      xm_cf

    ! --- Begin Code ---

    if ( .not. l_implemented ) then
      ! Only compute the Coriolis term if the model is running on it's own,
      ! and is not part of a larger, host model.

      select case ( trim( solve_type ) )

      case ( "um" )

        ixm_gf = ium_gf
        ixm_cf = ium_cf
        ixm_f  = ium_f

        xm_gf = - fcor * perp_wind_g(1:gr%nnzp)

        xm_cf = fcor * perp_wind_m(1:gr%nnzp)

      case ( "vm" )

        ixm_gf = ivm_gf
        ixm_cf = ivm_cf
        ixm_f  = ivm_f

        xm_gf = fcor * perp_wind_g(1:gr%nnzp)

        xm_cf = -fcor * perp_wind_m(1:gr%nnzp)

      case default

        ixm_gf = 0
        ixm_cf = 0
        ixm_f = 0

        xm_gf = 0.

        xm_cf = 0.

      end select

      xm_tndcy(1:gr%nnzp) = xm_gf(1:gr%nnzp) + xm_cf(1:gr%nnzp) + xm_forcing(1:gr%nnzp)

      if ( l_stats_samp ) then

        ! xm term gf is completely explicit; call stat_update_var.
        call stat_update_var( ixm_gf, xm_gf, zt )

        ! xm term cf is completely explicit; call stat_update_var.
        call stat_update_var( ixm_cf, xm_cf, zt )

        ! xm term F
        call stat_update_var( ixm_f, xm_forcing, zt )
      endif

    else   ! implemented in a host model.

      xm_tndcy = 0.0

    endif


    return
  end subroutine compute_uv_tndcy

!===============================================================================
  subroutine windm_edsclrm_lhs( dt, wm_zt, Kh_zm, wind_speed, u_star_sqd,   &
                                l_implemented, l_imp_sfc_momentum_flux,  &
                                lhs )

    ! Description:
    ! Calculate the implicit portion of the horizontal wind or eddy-scalar
    ! time-tendency equation.  See the description in subroutine
    ! windm_edsclrm_solve for more details.

    ! References:
    ! None
    !-----------------------------------------------------------------------

    use grid_class, only:  & 
        gr  ! Variable(s)

    use stats_precision, only:  & 
        time_precision ! Variable(s)

    use diffusion, only:  & 
        diffusion_zt_lhs ! Procedure(s)

    use mean_adv, only: & 
        term_ma_zt_lhs  ! Procedures

    use stats_variables, only: &
        ium_ma,  & ! Variable(s)
        ium_ta,  &
        ivm_ma,  &
        ivm_ta,  &
        ztscr01, &
        ztscr02, &
        ztscr03, &
        ztscr04, &
        ztscr05, &
        ztscr06, &
        l_stats_samp

    implicit none

    ! Constant parameters
    integer, parameter :: & 
      kp1_tdiag = 1,    & ! Thermodynamic superdiagonal index.
      k_tdiag   = 2,    & ! Thermodynamic main diagonal index.
      km1_tdiag = 3       ! Thermodynamic subdiagonal index.

    ! Input Variables
    real(kind=time_precision), intent(in) :: & 
      dt            ! Model timestep                           [s]

    real, dimension(gr%nnzp), intent(in) :: &
      wm_zt,      & ! w wind component on thermodynamic levels [m/s]
      Kh_zm,      & ! Eddy diffusivity on momentum levels      [m^2/s]
      wind_speed    ! wind speed; sqrt( u^2 + v^2 )            [m/s]

    real, intent(in) :: &
      u_star_sqd    ! Surface friction velocity, u_*, squared  [m/s]

    logical, intent(in) ::  & 
      l_implemented, & ! Flag for CLUBB being implemented in a larger model.
      l_imp_sfc_momentum_flux  ! Flag for implicit momentum surface fluxes.

    ! Output Variable
    real, dimension(3,gr%nnzp), intent(out) :: &
      lhs           ! Implicit contributions to xm (tridiagonal matrix)

    ! Local Variables
    integer :: k, km1  ! Array indices
    integer :: diff_k_in

    real, dimension(3) :: tmp

    ! --- Begin Code ---

    ! Initialize the LHS array.
    lhs = 0.0

    do k = 2, gr%nnzp, 1

      ! Define index
      km1 = max( k-1, 1 )

      ! LHS mean advection term.
      if ( .not. l_implemented ) then

        lhs(kp1_tdiag:km1_tdiag,k)  &
        = lhs(kp1_tdiag:km1_tdiag,k)  &
        + term_ma_zt_lhs( wm_zt(k), gr%dzt(k), k )

      else

        lhs(kp1_tdiag:km1_tdiag,k)  &
        = lhs(kp1_tdiag:km1_tdiag,k) + 0.0

      endif

      ! LHS turbulent advection term (solved as an eddy-diffusion term).
      if ( k == 2 ) then
        ! The lower boundary condition needs to be applied here at level 2.
        ! The lower boundary condition is a "fixed flux" boundary condition.
        ! The coding is the same as for a zero-flux boundary condition, but with
        ! an extra term added on the right-hand side at the boundary level.  For
        ! the rest of the model code, a zero-flux boundary condition is applied
        ! at level 1, and thus subroutine diffusion_zt_lhs is set-up to do that.
        ! In order to apply the same boundary condition code here at level 2, an
        ! adjuster needs to be used to tell diffusion_zt_lhs to use the code at
        ! level 2 that it normally uses at level 1.
        diff_k_in = 1
      else
        diff_k_in = k
      endif
      lhs(kp1_tdiag:km1_tdiag,k)  &
      = lhs(kp1_tdiag:km1_tdiag,k)  &
      + 0.5  &
      * diffusion_zt_lhs( Kh_zm(k), Kh_zm(km1), 0.0,  & 
                          gr%dzm(km1), gr%dzm(k), gr%dzt(k), diff_k_in )

      ! LHS time tendency.
      lhs(k_tdiag,k)  &
      = real( lhs(k_tdiag,k) + ( 1.0 / dt ) )

      if ( l_stats_samp ) then

        ! Statistics:  implicit contributions for um or vm.

        if ( ium_ma + ivm_ma > 0 ) then
          if ( .not. l_implemented ) then
            tmp(1:3) &
            = term_ma_zt_lhs( wm_zt(k), gr%dzt(k), k )
            ztscr01(k) = -tmp(3)
            ztscr02(k) = -tmp(2)
            ztscr03(k) = -tmp(1)
          else
            ztscr01(k) = 0.0
            ztscr02(k) = 0.0
            ztscr03(k) = 0.0
          endif
        endif

        if ( ium_ta + ivm_ta > 0 ) then
          tmp(1:3)  &
          = 0.5  &
          * diffusion_zt_lhs( Kh_zm(k), Kh_zm(km1), 0.0,  &
                              gr%dzm(km1), gr%dzm(k), gr%dzt(k), diff_k_in )
          ztscr04(k) = -tmp(3)
          ztscr05(k) = -tmp(2)
          ztscr06(k) = -tmp(1)
        endif

      endif  ! l_stats_samp

    enddo ! k = 2 .. gr%nnzp


    ! Boundary Conditions

    ! Lower Boundary

    ! The lower boundary condition is a fixed-flux boundary condition, which
    ! gets added into the time-tendency equation at level 2.
    ! The value of xm(1) is located below the model surface and does not effect
    ! the rest of the model.  Since xm can be either a horizontal wind component
    ! or a generic eddy scalar quantity, the value of xm(1) is simply set to the
    ! value of xm(2) after the equation matrix has been solved.

    ! k = 1
    lhs(k_tdiag,1) = 1.0

    ! k = 2; add implicit momentum surface flux.
    if ( l_imp_sfc_momentum_flux ) then

      ! LHS momentum surface flux.
      lhs(k_tdiag,2)  &
      = lhs(k_tdiag,2)  &
      + gr%dzt(2) * ( u_star_sqd / wind_speed(2) )

      if ( l_stats_samp ) then

        ! Statistics:  implicit contributions for um or vm.

        ! xm term ta is modified at level 2 to include the effects of the
        ! surface flux.  In this case, this effects the implicit portion of
        ! the term (after zmscr05, which handles the main diagonal for the
        ! turbulent advection term, has already been called at level 2).
        ! Modify zmscr05 accordingly.
        if ( ium_ta + ivm_ta > 0 ) then
          ztscr05(2)  &
          = ztscr05(2)  &
          - gr%dzt(2) * ( u_star_sqd / wind_speed(2) )
        endif

      endif ! l_stats_samp

    endif ! l_imp_sfc_momentum_flux


    return
  end subroutine windm_edsclrm_lhs

  !=============================================================================
  function windm_edsclrm_rhs( solve_type, dt, Kh_zm, xm, xm_tndcy,  &
                              l_imp_sfc_momentum_flux, xpwp_sfc )  &
  result( rhs )

    ! Description:
    ! Calculate the explicit portion of the horizontal wind or eddy-scalar
    ! time-tendency equation.  See the description in subroutine
    ! windm_edsclrm_solve for more details.

    ! References:
    ! None
    !-----------------------------------------------------------------------

    use stats_precision, only:  & 
        time_precision ! Variable(s)

    use diffusion, only:  & 
        diffusion_zt_lhs ! Procedure(s)

    use stats_variables, only: &
        ium_ta,  & ! Variable(s)
        ivm_ta,  &
        zt,      &
        l_stats_samp

    use stats_type, only: &
        stat_begin_update_pt,  & ! Procedure(s)
        stat_modify_pt

    use grid_class, only:  & 
        gr  ! Variable(s)

    implicit none

    ! External
    intrinsic :: max, min, real, trim

    ! Input Variables
    character(len=*), intent(in) :: &
      solve_type ! Description of what is being solved for

    real(kind=time_precision), intent(in) :: & 
      dt           ! Model timestep                                   [s]

    real, dimension(gr%nnzp), intent(in) :: &
      Kh_zm,     & ! Eddy diffusivity on momentum levels              [m^2/s]
      xm,        & ! Eddy-scalar variable, xm (thermodynamic levels)  [units vary]
      xm_tndcy     ! The explicit time-tendency acting on xm          [units vary]

    real, intent(in) :: &
      xpwp_sfc     ! x'w' at the surface                              [units vary]

    logical, intent(in) :: &
      l_imp_sfc_momentum_flux  ! Flag for implicit momentum surface fluxes.

    ! Output Variable
    real, dimension(gr%nnzp) :: &
      rhs          ! Right-hand side (explicit) contributions.

    ! Local Variables
    integer :: k, kp1, km1  ! Array indices
    integer :: diff_k_in

    ! For use in Crank-Nicholson eddy diffusion.
    real, dimension(3) :: rhs_diff

    integer :: ixm_ta

    ! --- Begin Code ---

    select case ( trim( solve_type ) )
    case ( "um" )
      ixm_ta = ium_ta
    case ( "vm" )
      ixm_ta = ivm_ta
    case default  ! Eddy scalars
      ixm_ta = 0
    end select


    ! Initialize the RHS vector.
    rhs = 0.0

    do k = 2, gr%nnzp-1, 1

      ! Define indices
      km1 = max( k-1, 1 )
      kp1 = min( k+1, gr%nnzp )

      ! RHS turbulent advection term (solved as an eddy-diffusion term).
      if ( k == 2 ) then
        ! The lower boundary condition needs to be applied here at level 2.
        ! The lower boundary condition is a "fixed flux" boundary condition.
        ! The coding is the same as for a zero-flux boundary condition, but with
        ! an extra term added on the right-hand side at the boundary level.  For
        ! the rest of the model code, a zero-flux boundary condition is applied
        ! at level 1, and thus subroutine diffusion_zt_lhs is set-up to do that.
        ! In order to apply the same boundary condition code here at level 2, an
        ! adjuster needs to be used to tell diffusion_zt_lhs to use the code at
        ! level 2 that it normally uses at level 1.
        diff_k_in = 1
      else
        diff_k_in = k
      endif
      rhs_diff(1:3)  & 
      = 0.5  &
      * diffusion_zt_lhs( Kh_zm(k), Kh_zm(km1), 0.0,  &
                          gr%dzm(km1), gr%dzm(k), gr%dzt(k), diff_k_in )
      rhs(k)   =   rhs(k) & 
                 - rhs_diff(3) * xm(km1) &
                 - rhs_diff(2) * xm(k)   &
                 - rhs_diff(1) * xm(kp1)

      ! RHS forcings.
      rhs(k) = rhs(k) + xm_tndcy(k)

      ! RHS time tendency
      rhs(k) = real( rhs(k) + ( 1.0 / dt ) * xm(k) )

      if ( l_stats_samp ) then

        ! Statistics:  explicit contributions for um or vm.

        ! xm term ta has both implicit and explicit components; call
        ! stat_begin_update_pt.  Since stat_begin_update_pt automatically
        ! subtracts the value sent in, reverse the sign on right-hand side
        ! turbulent advection component.
        if ( ixm_ta > 0 ) then
          call stat_begin_update_pt( ixm_ta, k, & 
                 rhs_diff(3) * xm(km1) &
               + rhs_diff(2) * xm(k)   &
               + rhs_diff(1) * xm(kp1), zt )
        endif

      endif  ! l_stats_samp

    enddo


    ! Boundary Conditions

    ! Lower Boundary

    ! The lower boundary condition is a fixed-flux boundary condition, which
    ! gets added into the time-tendency equation at level 2.
    ! The value of xm(1) is located below the model surface and does not effect
    ! the rest of the model.  Since xm can be either a horizontal wind component
    ! or a generic eddy scalar quantity, the value of xm(1) is simply set to the
    ! value of xm(2) after the equation matrix has been solved.  For purposes of
    ! the matrix equation, rhs(1) is simply set to 0.

    ! k = 1
    rhs(1) = 0.0

    ! k = 2; add generalized explicit surface flux.
    if ( .not. l_imp_sfc_momentum_flux ) then

      ! RHS generalized surface flux.
      rhs(2) = rhs(2) + gr%dzt(2) * xpwp_sfc

      if ( l_stats_samp ) then

        ! Statistics:  explicit contributions for um or vm.

        ! xm term ta is modified at level 2 to include the effects of the
        ! surface flux.  In this case, this effects the explicit portion of
        ! the term (after stat_begin_update_pt has already been called at
        ! level 2); call stat_modify_pt.
        if ( ixm_ta > 0 ) then
          call stat_modify_pt( ixm_ta, 2, &
               + gr%dzt(2) * xpwp_sfc, zt )
        endif

      endif  ! l_stats_samp

    endif ! l_imp_sfc_momentum_flux

    ! Upper Boundary

    ! A zero-flux boundary condition is used at the upper boundary, meaning that
    ! xm is not allowed to exit the model through the upper boundary.  This
    ! boundary condition is invoked by calling diffusion_zt_lhs at the uppermost
    ! level.
    k   = gr%nnzp
    km1 = max( k-1, 1 )

    ! RHS turbulent advection term (solved as an eddy-diffusion term) at the
    ! upper boundary.
    rhs_diff(1:3)  &
    = 0.5  &
    * diffusion_zt_lhs( Kh_zm(k), Kh_zm(km1), 0.0,  &
                        gr%dzm(km1), gr%dzm(k), gr%dzt(k), k )
    rhs(k)   =   rhs(k) &
               - rhs_diff(3) * xm(km1) &
               - rhs_diff(2) * xm(k)

    ! RHS forcing term at the upper boundary.
    rhs(k) = rhs(k) + xm_tndcy(k)

    ! RHS time tendency term at the upper boundary.
    rhs(k) = real( rhs(k) + ( 1.0 / dt ) * xm(k) )

    if ( l_stats_samp ) then

      ! Statistics:  explicit contributions for um or vm.

      ! xm term ta has both implicit and explicit components; call
      ! stat_begin_update_pt.  Since stat_begin_update_pt automatically
      ! subtracts the value sent in, reverse the sign on right-hand side
      ! turbulent advection component.
      if ( ixm_ta > 0 ) then
        call stat_begin_update_pt( ixm_ta, k, &
               rhs_diff(3) * xm(km1) &
             + rhs_diff(2) * xm(k), zt )
      endif

    endif  ! l_stats_samp


    return
  end function windm_edsclrm_rhs

!===============================================================================
  elemental function xpwp_fnc( Kh_zm, xm, xmp1, dzm )

    ! Description:
    ! Compute x'w' from x<k>, x<k+1>, Kh and dzm

    ! References:
    ! None
    !-----------------------------------------------------------------------

    implicit none

    ! Input variables
    real, intent(in) :: &
      Kh_zm, & ! Eddy diff. (k momentum level)                 [m^2/s]
      xm,    & ! x (k thermo level)                            [units vary]
      xmp1,  & ! x (k+1 thermo level)                          [units vary]
      dzm      ! Inverse of the grid spacing (k thermo level)  [1/m]
    ! Output variable
    real :: &
      xpwp_fnc ! x'w'   [(units vary)(m/s)]

    !-----------------------------------------------------------------------
    ! --- Begin Code ---

    ! Solve for x'w' at all intermediate model levels.
    xpwp_fnc = Kh_zm * dzm * ( xmp1 - xm )

    return
  end function xpwp_fnc


!===============================================================================

end module advance_windm_edsclrm_module
