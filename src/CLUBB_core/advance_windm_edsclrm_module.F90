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
               edsclrm_forcing, &
               rho_ds_zm, invrs_rho_ds_zt, &
               fcor, l_implemented, &
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

    use constants_clubb, only:  & 
        fstderr, &  ! Constant
        eps

    use sponge_layer_damping, only: &
      uv_sponge_damp_settings, &  
      uv_sponge_damp_profile, &  
      sponge_damp_xm ! Procedure(s)

    implicit none

    ! Input Variables
    real(kind=time_precision), intent(in) ::  &
      dt                 ! Model timestep                             [s]

    real, dimension(gr%nnzp), intent(in) ::  &
      wm_zt,           & ! w wind component on thermodynamic levels   [m/s]
      Kh_zm,           & ! Eddy diffusivity on momentum levels        [m^2/s]
      ug,              & ! u (west-to-east) geostrophic wind comp.    [m/s]
      vg,              & ! v (south-to-north) geostrophic wind comp.  [m/s]
      um_ref,          & ! Reference u wind component for nudging     [m/s]
      vm_ref,          & ! Reference v wind component for nudging     [m/s]
      wp2,             & ! w'^2 (momentum levels)                     [m^2/s^2]
      up2,             & ! u'^2 (momentum levels)                     [m^2/s^2]
      vp2,             & ! v'^2 (momentum levels)                     [m^2/s^2]
      um_forcing,      & ! u forcing                                  [m/s/s]
      vm_forcing,      & ! v forcing                                  [m/s/s]
      rho_ds_zm,       & ! Dry, static density on momentum levels     [kg/m^3]
      invrs_rho_ds_zt    ! Inv. dry, static density at thermo. levels [m^3/kg]

    real, dimension(gr%nnzp,edsclr_dim), intent(in) ::  &
      edsclrm_forcing  ! Eddy scalar large-scale forcing             [{units vary}/s]

    real, intent(in) ::  &
      fcor           ! Coriolis parameter                            [s^-1]

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
    real, dimension(gr%nnzp), intent(inout) ::  &
      upwp,        & ! u'w' (momentum levels)                        [m^2/s^2]
      vpwp           ! v'w' (momentum levels)                        [m^2/s^2]

    ! Output Variable for eddy-scalars
    real, dimension(gr%nnzp,edsclr_dim), intent(inout) ::  &
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
    u_star_sqd = sqrt( upwp(1)**2 + vpwp(1)**2 )

    ! Compute the explicit portion of the um equation.
    ! Build the right-hand side vector.
    rhs(1:gr%nnzp,1) = windm_edsclrm_rhs( "um", dt, Kh_zm, um, um_tndcy,  &  ! in
                                          rho_ds_zm, invrs_rho_ds_zt,  &     ! in
                                          l_imp_sfc_momentum_flux, upwp(1) ) ! in

    ! Compute the explicit portion of the vm equation.
    ! Build the right-hand side vector.
    rhs(1:gr%nnzp,2) = windm_edsclrm_rhs( "vm", dt, Kh_zm, vm, vm_tndcy,  &  ! in
                                          rho_ds_zm, invrs_rho_ds_zt,  &     ! in
                                          l_imp_sfc_momentum_flux, vpwp(1) ) ! in


    ! Store momentum flux (explicit component)

    ! The surface flux, x'w'(1) = x'w'|_sfc, is set elsewhere in the model.
!   upwp(1) = upwp_sfc
!   vpwp(1) = vpwp_sfc

    ! Solve for x'w' at all intermediate model levels.
    ! A Crank-Nicholson timestep is used.

    upwp(2:gr%nnzp-1) = - 0.5 * xpwp_fnc( Kh_zm(2:gr%nnzp-1), um(2:gr%nnzp-1), & ! in
                                          um(3:gr%nnzp), gr%invrs_dzm(2:gr%nnzp-1) )   ! in

    vpwp(2:gr%nnzp-1) = - 0.5 * xpwp_fnc( Kh_zm(2:gr%nnzp-1), vm(2:gr%nnzp-1), & ! in
                                          vm(3:gr%nnzp), gr%invrs_dzm(2:gr%nnzp-1) )   ! in

    ! A zero-flux boundary condition at the top of the model, d(xm)/dz = 0,
    ! means that x'w' at the top model level is 0,
    ! since x'w' = - K_zm * d(xm)/dz.
    upwp(gr%nnzp) = 0.
    vpwp(gr%nnzp) = 0.


    ! Compute the implicit portion of the um and vm equations.
    ! Build the left-hand side matrix.
    call windm_edsclrm_lhs( dt, wm_zt, Kh_zm, wind_speed, u_star_sqd,  & ! in
                            rho_ds_zm, invrs_rho_ds_zt,  &               ! in
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

    upwp(2:gr%nnzp-1) = upwp(2:gr%nnzp-1) - 0.5 * xpwp_fnc( Kh_zm(2:gr%nnzp-1), &
      um(2:gr%nnzp-1), um(3:gr%nnzp), gr%invrs_dzm(2:gr%nnzp-1) )!in

    vpwp(2:gr%nnzp-1) = vpwp(2:gr%nnzp-1) - 0.5 * xpwp_fnc( Kh_zm(2:gr%nnzp-1), &
      vm(2:gr%nnzp-1), vm(3:gr%nnzp), gr%invrs_dzm(2:gr%nnzp-1) )!in


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
                             rho_ds_zm, invrs_rho_ds_zt,  &                ! in
                             l_imp_sfc_momentum_flux, wpedsclrp(1,i) )     ! in
      enddo


      ! Store momentum flux (explicit component)

      ! The surface flux, x'w'(1) = x'w'|_sfc, is set elsewhere in the model.
!     wpedsclrp(1,1:edsclr_dim) =  wpedsclrp_sfc(1:edsclr_dim)

      ! Solve for x'w' at all intermediate model levels.
      ! A Crank-Nicholson timestep is used.
      ! Here we use a forall and high performance fortran directive to try to
      ! parallelize this computation.  Note that FORALL is more restrictive than DO.
!HPF$ INDEPENDENT, REDUCTION(wpedsclrp)
      forall( i = 1:edsclr_dim )
        wpedsclrp(2:gr%nnzp-1,i) = &
          - 0.5 * xpwp_fnc( Kh_zm(2:gr%nnzp-1), edsclrm(2:gr%nnzp-1,i), & ! in
                            edsclrm(3:gr%nnzp,i), gr%invrs_dzm(2:gr%nnzp-1) )   ! in
      end forall

      ! A zero-flux boundary condition at the top of the model, d(xm)/dz = 0,
      ! means that x'w' at the top model level is 0,
      ! since x'w' = - K_zm * d(xm)/dz.
      wpedsclrp(gr%nnzp,1:edsclr_dim) = 0.


      ! Compute the implicit portion of the xm (eddy-scalar) equations.
      ! Build the left-hand side matrix.
      call windm_edsclrm_lhs( dt, wm_zt, Kh_zm, wind_speed, u_star_sqd,  & ! in
                              rho_ds_zm, invrs_rho_ds_zt,  &               ! in
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
                            edsclrm(3:gr%nnzp,i), gr%invrs_dzm(2:gr%nnzp-1) )   ! in
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
      write(fstderr,*) "fcor = ", fcor
      write(fstderr,*) "l_implemented = ", l_implemented

      write(fstderr,*) "Intent(inout)"

      write(fstderr,*) "um = ", um
      write(fstderr,*) "vm = ", vm
      write(fstderr,*) "edsclrm = ", edsclrm

      write(fstderr,*) "upwp = ", upwp
      write(fstderr,*) "vpwp = ", vpwp
      write(fstderr,*) "wpedsclrp = ", wpedsclrp

      write(fstderr,*) "Intent(out)"


      return

    endif

    return
  end subroutine advance_windm_edsclrm

  !=============================================================================
  subroutine windm_edsclrm_solve( nrhs, ixm_matrix_condt_num, &
                                  lhs, rhs, solution, err_code )


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

    use constants_clubb, only:  & 
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

    use constants_clubb, only:  & 
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


    ! Upper boundary conditions
    k   = gr%nnzp
    km1 = max( k-1, 1 )

    ! xm mean advection
    ! xm term ma is completely implicit; call stat_update_var_pt.
    call stat_update_var_pt( ixm_ma, k, &
           ztscr01(k) * xm(km1) &
         + ztscr02(k) * xm(k), zt )

    ! xm turbulent transport (implicit component)
    ! xm term ta has both implicit and explicit components;
    ! call stat_end_update_pt.
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
  subroutine windm_edsclrm_lhs( dt, wm_zt, Kh_zm, wind_speed, u_star_sqd,  &
                                rho_ds_zm, invrs_rho_ds_zt,  &
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
      dt                 ! Model timestep                             [s]

    real, dimension(gr%nnzp), intent(in) :: &
      wm_zt,           & ! w wind component on thermodynamic levels   [m/s]
      Kh_zm,           & ! Eddy diffusivity on momentum levels        [m^2/s]
      wind_speed,      & ! wind speed; sqrt( u^2 + v^2 )              [m/s]
      rho_ds_zm,       & ! Dry, static density on momentum levels     [kg/m^3]
      invrs_rho_ds_zt    ! Inv. dry, static density at thermo. levels [m^3/kg]

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
        + term_ma_zt_lhs( wm_zt(k), gr%invrs_dzt(k), k )

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
      + 0.5 * invrs_rho_ds_zt(k)  &
      * diffusion_zt_lhs( rho_ds_zm(k) * Kh_zm(k),  &
                          rho_ds_zm(km1) * Kh_zm(km1), 0.0,  &
                          gr%invrs_dzm(km1), gr%invrs_dzm(k), gr%invrs_dzt(k), diff_k_in )

      ! LHS time tendency.
      lhs(k_tdiag,k)  &
      = real( lhs(k_tdiag,k) + ( 1.0 / dt ) )

      if ( l_stats_samp ) then

        ! Statistics:  implicit contributions for um or vm.

        if ( ium_ma + ivm_ma > 0 ) then
          if ( .not. l_implemented ) then
            tmp(1:3) &
            = term_ma_zt_lhs( wm_zt(k), gr%invrs_dzt(k), k )
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
          = 0.5 * invrs_rho_ds_zt(k)  &
          * diffusion_zt_lhs( rho_ds_zm(k) * Kh_zm(k),  &
                              rho_ds_zm(km1) * Kh_zm(km1), 0.0,  &
                              gr%invrs_dzm(km1), gr%invrs_dzm(k), gr%invrs_dzt(k), diff_k_in )
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
      + invrs_rho_ds_zt(2)  &
        * gr%invrs_dzt(2)  &
          * rho_ds_zm(1) * ( u_star_sqd / wind_speed(2) )

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
          - invrs_rho_ds_zt(2)  &
            * gr%invrs_dzt(2)  &
              * rho_ds_zm(1) * ( u_star_sqd / wind_speed(2) )
        endif

      endif ! l_stats_samp

    endif ! l_imp_sfc_momentum_flux


    return
  end subroutine windm_edsclrm_lhs

  !=============================================================================
  function windm_edsclrm_rhs( solve_type, dt, Kh_zm, xm, xm_tndcy,  &
                              rho_ds_zm, invrs_rho_ds_zt,  &
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
      dt                 ! Model timestep                             [s]

    real, dimension(gr%nnzp), intent(in) :: &
      Kh_zm,           & ! Eddy diffusivity on momentum levels        [m^2/s]
      xm,              & ! Eddy-scalar variable, xm (thermo. levels)  [units vary]
      xm_tndcy,        & ! The explicit time-tendency acting on xm    [units vary]
      rho_ds_zm,       & ! Dry, static density on momentum levels     [kg/m^3]
      invrs_rho_ds_zt    ! Inv. dry, static density at thermo. levels [m^3/kg]

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
      = 0.5 * invrs_rho_ds_zt(k)  &
      * diffusion_zt_lhs( rho_ds_zm(k) * Kh_zm(k),  &
                          rho_ds_zm(km1) * Kh_zm(km1), 0.0,  &
                          gr%invrs_dzm(km1), gr%invrs_dzm(k), gr%invrs_dzt(k), diff_k_in )
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
      rhs(2)  &
      = rhs(2)  &
      + invrs_rho_ds_zt(2)  &
        * gr%invrs_dzt(2)  &
          * rho_ds_zm(1) * xpwp_sfc

      if ( l_stats_samp ) then

        ! Statistics:  explicit contributions for um or vm.

        ! xm term ta is modified at level 2 to include the effects of the
        ! surface flux.  In this case, this effects the explicit portion of
        ! the term (after stat_begin_update_pt has already been called at
        ! level 2); call stat_modify_pt.
        if ( ixm_ta > 0 ) then
          call stat_modify_pt( ixm_ta, 2,  &
                               + invrs_rho_ds_zt(2)  &
                                 * gr%invrs_dzt(2)  &
                                   * rho_ds_zm(1) * xpwp_sfc,  &
                               zt )
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
    = 0.5 * invrs_rho_ds_zt(k)  &
    * diffusion_zt_lhs( rho_ds_zm(k) * Kh_zm(k),  &
                        rho_ds_zm(km1) * Kh_zm(km1), 0.0,  &
                        gr%invrs_dzm(km1), gr%invrs_dzm(k), gr%invrs_dzt(k), k )
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
  elemental function xpwp_fnc( Kh_zm, xm, xmp1, invrs_dzm )

    ! Description:
    ! Compute x'w' from x<k>, x<k+1>, Kh and invrs_dzm

    ! References:
    ! None
    !-----------------------------------------------------------------------

    implicit none

    ! Input variables
    real, intent(in) :: &
      Kh_zm, & ! Eddy diff. (k momentum level)                 [m^2/s]
      xm,    & ! x (k thermo level)                            [units vary]
      xmp1,  & ! x (k+1 thermo level)                          [units vary]
      invrs_dzm      ! Inverse of the grid spacing (k thermo level)  [1/m]
    ! Output variable
    real :: &
      xpwp_fnc ! x'w'   [(units vary)(m/s)]

    !-----------------------------------------------------------------------
    ! --- Begin Code ---

    ! Solve for x'w' at all intermediate model levels.
    xpwp_fnc = Kh_zm * invrs_dzm * ( xmp1 - xm )

    return
  end function xpwp_fnc


!===============================================================================

end module advance_windm_edsclrm_module
