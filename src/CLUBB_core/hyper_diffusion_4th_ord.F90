!-----------------------------------------------------------------------
! $Id$
!===============================================================================
module hyper_diffusion_4th_ord

  ! Description:
  ! Module hyper_diffusion_4th_ord computes the 4th-order numerical diffusion
  ! for any equation to which it is applied.  Hyper-diffusion will only be
  ! called if the model flag l_hyper_dfsn is set to true.  Function
  ! hyper_dfsn_4th_ord_zt_lhs handles 4th-order hyper-diffusion for variables
  ! that reside on thermodynamic levels.  Function hyper_dfsn_4th_ord_zm_lhs
  ! handles 4th-order hyper-diffusion for variables that reside on momentum
  ! levels.  A special constant coefficient of 4th-order numerical diffusion,
  ! nu_hd (which is sent in this module as nu), is used and has units of m^4/s.

  implicit none

  private ! Default Scope

  public :: hyper_dfsn_4th_ord_zt_lhs,  &
            hyper_dfsn_4th_ord_zm_lhs

contains

  !=============================================================================
  pure function hyper_dfsn_4th_ord_zt_lhs( boundary_cond, nu, invrs_dzt,  &
                                           invrs_dzm, invrs_dzmm1, invrs_dztp1,  &
                                           invrs_dztm1, invrs_dzmp1, invrs_dzmm2, level )  &
  result( lhs )


    use grid_class, only:  &
        gr  ! Variable(s)   gr%nnzp

    implicit none

    ! Constant parameters
    integer, parameter ::  &
      kp2_tdiag = 1,  & ! Thermodynamic super-super diagonal index.
      kp1_tdiag = 2,  & ! Thermodynamic super diagonal index.
      k_tdiag   = 3,  & ! Thermodynamic main diagonal index.
      km1_tdiag = 4,  & ! Thermodynamic sub diagonal index.
      km2_tdiag = 5     ! Thermodynamic sub-sub diagonal index.

    ! Input Variables
    character (len=*), intent(in) :: &
      boundary_cond   ! Type of boundary conditions being used
                      ! ('zero-flux' or 'fixed-point').

    real, intent(in) ::  &
      nu,     & ! Constant coefficient of 4th-order numerical diffusion  [m^4/s]
      invrs_dzt,    & ! Inverse of grid spacing over thermodynamic level (k)   [1/m]
      invrs_dzm,    & ! Inverse of grid spacing over momentum level (k)        [1/m]
      invrs_dzmm1,  & ! Inverse of grid spacing over momentum level (k-1)      [1/m]
      invrs_dztp1,  & ! Inverse of grid spacing over thermodynamic level (k+1) [1/m]
      invrs_dztm1,  & ! Inverse of grid spacing over thermodynamic level (k-1) [1/m]
      invrs_dzmp1,  & ! Inverse of grid spacing over momentum level (k+1)      [1/m]
      invrs_dzmm2     ! Inverse of grid spacing over momentum level (k-2)      [1/m]

    integer, intent(in) ::  & 
      level     ! Thermodynamic level where calculation occurs.          [-]

    ! Return Variable
    real, dimension(5) :: lhs


    if ( level == 1 ) then

       ! Lowest level
       ! k = 1; lower boundery level at surface.
       ! Only relevant if zero-flux boundary conditions are used.

       if ( trim( boundary_cond ) == 'zero-flux' ) then

          ! Zero-flux boundary conditions

          ! Thermodynamic sub-sub diagonal: [ x var_zt(k-2,<t+1>) ]
          lhs(km2_tdiag) = 0.0

          ! Thermodynamic sub diagonal: [ x var_zt(k-1,<t+1>) ]
          lhs(km1_tdiag) = 0.0

          ! Thermodynamic main diagonal: [ x var_zt(k,<t+1>) ]
          lhs(k_tdiag)   = +nu*invrs_dzt*invrs_dzm*(invrs_dztp1*invrs_dzm + invrs_dzt*invrs_dzm)

          ! Thermodynamic super diagonal: [ x var_zt(k+1,<t+1>) ]
          lhs(kp1_tdiag) = -nu*invrs_dzt*invrs_dzm*( invrs_dztp1*(invrs_dzmp1 + invrs_dzm)  &
                                        +invrs_dzt*invrs_dzm )

          ! Thermodynamic super-super diagonal: [ x var_zt(k+2,<t+1>) ]
          lhs(kp2_tdiag) = +nu*invrs_dzt*invrs_dzm*invrs_dztp1*invrs_dzmp1

       elseif ( trim( boundary_cond ) == 'fixed-point' ) then

          ! Fixed-point boundary conditions
          ! The left-hand side matrix contributions from level 1 are
          ! over-written or set in the parent subroutine.

          ! Thermodynamic sub-sub diagonal: [ x var_zt(k-2,<t+1>) ]
          lhs(km2_tdiag) = 0.0

          ! Thermodynamic sub diagonal: [ x var_zt(k-1,<t+1>) ]
          lhs(km1_tdiag) = 0.0

          ! Thermodynamic main diagonal: [ x var_zt(k,<t+1>) ]
          lhs(k_tdiag)   = 0.0

          ! Thermodynamic super diagonal: [ x var_zt(k+1,<t+1>) ]
          lhs(kp1_tdiag) = 0.0

          ! Thermodynamic super-super diagonal: [ x var_zt(k+2,<t+1>) ]
          lhs(kp2_tdiag) = 0.0

       endif


    elseif ( level == 2 ) then

       ! Second-lowest level

       if ( trim( boundary_cond ) == 'zero-flux' ) then

          ! Zero-flux boundary conditions

          ! Thermodynamic sub-sub diagonal: [ x var_zt(k-2,<t+1>) ]
          lhs(km2_tdiag) = 0.0

          ! Thermodynamic sub diagonal: [ x var_zt(k-1,<t+1>) ]
          lhs(km1_tdiag) = -nu*invrs_dzt*( invrs_dzm*invrs_dzt*invrs_dzmm1  &
                                    +invrs_dzmm1*( invrs_dzt*invrs_dzmm1  &
                                            +invrs_dztm1*invrs_dzmm1 ) )

          ! Thermodynamic main diagonal: [ x var_zt(k,<t+1>) ]
          lhs(k_tdiag)   = +nu*invrs_dzt*( invrs_dzm*( invrs_dztp1*invrs_dzm  &
                                          +invrs_dzt*(invrs_dzm + invrs_dzmm1) )  &
                                    +invrs_dzmm1*( invrs_dzt*(invrs_dzm + invrs_dzmm1)  &
                                            +invrs_dztm1*invrs_dzmm1 ) )

          ! Thermodynamic super diagonal: [ x var_zt(k+1,<t+1>) ]
          lhs(kp1_tdiag) = -nu*invrs_dzt*( invrs_dzm*( invrs_dztp1*(invrs_dzmp1 + invrs_dzm)  &
                                          +invrs_dzt*invrs_dzm )  &
                                    +invrs_dzmm1*invrs_dzt*invrs_dzm )

          ! Thermodynamic super-super diagonal: [ x var_zt(k+2,<t+1>) ]
          lhs(kp2_tdiag) = +nu*invrs_dzt*invrs_dzm*invrs_dztp1*invrs_dzmp1

       elseif ( trim( boundary_cond ) == 'fixed-point' ) then

          ! Fixed-point boundary conditions

          ! Thermodynamic sub-sub diagonal: [ x var_zt(k-2,<t+1>) ]
          lhs(km2_tdiag) = 0.0

          ! Thermodynamic sub diagonal: [ x var_zt(k-1,<t+1>) ]
          lhs(km1_tdiag) = -nu*invrs_dzt*( invrs_dzm*invrs_dzt*invrs_dzmm1  &
                                    +invrs_dzmm1*invrs_dzt*invrs_dzmm1 )

          ! Thermodynamic main diagonal: [ x var_zt(k,<t+1>) ]
          lhs(k_tdiag)   = +nu*invrs_dzt*( invrs_dzm*( invrs_dztp1*invrs_dzm  &
                                          +invrs_dzt*(invrs_dzm + invrs_dzmm1) )  &
                                    +invrs_dzmm1*( invrs_dzt*(invrs_dzm + invrs_dzmm1) ) )

          ! Thermodynamic super diagonal: [ x var_zt(k+1,<t+1>) ]
          lhs(kp1_tdiag) = -nu*invrs_dzt*( invrs_dzm*( invrs_dztp1*(invrs_dzmp1 + invrs_dzm)  &
                                          +invrs_dzt*invrs_dzm )  &
                                    +invrs_dzmm1*invrs_dzt*invrs_dzm )

          ! Thermodynamic super-super diagonal: [ x var_zt(k+2,<t+1>) ]
          lhs(kp2_tdiag) = +nu*invrs_dzt*invrs_dzm*invrs_dztp1*invrs_dzmp1

       endif


    elseif ( level > 2 .and. level < gr%nnzp-1 ) then

       ! k > 2 and k < num_levels-1
       ! These interior level are not effected by boundary conditions.

       ! Thermodynamic sub-sub diagonal: [ x var_zt(k-2,<t+1>) ]
       lhs(km2_tdiag) = +nu*invrs_dzt*invrs_dzmm1*invrs_dztm1*invrs_dzmm2

       ! Thermodynamic sub diagonal: [ x var_zt(k-1,<t+1>) ]
       lhs(km1_tdiag) = -nu*invrs_dzt*( invrs_dzm*invrs_dzt*invrs_dzmm1  &
                                 +invrs_dzmm1*( invrs_dzt*invrs_dzmm1  &
                                         +invrs_dztm1*(invrs_dzmm1 + invrs_dzmm2) ) )

       ! Thermodynamic main diagonal: [ x var_zt(k,<t+1>) ]
       lhs(k_tdiag)   = +nu*invrs_dzt*( invrs_dzm*( invrs_dztp1*invrs_dzm  &
                                       +invrs_dzt*(invrs_dzm + invrs_dzmm1) )  &
                                 +invrs_dzmm1*( invrs_dzt*(invrs_dzm + invrs_dzmm1)  &
                                         +invrs_dztm1*invrs_dzmm1 ) )

       ! Thermodynamic super diagonal: [ x var_zt(k+1,<t+1>) ]
       lhs(kp1_tdiag) = -nu*invrs_dzt*( invrs_dzm*( invrs_dztp1*(invrs_dzmp1 + invrs_dzm)  &
                                       +invrs_dzt*invrs_dzm )  &
                                 +invrs_dzmm1*invrs_dzt*invrs_dzm )

       ! Thermodynamic super-super diagonal: [ x var_zt(k+2,<t+1>) ]
       lhs(kp2_tdiag) = +nu*invrs_dzt*invrs_dzm*invrs_dztp1*invrs_dzmp1


    elseif ( level == gr%nnzp-1 ) then

       ! Second-highest level

       if ( trim( boundary_cond ) == 'zero-flux' ) then

          ! Zero-flux boundary conditions

          ! Thermodynamic sub-sub diagonal: [ x var_zt(k-2,<t+1>) ]
          lhs(km2_tdiag) = +nu*invrs_dzt*invrs_dzmm1*invrs_dztm1*invrs_dzmm2

          ! Thermodynamic sub diagonal: [ x var_zt(k-1,<t+1>) ]
          lhs(km1_tdiag) = -nu*invrs_dzt*( invrs_dzm*invrs_dzt*invrs_dzmm1  &
                                    +invrs_dzmm1*( invrs_dzt*invrs_dzmm1  &
                                            +invrs_dztm1*(invrs_dzmm1 + invrs_dzmm2) ) )

          ! Thermodynamic main diagonal: [ x var_zt(k,<t+1>) ]
          lhs(k_tdiag)   = +nu*invrs_dzt*( invrs_dzm*( invrs_dztp1*invrs_dzm  &
                                          +invrs_dzt*(invrs_dzm + invrs_dzmm1) )  &
                                    +invrs_dzmm1*( invrs_dzt*(invrs_dzm + invrs_dzmm1)  &
                                            +invrs_dztm1*invrs_dzmm1 ) )

          ! Thermodynamic super diagonal: [ x var_zt(k+1,<t+1>) ]
          lhs(kp1_tdiag) = -nu*invrs_dzt*( invrs_dzm*( invrs_dztp1*invrs_dzm  &
                                          +invrs_dzt*invrs_dzm )  &
                                    +invrs_dzmm1*invrs_dzt*invrs_dzm )

          ! Thermodynamic super-super diagonal: [ x var_zt(k+2,<t+1>) ]
          lhs(kp2_tdiag) = 0.0

       elseif ( trim( boundary_cond ) == 'fixed-point' ) then

          ! Fixed-point boundary conditions

          ! Thermodynamic sub-sub diagonal: [ x var_zt(k-2,<t+1>) ]
          lhs(km2_tdiag) = +nu*invrs_dzt*invrs_dzmm1*invrs_dztm1*invrs_dzmm2

          ! Thermodynamic sub diagonal: [ x var_zt(k-1,<t+1>) ]
          lhs(km1_tdiag) = -nu*invrs_dzt*( invrs_dzm*invrs_dzt*invrs_dzmm1  &
                                    +invrs_dzmm1*( invrs_dzt*invrs_dzmm1  &
                                            +invrs_dztm1*(invrs_dzmm1 + invrs_dzmm2) ) )

          ! Thermodynamic main diagonal: [ x var_zt(k,<t+1>) ]
          lhs(k_tdiag)   = +nu*invrs_dzt*( invrs_dzm*( invrs_dzt*(invrs_dzm + invrs_dzmm1) )  &
                                    +invrs_dzmm1*( invrs_dzt*(invrs_dzm + invrs_dzmm1)  &
                                            +invrs_dztm1*invrs_dzmm1 ) )

          ! Thermodynamic super diagonal: [ x var_zt(k+1,<t+1>) ]
          lhs(kp1_tdiag) = -nu*invrs_dzt*( invrs_dzm*invrs_dzt*invrs_dzm  &
                                    +invrs_dzmm1*invrs_dzt*invrs_dzm )

          ! Thermodynamic super-super diagonal: [ x var_zt(k+2,<t+1>) ]
          lhs(kp2_tdiag) = 0.0

       endif


    elseif ( level == gr%nnzp ) then

       ! Highest level
       ! k = gr%nnzp; upper boundery level at model top.
       ! Only relevant if zero-flux boundary conditions are used.

       if ( trim( boundary_cond ) == 'zero-flux' ) then

          ! Zero-flux boundary conditions

          ! Thermodynamic sub-sub diagonal: [ x var_zt(k-2,<t+1>) ]
          lhs(km2_tdiag) = +nu*invrs_dzt*invrs_dzmm1*invrs_dztm1*invrs_dzmm2

          ! Thermodynamic sub diagonal: [ x var_zt(k-1,<t+1>) ]
          lhs(km1_tdiag) = -nu*invrs_dzt*invrs_dzmm1*( invrs_dzt*invrs_dzmm1  &
                                          +invrs_dztm1*(invrs_dzmm1 + invrs_dzmm2) )

          ! Thermodynamic main diagonal: [ x var_zt(k,<t+1>) ]
          lhs(k_tdiag)   = +nu*invrs_dzt*invrs_dzmm1* &
                             (invrs_dzt*invrs_dzmm1 + invrs_dztm1*invrs_dzmm1)

          ! Thermodynamic super diagonal: [ x var_zt(k+1,<t+1>) ]
          lhs(kp1_tdiag) = 0.0

          ! Thermodynamic super-super diagonal: [ x var_zt(k+2,<t+1>) ]
          lhs(kp2_tdiag) = 0.0

       elseif ( trim( boundary_cond ) == 'fixed-point' ) then

          ! Fixed-point boundary conditions
          ! The left-hand side matrix contributions from level gr%nnzp are
          ! over-written or set in the parent subroutine.

          ! Thermodynamic sub-sub diagonal: [ x var_zt(k-2,<t+1>) ]
          lhs(km2_tdiag) = 0.0

          ! Thermodynamic sub diagonal: [ x var_zt(k-1,<t+1>) ]
          lhs(km1_tdiag) = 0.0

          ! Thermodynamic main diagonal: [ x var_zt(k,<t+1>) ]
          lhs(k_tdiag)   = 0.0

          ! Thermodynamic super diagonal: [ x var_zt(k+1,<t+1>) ]
          lhs(kp1_tdiag) = 0.0

          ! Thermodynamic super-super diagonal: [ x var_zt(k+2,<t+1>) ]
          lhs(kp2_tdiag) = 0.0

       endif

    endif

    return

  end function hyper_dfsn_4th_ord_zt_lhs

  !=============================================================================
  pure function hyper_dfsn_4th_ord_zm_lhs( boundary_cond, nu, invrs_dzm,  &
                                           invrs_dztp1, invrs_dzt, invrs_dzmp1,  &
                                           invrs_dzmm1, dztp2, invrs_dztm1, level )  &
  result( lhs )


    use grid_class, only:  &
        gr  ! Variable(s)   gr%nnzp

    implicit none

    ! Constant parameters
    integer, parameter ::  &
      kp2_mdiag = 1,  & ! Momentum super-super diagonal index.
      kp1_mdiag = 2,  & ! Momentum super diagonal index.
      k_mdiag   = 3,  & ! Momentum main diagonal index.
      km1_mdiag = 4,  & ! Momentum sub diagonal index.
      km2_mdiag = 5     ! Momentum sub-sub diagonal index.

    ! Input Variables
    character (len=*), intent(in) :: &
      boundary_cond   ! Type of boundary conditions being used
                      ! ('zero-flux' or 'fixed-point').

    real, intent(in) ::  &
      nu,     & ! Constant coefficient of 4th-order numerical diffusion  [m^4/s]
      invrs_dzm,    & ! Inverse of grid spacing over momentum level (k)        [1/m]
      invrs_dztp1,  & ! Inverse of grid spacing over thermodynamic level (k+1) [1/m]
      invrs_dzt,    & ! Inverse of grid spacing over thermodynamic level (k)   [1/m]
      invrs_dzmp1,  & ! Inverse of grid spacing over momentum level (k+1)      [1/m]
      invrs_dzmm1,  & ! Inverse of grid spacing over momentum level (k-1)      [1/m]
      dztp2,  & ! Inverse of grid spacing over thermodynamic level (k+2) [1/m]
      invrs_dztm1     ! Inverse of grid spacing over thermodynamic level (k-1) [1/m]

    integer, intent(in) ::  & 
      level     ! Momentum level where calculation occurs.               [-]

    ! Return Variable
    real, dimension(5) :: lhs


    if ( level == 1 ) then

       ! Lowest level
       ! k = 1; lower boundery level at surface.
       ! Only relevant if zero-flux boundary conditions are used.

       if ( trim( boundary_cond ) == 'zero-flux' ) then

          ! Zero-flux boundary conditions

          ! Momentum sub-sub diagonal: [ x var_zm(k-2,<t+1>) ]
          lhs(km2_mdiag) = 0.0

          ! Momentum sub diagonal: [ x var_zm(k-1,<t+1>) ]
          lhs(km1_mdiag) = 0.0

          ! Momentum main diagonal: [ x var_zm(k,<t+1>) ]
          lhs(k_mdiag)   = +nu*invrs_dzm*invrs_dztp1* &
                             (invrs_dzmp1*invrs_dztp1 + invrs_dzm*invrs_dztp1)

          ! Momentum super diagonal: [ x var_zm(k+1,<t+1>) ]
          lhs(kp1_mdiag) = -nu*invrs_dzm*invrs_dztp1*( invrs_dzmp1*(dztp2 + invrs_dztp1)  &
                                          +invrs_dzm*invrs_dztp1 )

          ! Momentum super-super diagonal: [ x var_zm(k+2,<t+1>) ]
          lhs(kp2_mdiag) = +nu*invrs_dzm*invrs_dztp1*invrs_dzmp1*dztp2

       elseif ( trim( boundary_cond ) == 'fixed-point' ) then

          ! Fixed-point boundary conditions
          ! The left-hand side matrix contributions from level 1 are
          ! over-written or set in the parent subroutine.

          ! Momentum sub-sub diagonal: [ x var_zm(k-2,<t+1>) ]
          lhs(km2_mdiag) = 0.0

          ! Momentum sub diagonal: [ x var_zm(k-1,<t+1>) ]
          lhs(km1_mdiag) = 0.0

          ! Momentum main diagonal: [ x var_zm(k,<t+1>) ]
          lhs(k_mdiag)   = 0.0

          ! Momentum super diagonal: [ x var_zm(k+1,<t+1>) ]
          lhs(kp1_mdiag) = 0.0

          ! Momentum super-super diagonal: [ x var_zm(k+2,<t+1>) ]
          lhs(kp2_mdiag) = 0.0

       endif


    elseif ( level == 2 ) then

       ! Second-lowest level

       if ( trim( boundary_cond ) == 'zero-flux' ) then

          ! Zero-flux boundary conditions

          ! Momentum sub-sub diagonal: [ x var_zm(k-2,<t+1>) ]
          lhs(km2_mdiag) = 0.0

          ! Momentum sub diagonal: [ x var_zm(k-1,<t+1>) ]
          lhs(km1_mdiag) = -nu*invrs_dzm*( invrs_dztp1*invrs_dzm*invrs_dzt  &
                                    +invrs_dzt*( invrs_dzm*invrs_dzt  &
                                          +invrs_dzmm1*invrs_dzt ) )

          ! Momentum main diagonal: [ x var_zm(k,<t+1>) ]
          lhs(k_mdiag)   = +nu*invrs_dzm*( invrs_dztp1*( invrs_dzmp1*invrs_dztp1  &
                                            +invrs_dzm*(invrs_dztp1 + invrs_dzt) )  &
                                    +invrs_dzt*( invrs_dzm*(invrs_dztp1 + invrs_dzt)  &
                                          +invrs_dzmm1*invrs_dzt ) )

          ! Momentum super diagonal: [ x var_zm(k+1,<t+1>) ]
          lhs(kp1_mdiag) = -nu*invrs_dzm*( invrs_dztp1*( invrs_dzmp1*(dztp2 + invrs_dztp1)  &
                                            +invrs_dzm*invrs_dztp1 )  &
                                    +invrs_dzt*invrs_dzm*invrs_dztp1 )

          ! Momentum super-super diagonal: [ x var_zm(k+2,<t+1>) ]
          lhs(kp2_mdiag) = +nu*invrs_dzm*invrs_dztp1*invrs_dzmp1*dztp2

       elseif ( trim( boundary_cond ) == 'fixed-point' ) then

          ! Fixed-point boundary conditions

          ! Momentum sub-sub diagonal: [ x var_zm(k-2,<t+1>) ]
          lhs(km2_mdiag) = 0.0

          ! Momentum sub diagonal: [ x var_zm(k-1,<t+1>) ]
          lhs(km1_mdiag) = -nu*invrs_dzm*( invrs_dztp1*invrs_dzm*invrs_dzt  &
                                    +invrs_dzt*invrs_dzm*invrs_dzt )

          ! Momentum main diagonal: [ x var_zm(k,<t+1>) ]
          lhs(k_mdiag)   = +nu*invrs_dzm*( invrs_dztp1*( invrs_dzmp1*invrs_dztp1  &
                                            +invrs_dzm*(invrs_dztp1 + invrs_dzt) )  &
                                    +invrs_dzt*invrs_dzm*(invrs_dztp1 + invrs_dzt) )

          ! Momentum super diagonal: [ x var_zm(k+1,<t+1>) ]
          lhs(kp1_mdiag) = -nu*invrs_dzm*( invrs_dztp1*( invrs_dzmp1*(dztp2 + invrs_dztp1)  &
                                            +invrs_dzm*invrs_dztp1 )  &
                                    +invrs_dzt*invrs_dzm*invrs_dztp1 )

          ! Momentum super-super diagonal: [ x var_zm(k+2,<t+1>) ]
          lhs(kp2_mdiag) = +nu*invrs_dzm*invrs_dztp1*invrs_dzmp1*dztp2
   
       endif


    elseif ( level > 2 .and. level < gr%nnzp-1 ) then

       ! k > 2 and k < num_levels-1
       ! These interior level are not effected by boundary conditions.

       ! Momentum sub-sub diagonal: [ x var_zm(k-2,<t+1>) ]
       lhs(km2_mdiag) = +nu*invrs_dzm*invrs_dzt*invrs_dzmm1*invrs_dztm1

       ! Momentum sub diagonal: [ x var_zm(k-1,<t+1>) ]
       lhs(km1_mdiag) = -nu*invrs_dzm*( invrs_dztp1*invrs_dzm*invrs_dzt  &
                                 +invrs_dzt*( invrs_dzm*invrs_dzt  &
                                       +invrs_dzmm1*(invrs_dzt + invrs_dztm1) ) )

       ! Momentum main diagonal: [ x var_zm(k,<t+1>) ]
       lhs(k_mdiag)   = +nu*invrs_dzm*( invrs_dztp1*( invrs_dzmp1*invrs_dztp1  &
                                         +invrs_dzm*(invrs_dztp1 + invrs_dzt) )  &
                                 +invrs_dzt*( invrs_dzm*(invrs_dztp1 + invrs_dzt)  &
                                       +invrs_dzmm1*invrs_dzt ) )

       ! Momentum super diagonal: [ x var_zm(k+1,<t+1>) ]
       lhs(kp1_mdiag) = -nu*invrs_dzm*( invrs_dztp1*( invrs_dzmp1*(dztp2 + invrs_dztp1)  &
                                         +invrs_dzm*invrs_dztp1 )  &
                                 +invrs_dzt*invrs_dzm*invrs_dztp1 )

       ! Momentum super-super diagonal: [ x var_zm(k+2,<t+1>) ]
       lhs(kp2_mdiag) = +nu*invrs_dzm*invrs_dztp1*invrs_dzmp1*dztp2
   

    elseif ( level == gr%nnzp-1 ) then

       ! Second-highest level

       if ( trim( boundary_cond ) == 'zero-flux' ) then

          ! Zero-flux boundary conditions

          ! Momentum sub-sub diagonal: [ x var_zm(k-2,<t+1>) ]
          lhs(km2_mdiag) = +nu*invrs_dzm*invrs_dzt*invrs_dzmm1*invrs_dztm1

          ! Momentum sub diagonal: [ x var_zm(k-1,<t+1>) ]
          lhs(km1_mdiag) = -nu*invrs_dzm*( invrs_dztp1*invrs_dzm*invrs_dzt  &
                                    +invrs_dzt*( invrs_dzm*invrs_dzt  &
                                          +invrs_dzmm1*(invrs_dzt + invrs_dztm1) ) )

          ! Momentum main diagonal: [ x var_zm(k,<t+1>) ]
          lhs(k_mdiag)   = +nu*invrs_dzm*( invrs_dztp1*( invrs_dzmp1*invrs_dztp1  &
                                            +invrs_dzm*(invrs_dztp1 + invrs_dzt) )  &
                                    +invrs_dzt*( invrs_dzm*(invrs_dztp1 + invrs_dzt)  &
                                          +invrs_dzmm1*invrs_dzt ) )

          ! Momentum super diagonal: [ x var_zm(k+1,<t+1>) ]
          lhs(kp1_mdiag) = -nu*invrs_dzm*( invrs_dztp1*( invrs_dzmp1*invrs_dztp1  &
                                            +invrs_dzm*invrs_dztp1 )  &
                                    +invrs_dzt*invrs_dzm*invrs_dztp1 )

          ! Momentum super-super diagonal: [ x var_zm(k+2,<t+1>) ]
          lhs(kp2_mdiag) = 0.0

       elseif ( trim( boundary_cond ) == 'fixed-point' ) then

          ! Fixed-point boundary conditions

          ! Momentum sub-sub diagonal: [ x var_zm(k-2,<t+1>) ]
          lhs(km2_mdiag) = +nu*invrs_dzm*invrs_dzt*invrs_dzmm1*invrs_dztm1

          ! Momentum sub diagonal: [ x var_zm(k-1,<t+1>) ]
          lhs(km1_mdiag) = -nu*invrs_dzm*( invrs_dztp1*invrs_dzm*invrs_dzt  &
                                    +invrs_dzt*( invrs_dzm*invrs_dzt  &
                                          +invrs_dzmm1*(invrs_dzt + invrs_dztm1) ) )

          ! Momentum main diagonal: [ x var_zm(k,<t+1>) ]
          lhs(k_mdiag)   = +nu*invrs_dzm*( invrs_dztp1*( invrs_dzm*(invrs_dztp1 + invrs_dzt) )  &
                                    +invrs_dzt*( invrs_dzm*(invrs_dztp1 + invrs_dzt)  &
                                          +invrs_dzmm1*invrs_dzt ) )

          ! Momentum super diagonal: [ x var_zm(k+1,<t+1>) ]
          lhs(kp1_mdiag) = -nu*invrs_dzm*( invrs_dztp1*invrs_dzm*invrs_dztp1  &
                                    +invrs_dzt*invrs_dzm*invrs_dztp1 )

          ! Momentum super-super diagonal: [ x var_zm(k+2,<t+1>) ]
          lhs(kp2_mdiag) = 0.0

       endif


    elseif ( level == gr%nnzp ) then

       ! Highest level
       ! k = gr%nnzp; upper boundery level at model top.
       ! Only relevant if zero-flux boundary conditions are used.

       if ( trim( boundary_cond ) == 'zero-flux' ) then

          ! Zero-flux boundary conditions

          ! Momentum sub-sub diagonal: [ x var_zm(k-2,<t+1>) ]
          lhs(km2_mdiag) = +nu*invrs_dzm*invrs_dzt*invrs_dzmm1*invrs_dztm1

          ! Momentum sub diagonal: [ x var_zm(k-1,<t+1>) ]
          lhs(km1_mdiag) = -nu*invrs_dzm*invrs_dzt*( invrs_dzm*invrs_dzt  &
                                        +invrs_dzmm1*(invrs_dzt + invrs_dztm1) )

          ! Momentum main diagonal: [ x var_zm(k,<t+1>) ]
          lhs(k_mdiag)   = +nu*invrs_dzm*invrs_dzt*(invrs_dzm*invrs_dzt + invrs_dzmm1*invrs_dzt)

          ! Momentum super diagonal: [ x var_zm(k+1,<t+1>) ]
          lhs(kp1_mdiag) = 0.0

          ! Momentum super-super diagonal: [ x var_zm(k+2,<t+1>) ]
          lhs(kp2_mdiag) = 0.0

       elseif ( trim( boundary_cond ) == 'fixed-point' ) then

          ! Fixed-point boundary conditions
          ! The left-hand side matrix contributions from level gr%nnzp are
          ! over-written or set in the parent subroutine.

          ! Momentum sub-sub diagonal: [ x var_zm(k-2,<t+1>) ]
          lhs(km2_mdiag) = 0.0

          ! Momentum sub diagonal: [ x var_zm(k-1,<t+1>) ]
          lhs(km1_mdiag) = 0.0

          ! Momentum main diagonal: [ x var_zm(k,<t+1>) ]
          lhs(k_mdiag)   = 0.0

          ! Momentum super diagonal: [ x var_zm(k+1,<t+1>) ]
          lhs(kp1_mdiag) = 0.0

          ! Momentum super-super diagonal: [ x var_zm(k+2,<t+1>) ]
          lhs(kp2_mdiag) = 0.0

       endif


    endif


    return

  end function hyper_dfsn_4th_ord_zm_lhs

!===============================================================================

end module hyper_diffusion_4th_ord
