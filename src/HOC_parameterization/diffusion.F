        module diffusion

        contains

!===============================================================================

        pure function diffusion_tlev( Km, Kmm1, nu, 
     .                                dzmm1, dzm, dzt, level )
     .  result( lhs )

!       Description:  Vertical diffusion of a thermodynamic-level variable 
!                     using an implicit scheme.

!       References:
!       None
!-----------------------------------------------------------------------

        use grid_class

        implicit none

        ! Constant parameters
        integer, parameter ::
     .  kp1diag = 1,
     .  kdiag   = 2,
     .  km1diag = 3

        ! Input Variables
        real, intent(in) :: Km, Kmm1, nu ! Diffusion coefficients.
                                         ! Km & Kmm1 are at momentum levels.
                                         ! nu is a background constant.

        real, intent(in) ::
     .  dzmm1, dzm, dzt ! Inverse of the grid spacing [m]

        integer, intent(in) :: level ! thermo. level that the calc. is for.

        ! Output
        real, dimension(3) :: lhs

!-----------------------------------------------------------------------

        ! DESCRIPTION
        !
        ! This function handles the portion of any equation that is concerned 
        ! with eddy diffusion for a variable that is located on thermodynamic
        ! levels:
        !
        ! d(xx)/dt = ... + d/dz [ (K+nu) * ( dxx/dz ) ] + ...
        !
        ! The normal discretization is to use the values of the variable that 
        ! is located on thermodynamic levels, xx, from three levels -- the
        ! central level and one-level-above and one-level-below the central
        ! level.  These three values lead to two first derivatives centered at
        ! the intermediate momentum levels.  The coefficient of eddy 
        ! diffusivity, K, is also located on the intermediate momentum levels,
        ! and nu is a background constant coefficient of eddy diffusivity.   
        ! These two coefficent-times-first-derivative factors lead to one second
        ! derivative centered at the central thermodynamic level.  This is 
        ! illustrated in the diagram below:
        !
        ! ---xx(k+1)------------------------------------------------  t(k+1)
        !
        ! ===============(K+nu)*(dxx/dz)============================  m(k)
        !
        ! ---xx( k )-------------------------d/dz[(K+nu)*(dxx/dz)]--  t(k)
        !
        ! ===============(K+nu)*(dxx/dz)============================  m(k-1)
        !
        ! ---xx(k-1)------------------------------------------------  t(k-1)
        !
        ! In equation form, this looks like:
        !
        ! dzt(k) * [   ( K(k) + nu ) * dzm(k) * ( xx(k+1) - xx(k) ) 
        !            - ( K(k-1) + nu ) * dzm(k-1) * ( xx(k) - xx(k-1) )   ];
        !
        ! where dzt(k) and dzm(k) are 1/(delta z) at level k.
        !
        ! This is then moved to the left-hand side of the equation and
        ! rearranged to be placed onto a tridiagonal matrix.
        !
        ! The boundary conditions used in this function are zero-flux boundary
        ! conditions, meaning the derivative of xx, dxx/dz, equals 0 at both
        ! the upper boundary and the lower boundary.  Thus, the value of xx
        ! at the boundary level has the same value of xx at a "phantom point"
        ! just outside the boundary.  This is illustrated in the diagram below:
        !
        ! ----------------------xx(2)----------------------  level 2
        !
        !
        ! ----------------------xx(1)----------------------  level 1 (boundary)
        !
        !
        ! ----------"phantom point" xx(0) = xx(1)----------  
        !
        ! At the boundary, this is still discretized the same way as before, 
        ! yielding (as an example for the lower boundary):
        !
        ! dzt(1) * [   ( K(1) + nu ) * dzm(1) * ( xx(2) - xx(1) ) 
        !            - ( K(0) + nu ) * dzm(0) * ( xx(1) - xx(0) )   ];
        !
        ! but since xx(1) = xx(0):
        !
        ! dzt(1) * [   ( K(1) + nu ) * dzm(1) * ( xx(2) - xx(1) )   ];
        !
        ! This is then moved to the left-hand side of the equation and
        ! rearranged to be placed onto a tridiagonal matrix.  The same thing 
        ! is done for the upper boundary.
        !
        ! Some equations in the model use fixed-point boundary conditions
        ! rather than zero-flux boundary conditions, such as wp3.  For 
        ! equations requiring fixed-point boundary conditions, either:
        ! a) loop over all levels from the second-lowest to the second-highest,
        !    ignoring the boundary levels, and then set the fixed-point values
        !    in the subroutine where the equation is computed, or
        ! b) loop over all levels and then overwrite the results for the 
        !    boundary levels with fixed-point values in the subroutine where 
        !    the equation is computed.
        !
        ! Finally, it should be pointed out that this technique of discretizing
        ! the equations leads to conservative differencing (for the diffusional
        ! term), as long as the grid levels are equally spaced.  When
        ! conservative differencing is in place, the column totals in the LHS
        ! matrix (for the diffusional term) should be equal to 0.  This insures
        ! that the total amount of the quantity over the entire domain is being
        ! conserved, meaning that nothing is lost due to diffusional effects.
        ! This is illustrated in the matrix setup below:
        !
        ! Left-hand side matrix, diffusional term only; first four levels:
        !
        !   ----->
        !   | +dzt(k)*(K(k)+nu)*dzm(k)       -dzt(k)*(K(k)+nu)*dzm(k)                        0
        !   |
        !   | -dzt(k)*(K(k-1)+nu)*dzm(k-1)   +dzt(k)*[ (K(k)+nu)*dzm(k)
        !   |                                         +(K(k-1)+nu)dzm*(k-1) ]   -dzt(k)*(K(k)+nu)*dzm(k)
        !   |
        !   |            0                   -dzt(k)*(K(k-1)+nu)*dzm(k-1)       +dzt(k)*[ (K(k)+nu)*dzm(k)
        !   |                                                                            +(K(k-1)+nu)*dzm(k-1) ]
        !   |
        !   |            0                               0                      -dzt(k)*(K(k-1)+nu)*dzm(k-1)
        !   |
        !  \ /
        !
        ! Brian Griffin.  November 24, 2007.

!-----------------------------------------------------------------------

        if ( level == 1 ) then
           ! k = 1 (bottom level); lower boundary level.
           lhs(kp1diag) = - dzt * (Km+nu) * dzm
           lhs(kdiag)   = + dzt * (Km+nu) * dzm
           lhs(km1diag) = 0.0

        elseif ( level > 1 .and. level < gr%nnzp ) then
           ! Most of the interior model; normal conditions.
           lhs(kp1diag) = - dzt * (Km+nu) * dzm
           lhs(kdiag)   = + dzt * ( (Km+nu) * dzm + (Kmm1+nu) * dzmm1 )
           lhs(km1diag) = - dzt * (Kmm1+nu) * dzmm1

        elseif ( level == gr%nnzp ) then
           ! k = gr%nnzp (top level); upper boundary level.
           lhs(kp1diag) = 0.0
           lhs(kdiag)   = + dzt * (Kmm1+nu) * dzmm1
           lhs(km1diag) = - dzt * (Kmm1+nu) * dzmm1

        endif

        end function diffusion_tlev

!===============================================================================

        pure function diffusion_mlev( Kt, Ktp1, nu, 
     .                                dztp1, dzt, dzm, level )
     .  result( lhs )

!       Description:  Vertical diffusion of a momentum-level variable 
!                     using an implicit scheme.

!       References:
!       None
!-----------------------------------------------------------------------

        use grid_class

        implicit none

        ! Constant parameters
        integer, parameter ::
     .  kp1diag = 1,
     .  kdiag   = 2,
     .  km1diag = 3

        ! Input Variables
        real, intent(in) :: Kt, Ktp1, nu ! Diffusion coefficients.
                                         ! Kt & Ktp1 are at thermo. levels.
                                         ! nu is a background constant.

        real, intent(in) ::
     .  dztp1, dzt, dzm ! Inverse of the grid spacing [m]

        integer, intent(in) :: level ! momentum level that the calc. is for.

        ! Output
        real, dimension(3) :: lhs

!-----------------------------------------------------------------------

        ! DESCRIPTION
        !
        ! This function handles the portion of any equation that is concerned 
        ! with eddy diffusion for a variable that is located on momentum
        ! levels:
        !
        ! d(xx)/dt = ... + d/dz [ (K+nu) * ( dxx/dz ) ] + ...
        !
        ! The normal discretization is to use the values of the variable that 
        ! is located on momentum levels, xx, from three levels -- the central 
        ! level and one-level-above and one-level-below the central level.  
        ! These three values lead to two first derivatives centered at the 
        ! intermediate thermodynamic levels.  The coefficient of eddy 
        ! diffusivity, K, is also located on the intermediate thermodynamic 
        ! levels, and nu is a background constant coefficient of eddy 
        ! diffusivity.  These two coefficent-times-first-derivative factors 
        ! lead to one second derivative centered at the central momentum level. 
        ! This is illustrated in the diagram below:
        !
        ! ===xx(k+1)================================================  m(k+1)
        !
        ! ---------------(K+nu)*(dxx/dz)----------------------------  t(k+1)
        !
        ! ===xx( k )=========================d/dz[(K+nu)*(dxx/dz)]==  m(k)
        !
        ! ---------------(K+nu)*(dxx/dz)----------------------------  t(k)
        !
        ! ===xx(k-1)================================================  m(k-1)
        !
        ! In equation form, this looks like:
        !
        ! dzm(k) * [   ( K(k+1) + nu ) * dzt(k+1) * ( xx(k+1) - xx(k) ) 
        !            - ( K(k) + nu ) * dzt(k) * ( xx(k) - xx(k-1) )   ];
        !
        ! where dzt(k) and dzm(k) are 1/(delta z) at level k.
        !
        ! This is then moved to the left-hand side of the equation and
        ! rearranged to be placed onto a tridiagonal matrix.
        !
        ! The boundary conditions used in this function are zero-flux boundary
        ! conditions, meaning the derivative of xx, dxx/dz, equals 0 at both
        ! the upper boundary and the lower boundary.  Thus, the value of xx
        ! at the boundary level has the same value of xx at a "phantom point"
        ! just outside the boundary.  This is illustrated in the diagram below:
        !
        ! ----------------------xx(2)----------------------  level 2
        !
        !
        ! ----------------------xx(1)----------------------  level 1 (boundary)
        !
        !
        ! ----------"phantom point" xx(0) = xx(1)----------  
        !
        ! At the boundary, this is still discretized the same way as before, 
        ! yielding (as an example for the lower boundary):
        !
        ! dzm(1) * [   ( K(2) + nu ) * dzt(2) * ( xx(2) - xx(1) ) 
        !            - ( K(1) + nu ) * dzt(1) * ( xx(1) - xx(0) )   ];
        !
        ! but since xx(1) = xx(0):
        !
        ! dzm(1) * [   ( K(2) + nu ) * dzt(2) * ( xx(2) - xx(1) )   ];
        !
        ! This is then moved to the left-hand side of the equation and
        ! rearranged to be placed onto a tridiagonal matrix.  The same thing 
        ! is done for the upper boundary.
        !
        ! Some equations in the model use fixed-point boundary conditions
        ! rather than zero-flux boundary conditions, such as wpthlp, wprtp, wp2,
        ! thlp2, rtp2, rtpthlp, etc.  For equations requiring fixed-point 
        ! boundary conditions, either:
        ! a) loop over all levels from the second-lowest to the second-highest,
        !    ignoring the boundary levels, and then set the fixed-point values
        !    in the subroutine where the equation is computed, or
        ! b) loop over all levels and then overwrite the results for the 
        !    boundary levels with fixed-point values in the subroutine where 
        !    the equation is computed.
        !
        ! Finally, it should be pointed out that this technique of discretizing
        ! the equations leads to conservative differencing (for the diffusional
        ! term), as long as the grid levels are equally spaced.  When
        ! conservative differencing is in place, the column totals in the LHS
        ! matrix (for the diffusional term) should be equal to 0.  This insures
        ! that the total amount of the quantity over the entire domain is being
        ! conserved, meaning that nothing is lost due to diffusional effects.
        ! This is illustrated in the matrix setup below:
        !
        ! Left-hand side matrix, diffusional term only; first four levels:
        !
        !   ----->
        !   | +dzm(k)*(K(k+1)+nu)*dzt(k+1)   -dzm(k)*(K(k+1)+nu)*dzt(k+1)                    0
        !   |
        !   | -dzm(k)*(K(k)+nu)*dzt(k)       +dzm(k)*[ (K(k+1)+nu)*dzt(k+1)
        !   |                                         +(K(k)+nu)*dzt(k) ]       -dzm(k)*(K(k+1)+nu)*dzt(k+1)
        !   |
        !   |            0                   -dzm(k)*(K(k)+nu)*dzt(k)           +dzm(k)*[ (K(k+1)+nu)*dzt(k+1)
        !   |                                                                            +(K(k)+nu)*dzt(k) ]
        !   |
        !   |            0                               0                      -dzm(k)*(K(k)+nu)*dzt(k)
        !   |
        !  \ /
        !
        ! Brian Griffin.  November 25, 2007.

!-----------------------------------------------------------------------

        if ( level == 1 ) then
           ! k = 1; lower boundery level at surface.
           lhs(kp1diag) = + dzm * (Ktp1+nu) * dztp1
           lhs(kdiag)   = - dzm * (Ktp1+nu) * dztp1
           lhs(km1diag) = 0.0

        elseif ( level > 1 .and. level < gr%nnzp ) then
           ! Most of the interior model; normal conditions.
           lhs(kp1diag) = - dzm * (Ktp1+nu) * dztp1
           lhs(kdiag)   = + dzm * ( (Ktp1+nu) * dztp1 + (Kt+nu) * dzt )
           lhs(km1diag) = - dzm * (Kt+nu) * dzt

        elseif ( level == gr%nnzp ) then
           ! k = gr%nnzp (top level); upper boundary level.
           lhs(kp1diag) = 0.0
           lhs(kdiag)   = - dzm * (Kt+nu) * dzt
           lhs(km1diag) = + dzm * (Kt+nu) * dzt

        endif

        end function diffusion_mlev

!===============================================================================

        end module diffusion
