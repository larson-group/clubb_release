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

!        if ( level == 1 ) then
!           ! k = 1; below surface; no effects.
!           lhs(kp1diag) = 0.0
!           lhs(kdiag)   = 0.0
!           lhs(km1diag) = 0.0
!
        if ( level == 1 ) then
!        elseif ( level == 2 ) then
           ! k = 2; lower boundary level.
           lhs(kp1diag) = - dzt * (Km+nu) * dzm
           lhs(kdiag)   = + dzt * (Km+nu) * dzm
           lhs(km1diag) = 0.0

!        elseif ( level > 2 .and. level < gr%nnzp ) then
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
