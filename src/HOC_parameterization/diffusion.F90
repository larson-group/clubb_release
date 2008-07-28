! $Id: diffusion.F90,v 1.3 2008-07-28 19:34:42 faschinj Exp $
!===============================================================================
module diffusion

!       Description:
!       Module diffusion computes the eddy diffusion terms for all of the
!       time-tendency (prognostic) equations in the CLUBB parameterization.
!       Most of the eddy diffusion terms are solved for completely implicitly,
!       and therefore become part of the left-hand side of their respective
!       equations.  However, wp2 and wp3 have an option to use a Crank-Nicholson
!       eddy diffusion scheme, which has both implicit and explicit components.
!
!       Function diffusion_zt_lhs handles the eddy diffusion terms for the
!       variables located at thermodynamic grid levels.  These variables are:
!       wp3 and all hydrometeor species.
!
!       Function diffusion_zm_lhs handles the eddy diffusion terms for the
!       variables located at momentum grid levels.  The variables are:
!       wprtp, wpthlp, wp2, rtp2, thlp2, rtpthlp, up2, vp2, wpsclrp, sclrprtp,
!       sclrpthlp, and sclrp2.

implicit none

private ! Default Scope

public :: diffusion_zt_lhs, & 
          diffusion_zm_lhs

contains

!===============================================================================
pure function diffusion_zt_lhs( K_m, K_mm1, nu,  & 
                                dzmm1, dzm, dzt, level ) & 
result( lhs )

!       Description:
!       Vertical eddy diffusion of var_t:  implicit portion of the code.
!
!       The variable "var_t" stands for a variable that is located at
!       thermodynamic grid levels.
!
!       The d(var_t)/dt equation contains an eddy diffusion term:
!
!       + d [ ( K_m + nu ) * d(var_t)/dz ] / dz.
!
!       This term is usually solved for completely implicitly, such that:
!
!       + d [ ( K_m + nu ) * d( var_t(t+1) )/dz ] / dz.
!
!       However, when a Crank-Nicholson scheme is used, the eddy 
!       diffusion term has both implicit and explicit components, 
!       such that:
!
!       + (1/2) * d [ ( K_m + nu ) * d( var_t(t+1) )/dz ] / dz
!          + (1/2) * d [ ( K_m + nu ) * d( var_t(t) )/dz ] / dz;
!
!       for which the implicit component is:
!
!       + (1/2) * d [ ( K_m + nu ) * d( var_t(t+1) )/dz ] / dz.
!
!       Note:  When the implicit term is brought over to the left-hand
!              side, the sign is reversed and the leading "+" in front
!              of the term is changed to a "-".
!
!       Timestep index (t) stands for the index of the current timestep,
!       while timestep index (t+1) stands for the index of the next
!       timestep, which is being advanced to in solving the d(var_t)/dt
!       equation.
!
!       The implicit portion of this term is discretized as follows:
!
!       The values of var_t are found on the thermodynamic levels, while
!       the values of K_m are found on the momentum levels.  The 
!       derivatives (d/dz) of var_t are taken over the intermediate 
!       momentum levels.  At the intermediate momentum levels, 
!       d(var_t)/dz is multiplied by ( K_m + nu ).  Then, the derivative
!       of the whole mathematical expression is taken over the central 
!       thermodynamic level, which yields the desired result.
!
!       ---var_tp1----------------------------------------------- t(k+1)
!
!       ===========d(var_t)/dz==(K_m+nu)========================= m(k)
!
!       ---var_t--------------------d[(K_m+nu)*d(var_t)/dz]/dz--- t(k)
!
!       ===========d(var_t)/dz==(K_mm1+nu)======================= m(k-1)
!
!       ---var_tm1----------------------------------------------- t(k-1)
!
!       The vertical indices t(k+1), m(k), t(k), m(k-1), and t(k-1)
!       correspond with altitudes zt(k+1), zm(k), zt(k), zm(k-1),
!       and zt(k-1), respectively.  The letter "t" is used for
!       thermodynamic levels and the letter "m" is used for momentum
!       levels.
!
!       dzt(k)   = 1 / ( zm(k) - zm(k-1) )
!       dzm(k)   = 1 / ( zt(k+1) - zt(k) )
!       dzm(k-1) = 1 / ( zt(k) - zt(k-1) )
!
!       Note:  This function only computes the general implicit form:
!              + d [ ( K_m + nu ) * d( var_t(t+1) )/dz ] / dz.
!              For a Crank-Nicholson scheme, the left-hand side result
!              of this function will have to be multiplied by (1/2).
!              For a Crank-Nicholson scheme, the right-hand side 
!              (explicit) component needs to be computed by multiplying 
!              the left-hand side results by (1/2), reversing the sign
!              on each left-hand side element, and then multiplying 
!              each element by the appropriate var_t(t) value from the
!              appropriate vertical level.
!
!
!       Boundary Conditions:
!
!       1) Zero-flux boundary conditions.
!          This function is set up to use a zero-flux boundary condition 
!          at both the lower boundary level (k=1) and the upper boundary 
!          level (k=gr%nnzp).  This means that the derivative of var_t, 
!          d(var_t)/dz, equals 0 at both the lower boundary and the upper
!          boundary.
!
!          In order to discretize the lower boundary condition, consider
!          a new level outside the model (thermodynamic level 0) just 
!          below the lower boundary level (thermodynamic level 1).  The 
!          value of var_t at the level just outside the model is defined 
!          to be the same as the value of var_t at the lower boundary 
!          level.  Therefore, the value of d(var_t)/dz between the level 
!          just outside the model and the lower boundary level is 0, 
!          satisfying the zero-flux boundary condition.  The other value
!          for d(var_t)/dz (between thermodynamic level 2 and 
!          thermodynamic level 1) is taken over the intermediate momentum
!          level (momentum level 1), where it is multiplied by the factor
!          ( K_m(1) + nu ).  Then, the derivative of the whole expression
!          is taken over the central thermodynamic level.
!
!          ---var_t(2)------------------------------------------- t(2)
!
!          ============d(var_t)/dz==(K_m(1)+nu)================== m(1)
!
!          ---var_t(1)---------------d[(K_m+nu)*d(var_t)/dz]/dz-- t(1) Boundary
!
!                      [d(var_t)/dz = 0]
!
!          ---[var_t(0) = var_t(1)]----(level outside model)----- t(0)
!
!          The result is dependent only on values of K_m found at 
!          momentum level 1 and values of var_t found at thermodynamic 
!          levels 1 and 2.  Thus, it only affects 2 diagonals on the
!          left-hand side matrix.
!
!          The same method can be used to discretize the upper boundary
!          by considering a new level outside the model just above the 
!          upper boundary level.
!
!       2) Fixed-point boundary conditions.
!          Many equations in the model use fixed-point boundary 
!          conditions rather than zero-flux boundary conditions.  This
!          means that the value of var_t stays the same over the course
!          of the timestep at the lower boundary, as well as at the 
!          upper boundary.
!
!          In order to discretize the boundary conditions for equations 
!          requiring fixed-point boundary conditions, either:
!          a) in the parent subroutine or function (that calls this 
!             function), loop over all vertical levels from the 
!             second-lowest to the second-highest, ignoring the boundary
!             levels.  Then set the values at the boundary levels in the
!             parent subroutine; or
!          b) in the parent subroutine or function, loop over all 
!             vertical levels and then overwrite the results at the 
!             boundary levels.
!
!          Either way, at the boundary levels, an array with a value 
!          of 1 at the main diagonal on the left-hand side and with 
!          values of 0 at all other diagonals on the left-hand side will
!          preserve the right-hand side value at that level, thus 
!          satisfying the fixed-point boundary conditions.
!
!
!       Conservation Properties:
!
!       When zero-flux boundary conditions are used, this technique of
!       discretizing the eddy diffusion term leads to conservative 
!       differencing.  When conservative differencing is in place, the
!       column totals for each column in the left-hand side matrix 
!       (for the eddy diffusion term) should be equal to 0.  This 
!       ensures that the total amount of the quantity var_t over the 
!       entire vertical domain is being conserved, meaning that 
!       nothing is lost due to diffusional effects.
!
!       To see that this conservation law is satisfied, compute the 
!       eddy diffusion of var_t and integrate vertically.  In 
!       discretized matrix notation (where "i" stands for the matrix 
!       column and "j" stands for the matrix row):
!
!        0 =
!        Sum_j Sum_i ( 1/dzt )_i ( dzt * ((K_m+nu)*dzm) )_ij (var_t)_j.
!
!       The left-hand side matrix, ( dzt * ((K_m+nu)*dzm) )_ij, is
!       partially written below.  The sum over i in the above equation 
!       removes dzt everywhere from the matrix below.  The sum over j
!       leaves the column totals that are desired.
!
!       Left-hand side matrix contributions from eddy diffusion term;
!       first four vertical levels:
!
!     -------------------------------------------------------------------------------------------->
!k=1 | +dzt(k)*                   -dzt(k)*(K_m(k)+nu)*dzm(k)                     0
!    |   (K_m(k)+nu)*dzm(k)
!    | 
!k=2 | -dzt(k)*                   +dzt(k)*[ (K_m(k)+nu)*dzm(k)      -dzt(k)*(K_m(k)+nu)*dzm(k)
!    |    (K_m(k-1)+nu)*dzm(k-1)                                      +(K_m(k-1)+nu)*dzm(k-1) ]
!    |
!k=3 |              0             -dzt(k)*(K_m(k-1)+nu)*dzm(k-1)    +dzt(k)*[ (K_m(k)+nu)*dzm(k)
!    |                                                                +(K_m(k-1)+nu)*dzm(k-1) ]
!    |
!k=4 |              0                           0                    -dzt(k)*(K_m(k-1)+nu)*dzm(k-1)
!    |
!   \ /
!
!       Note:  The superdiagonal term from level 3 and both the main 
!              diagonal and superdiagonal terms from level 4 are not 
!              shown on this diagram.
!
!       Note:  The matrix shown is a tridiagonal matrix.  For a band 
!              diagonal matrix (with 5 diagonals), there would be an 
!              extra row between each of the rows shown and an extra
!              column between each of the columns shown.  However, 
!              for the purposes of the var_t eddy diffusion term, those
!              extra row and column values are all 0, and the 
!              conservation properties of the matrix aren't effected.
!              
!       If fixed-point boundary conditions are used, the matrix 
!       entries at level 1 (k=1) read:  1   0   0; which means that 
!       conservative differencing is not in play.  The total amount of 
!       var_t over the entire vertical domain is not being conserved, 
!       as amounts of var_t may be fluxed out through the upper boundary
!       or lower boundary through the effects of diffusion.
!
!       Brian Griffin.  April 26, 2008.

!       References:
!       None
!-----------------------------------------------------------------------

use grid_class, only: & 
    gr ! Variable(s)

implicit none

! Constant parameters
integer, parameter :: & 
kp1_tdiag = 1,    & ! Thermodynamic superdiagonal index.
k_tdiag   = 2,    & ! Thermodynamic main diagonal index.
km1_tdiag = 3    ! Thermodynamic subdiagonal index.

! Input Variables
real, intent(in) ::  & 
K_m,    & ! Coefficient of eddy diffusivity at momentum level (k)   [m^2/s]
K_mm1,  & ! Coefficient of eddy diffusivity at momentum level (k-1) [m^2/s]
nu,     & ! Background constant coefficient of eddy diffusivity     [m^2/s]
dzt,    & ! Inverse of grid spacing over thermodynamic level (k)    [1/m]
dzm,    & ! Inverse of grid spacing over momentum level (k)         [1/m]
dzmm1  ! Inverse of grid spacing over momentum level (k-1)       [1/m]

integer, intent(in) ::  & 
level ! Thermodynamic level where calculation occurs.           [-]

! Return Variable
real, dimension(3) :: lhs

if ( level == 1 ) then

   ! k = 1 (bottom level); lower boundary level.
   ! Only relevant if zero-flux boundary conditions are used.

   ! Thermodynamic superdiagonal: [ x var_t(k+1,<t+1>) ]
   lhs(kp1_tdiag) = - dzt * (K_m+nu) * dzm

   ! Thermodynamic main diagonal: [ x var_t(k,<t+1>) ]
   lhs(k_tdiag)   = + dzt * (K_m+nu) * dzm

   ! Thermodynamic subdiagonal: [ x var_t(k-1,<t+1>) ]
   lhs(km1_tdiag) = 0.0


elseif ( level > 1 .and. level < gr%nnzp ) then

   ! Most of the interior model; normal conditions.

   ! Thermodynamic superdiagonal: [ x var_t(k+1,<t+1>) ]
   lhs(kp1_tdiag) = - dzt * (K_m+nu) * dzm

   ! Thermodynamic main diagonal: [ x var_t(k,<t+1>) ]
   lhs(k_tdiag)   = + dzt * ( (K_m+nu)*dzm + (K_mm1+nu)*dzmm1 )

   ! Thermodynamic subdiagonal: [ x var_t(k-1,<t+1>) ]
   lhs(km1_tdiag) = - dzt * (K_mm1+nu) * dzmm1


elseif ( level == gr%nnzp ) then

   ! k = gr%nnzp (top level); upper boundary level.
   ! Only relevant if zero-flux boundary conditions are used.

   ! Thermodynamic superdiagonal: [ x var_t(k+1,<t+1>) ]
   lhs(kp1_tdiag) = 0.0

   ! Thermodynamic main diagonal: [ x var_t(k,<t+1>) ]
   lhs(k_tdiag)   = + dzt * (K_mm1+nu) * dzmm1

   ! Thermodynamic subdiagonal: [ x var_t(k-1,<t+1>) ]
   lhs(km1_tdiag) = - dzt * (K_mm1+nu) * dzmm1


endif

end function diffusion_zt_lhs

!===============================================================================
pure function diffusion_zm_lhs( K_t, K_tp1, nu,  & 
                                dztp1, dzt, dzm, level ) & 
result( lhs )

!       Description:
!       Vertical eddy diffusion of var_m:  implicit portion of the code.
!
!       The variable "var_m" stands for a variable that is located at
!       momentum grid levels.
!
!       The d(var_m)/dt equation contains an eddy diffusion term:
!
!       + d [ ( K_t + nu ) * d(var_m)/dz ] / dz.
!
!       This term is usually solved for completely implicitly, such that:
!
!       + d [ ( K_t + nu ) * d( var_m(t+1) )/dz ] / dz.
!
!       However, when a Crank-Nicholson scheme is used, the eddy
!       diffusion term has both implicit and explicit components,
!       such that:
!
!       + (1/2) * d [ ( K_t + nu ) * d( var_m(t+1) )/dz ] / dz
!          + (1/2) * d [ ( K_t + nu ) * d( var_m(t) )/dz ] / dz;
!
!       for which the implicit component is:
!
!       + (1/2) * d [ ( K_t + nu ) * d( var_m(t+1) )/dz ] / dz.
!
!       Note:  When the implicit term is brought over to the left-hand
!              side, the sign is reversed and the leading "+" in front
!              of the term is changed to a "-".
!
!       Timestep index (t) stands for the index of the current timestep,
!       while timestep index (t+1) stands for the index of the next
!       timestep, which is being advanced to in solving the d(var_m)/dt
!       equation.
!
!       The implicit portion of this term is discretized as follows:
!
!       The values of var_m are found on the momentum levels, while
!       the values of K_t are found on the thermodynamic levels.  The
!       derivatives (d/dz) of var_m are taken over the intermediate
!       thermodynamic levels.  At the intermediate thermodynamic levels,
!       d(var_m)/dz is multiplied by ( K_t + nu ).  Then, the derivative
!       of the whole mathematical expression is taken over the central
!       momentum level, which yields the desired result.
!
!       ===var_mp1=============================================== m(k+1)
!
!       -----------d(var_m)/dz--(K_tp1+nu)----------------------- t(k+1)
!
!       ===var_m====================d[(K_t+nu)*d(var_m)/dz]/dz=== m(k)
!
!       -----------d(var_m)/dz--(K_t+nu)------------------------- t(k)
!
!       ===var_mm1=============================================== m(k-1)
!
!       The vertical indices m(k+1), t(k+1), m(k), t(k), and m(k-1)
!       correspond with altitudes zm(k+1), zt(k+1), zm(k), zt(k),
!       and zm(k-1), respectively.  The letter "t" is used for
!       thermodynamic levels and the letter "m" is used for momentum
!       levels.
!
!       dzm(k)   = 1 / ( zt(k+1) - zt(k) )
!       dzt(k+1) = 1 / ( zm(k+1) - zm(k) )
!       dzt(k)   = 1 / ( zm(k) - zm(k-1) )
!
!       Note:  This function only computes the general implicit form:
!              + d [ ( K_t + nu ) * d( var_m(t+1) )/dz ] / dz.
!              For a Crank-Nicholson scheme, the left-hand side result
!              of this function will have to be multiplied by (1/2).
!              For a Crank-Nicholson scheme, the right-hand side 
!              (explicit) component needs to be computed by multiplying 
!              the left-hand side results by (1/2), reversing the sign
!              on each left-hand side element, and then multiplying 
!              each element by the appropriate var_m(t) value from the
!              appropriate vertical level.
!
!
!       Boundary Conditions:
!
!       1) Zero-flux boundary conditions.
!          This function is set up to use a zero-flux boundary condition
!          at both the lower boundary level (k=1) and the upper boundary
!          level (k=gr%nnzp).  This means that the derivative of var_m,
!          d(var_m)/dz, equals 0 at both the lower boundary and the upper
!          boundary.
!
!          In order to discretize the lower boundary condition, consider
!          a new level outside the model (momentum level 0) just below 
!          the lower boundary level (momentum level 1).  The value of 
!          var_m at the level just outside the model is defined
!          to be the same as the value of var_m at the lower boundary
!          level.  Therefore, the value of d(var_m)/dz between the level
!          just outside the model and the lower boundary level is 0,
!          satisfying the zero-flux boundary condition.  The other value
!          for d(var_m)/dz (between momentum level 2 and momentum 
!          level 1) is taken over the intermediate thermodynamic level 
!          (thermodynamic level 2), where it is multiplied by the factor
!          ( K_t(2) + nu ).  Then, the derivative of the whole expression
!          is taken over the central momentum level.
!
!          ===var_m(2)=========================================== m(2)
!
!          ------------d(var_m)/dz==(K_t(2)+nu)------------------ t(2)
!
!          ===var_m(1)===============d[(K_t+nu)*d(var_m)/dz]/dz== m(1) Boundary
!
!          ------------[d(var_m)/dz = 0]------------------------- t(1)
!
!          ===[var_m(0) = var_m(1)]====(level outside model)===== m(0)
!
!          The result is dependent only on values of K_t found at
!          thermodynamic level 2 and values of var_m found at momentum
!          levels 1 and 2.  Thus, it only affects 2 diagonals on the
!          left-hand side matrix.
!
!          The same method can be used to discretize the upper boundary
!          by considering a new level outside the model just above the
!          upper boundary level.
!
!       2) Fixed-point boundary conditions.
!          Many equations in the model use fixed-point boundary
!          conditions rather than zero-flux boundary conditions.  This
!          means that the value of var_m stays the same over the course
!          of the timestep at the lower boundary, as well as at the 
!          upper boundary.
!
!          In order to discretize the boundary conditions for equations
!          requiring fixed-point boundary conditions, either:
!          a) in the parent subroutine or function (that calls this
!             function), loop over all vertical levels from the
!             second-lowest to the second-highest, ignoring the boundary
!             levels.  Then set the values at the boundary levels in the
!             parent subroutine; or
!          b) in the parent subroutine or function, loop over all
!             vertical levels and then overwrite the results at the
!             boundary levels.
!
!          Either way, at the boundary levels, an array with a value 
!          of 1 at the main diagonal on the left-hand side and with 
!          values of 0 at all other diagonals on the left-hand side will
!          preserve the right-hand side value at that level, thus
!          satisfying the fixed-point boundary conditions.
!
!
!       Conservation Properties:
!
!       When zero-flux boundary conditions are used, this technique of
!       discretizing the eddy diffusion term leads to conservative
!       differencing.  When conservative differencing is in place, the
!       column totals for each column in the left-hand side matrix
!       (for the eddy diffusion term) should be equal to 0.  This
!       ensures that the total amount of the quantity var_m over the
!       entire vertical domain is being conserved, meaning that
!       nothing is lost due to diffusional effects.
!
!       To see that this conservation law is satisfied, compute the 
!       eddy diffusion of var_m and integrate vertically.  In 
!       discretized matrix notation (where "i" stands for the matrix 
!       column and "j" stands for the matrix row):
!
!        0 =
!        Sum_j Sum_i ( 1/dzm )_i ( dzm * ((K_t+nu)*dzt) )_ij (var_m)_j.
!
!       The left-hand side matrix, ( dzm * ((K_t+nu)*dzt) )_ij, is
!       partially written below.  The sum over i in the above equation 
!       removes dzm everywhere from the matrix below.  The sum over j
!       leaves the column totals that are desired.
!
!       Left-hand side matrix contributions from eddy diffusion term;
!       first four vertical levels:
!
!     -------------------------------------------------------------------------------------------->
!k=1 | +dzm(k)                  -dzm(k)*(K_t(k+1)+nu)*dzt(k+1)                    0
!    |   (K_t(k+1)+nu)*dzt(k+1) 
!    |          
!k=2 | -dzm(k)*                 +dzm(k)*[ (K_t(k+1)+nu)*dzt(k+1)  -dzm(k)*(K_t(k+1)+nu)*dzt(k+1)
!    |  (K_t(k)+nu)*dzt(k)       +(K_t(k)+nu)*dzt(k) ]
!    |
!k=3 |              0          -dzm(k)*(K_t(k)+nu)*dzt(k)         +dzm(k)*[ (K_t(k+1)+nu)*dzt(k+1)
!    |                                                                    +(K_t(k)+nu)*dzt(k) ]
!    |
!k=4 |              0                       0                       -dzm(k)*(K_t(k)+nu)*dzt(k)
!    |
!   \ /
!
!       Note:  The superdiagonal term from level 3 and both the main
!              diagonal and superdiagonal terms from level 4 are not
!              shown on this diagram.
!
!       Note:  The matrix shown is a tridiagonal matrix.  For a band 
!              diagonal matrix (with 5 diagonals), there would be an 
!              extra row between each of the rows shown and an extra
!              column between each of the columns shown.  However, 
!              for the purposes of the var_m eddy diffusion term, those
!              extra row and column values are all 0, and the 
!              conservation properties of the matrix aren't effected.
!              
!       If fixed-point boundary conditions are used, the matrix
!       entries at level 1 (k=1) read:  1   0   0; which means that
!       conservative differencing is not in play.  The total amount of
!       var_m over the entire vertical domain is not being conserved,
!       as amounts of var_m may be fluxed out through the upper boundary
!       or lower boundary through the effects of diffusion.
!
!       Brian Griffin.  April 26, 2008.

!       References:
!       None
!-----------------------------------------------------------------------

use grid_class, only: & 
gr       ! Variable(s)

implicit none

! Constant parameters
integer, parameter :: & 
kp1_mdiag = 1,    & ! Momentum superdiagonal index.
k_mdiag   = 2,    & ! Momentum main diagonal index.
km1_mdiag = 3    ! Momentum subdiagonal index.

! Input Variables
real, intent(in) ::  & 
K_t,    & ! Coefficient of eddy diffusivity at thermo. level (k)   [m^2/s]
K_tp1,  & ! Coefficient of eddy diffusivity at thermo. level (k+1) [m^2/s]
nu,     & ! Background constant coefficient of eddy diffusivity    [m^2/s]
dzm,    & ! Inverse of grid spacing over momentum level (k)        [1/m]
dzt,    & ! Inverse of grid spacing over thermodynamic level (k)   [1/m]
dztp1  ! Inverse of grid spacing over thermodynamic level (k+1) [1/m]

integer, intent(in) ::  & 
level ! Momentum level where calculation occurs.               [-]

! Return Variable
real, dimension(3) :: lhs

if ( level == 1 ) then

   ! k = 1; lower boundery level at surface.
   ! Only relevant if zero-flux boundary conditions are used.

   ! Momentum superdiagonal: [ x var_m(k+1,<t+1>) ]
   lhs(kp1_mdiag) = - dzm * (K_tp1+nu) * dztp1

   ! Momentum main diagonal: [ x var_m(k,<t+1>) ]
   lhs(k_mdiag)   = + dzm * (K_tp1+nu) * dztp1

   ! Momentum subdiagonal: [ x var_m(k-1,<t+1>) ]
   lhs(km1_mdiag) = 0.0


elseif ( level > 1 .and. level < gr%nnzp ) then

   ! Most of the interior model; normal conditions.

   ! Momentum superdiagonal: [ x var_m(k+1,<t+1>) ]
   lhs(kp1_mdiag) = - dzm * (K_tp1+nu) * dztp1

   ! Momentum main diagonal: [ x var_m(k,<t+1>) ]
   lhs(k_mdiag)   = + dzm * ( (K_tp1+nu)*dztp1 + (K_t+nu)*dzt )

   ! Momentum subdiagonal: [ x var_m(k-1,<t+1>) ]
   lhs(km1_mdiag) = - dzm * (K_t+nu) * dzt


elseif ( level == gr%nnzp ) then

   ! k = gr%nnzp (top level); upper boundary level.
   ! Only relevant if zero-flux boundary conditions are used.

   ! Momentum superdiagonal: [ x var_m(k+1,<t+1>) ]
   lhs(kp1_mdiag) = 0.0

   ! Momentum main diagonal: [ x var_m(k,<t+1>) ]
   lhs(k_mdiag)   = + dzm * (K_t+nu) * dzt

   ! Momentum subdiagonal: [ x var_m(k-1,<t+1>) ]
   lhs(km1_mdiag) = - dzm * (K_t+nu) * dzt


endif

end function diffusion_zm_lhs

!===============================================================================

end module diffusion
