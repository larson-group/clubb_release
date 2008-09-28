!-----------------------------------------------------------------------
! $Id$
!===============================================================================
module hyper_diffusion_4th_ord

  implicit none

  private ! Default Scope

  public :: hyper_dfsn_4th_ord_zm_lhs

contains

  !=============================================================================
  pure function hyper_dfsn_4th_ord_zm_lhs

    ! Description:
    ! Vertical 4th-order numerical diffusion of var_zm:  implicit portion of the
    ! code.
    !
    ! Fourth-order numerical diffusion, or fourth-order hyper-diffusion, is used
    ! to help eliminate small-scale noise without altering larger-scale
    ! features.
    !
    ! The variable "var_zm" stands for a variable that is located at momentum
    ! grid levels.
    !
    ! The d(var_zm)/dt equation contains a 4th-order numerical diffusion term:
    !
    ! - nu * d^4(var_zm)/dz^4.
    !
    ! This term is solved for completely implicitly, such that:
    !
    ! - nu * d^4( var_zm(t+1) )/dz^4.
    !
    ! Note:  When the term is brought over to the left-hand side, the sign 
    !        is reversed and the leading "-" in front of the term is changed 
    !        to a "+".
    !
    ! The timestep index (t+1) means that the value of var_zm being used is from
    ! the next timestep, which is being advanced to in solving the d(var_zm)/dt 
    ! equation.
    !
    ! The term is discretized as follows:
    !
    ! The five values of var_zm are found on the momentum levels.  All four
    ! derivatives (d/dz) of var_zm are taken over all the intermediate
    ! thermodynamic levels.  Then, all three derivatives (d/dz) of d(var_zm)/dz
    ! are taken over all the intermediate momentum levels, which results in the
    ! second derivatives.  Then, both derivatives (d/dz) of d^2(var_zm)/dz^2 are
    ! taken over the intermediate thermodynamic levels, which results in the
    ! third derivatives.  Finally, the derivative (d/dz) of d^3(var_zm)/dz^3 is
    ! taken over the intermediate (central) momentum level, which results in the
    ! fourth derivative.  At the central momentum level, d^4(var_zm)/dz^4 is
    ! multiplied by constant coefficient nu.
    !
    ! ==var_zmp2=============================================== m(k+2)
    !
    ! ------d(var_zm)/dz--------------------------------------- t(k+2)
    !
    ! ==var_zmp1====d^2(var_zm)/dz^2=========================== m(k+1)
    !
    ! ------d(var_zm)/dz--------d^3(var_zm)/dz^3--------------- t(k+1)
    !
    ! ==var_zm======d^2(var_zm)/dz^2========d^4(var_zm)/dz^4=== m(k)
    !
    ! ------d(var_zm)/dz--------d^3(var_zm)/dz^3--------------- t(k)
    !
    ! ==var_zmm1====d^2(var_zm)/dz^2=========================== m(k-1)
    !
    ! ------d(var_zm)/dz--------------------------------------- t(k-1)
    !
    ! ==var_zmm2=============================================== m(k-2)
    !
    ! The vertical indices m(k+2), t(k+2), m(k+1), t(k+1), m(k), t(k), m(k-1),
    ! t(k-1), and m(k-2) correspond with altitudes zm(k+2), zt(k+2), zm(k+1),
    ! zt(k+1), zm(k), zt(k), zm(k-1), zt(k-1), and zm(k-2) respectively.  The
    ! letter "t" is used for thermodynamic levels and the letter "m" is used for
    ! momentum levels.
    !
    ! dzm(k)   = 1 / ( zt(k+1) - zt(k) )
    ! dzt(k+1) = 1 / ( zm(k+1) - zm(k) )
    ! dzt(k)   = 1 / ( zm(k) - zm(k-1) )
    ! dzm(k+1) = 1 / ( zt(k+2) - zt(k+1) )
    ! dzm(k-1) = 1 / ( zt(k) - zt(k-1) )
    ! dzt(k+2) = 1 / ( zm(k+2) - zm(k+1) )
    ! dzt(k-1) = 1 / ( zm(k-1) - zm(k-2) )
    !
    ! The discretization of -nu*d^4(var_zm)/dz^4 at momentum level (k) is
    ! written out as follows:
    !
    ! -nu*dzm(k)*[ dzt(k+1)*{ dzm(k+1)*( dzt(k+2)*(var_zm(k+2)-var_zm(k+1))
    !                                   -dzt(k+1)*(var_zm(k+1)-var_zm(k)) )
    !                        -dzm(k)*( dzt(k+1)*(var_zm(k+1)-var_zm(k))
    !                                 -dzt(k)*(var_zm(k)-var_zm(k-1)) ) }
    !             -dzt(k)*{ dzm(k)*( dzt(k+1)*(var_zm(k+1)-var_zm(k))
    !                               -dzt(k)*(var_zm(k)-var_zm(k-1)) )
    !                      -dzm(k-1)*( dzt(k)*(var_zm(k)-var_zm(k-1))
    !                                 -dzt(k-1)*(var_zm(k-1)-var_zm(k-2)) ) } ].
    !
    ! Again, the term is treated completely implicitly, so the leading "-" sign
    ! changes to a "+" sign when the term is brought over to the left-hand side,
    ! and var_zm is considered to be at timestep (t+1).
    !
    !
    ! Boundary Conditions:
    !
    ! 1) Zero-flux boundary conditions.
    !    This function is set up to use zero-flux boundary conditions at both
    !    the lower boundary level and the upper boundary level.  The flux, F,
    !    is the amount of var_zm flowing normal through the boundary per unit
    !    time per unit surface area.  The derivative of the flux effects the
    !    time-tendency of var_zm, such that:
    !
    !    d(var_zm)/dt = -dF/dz.
    !
    !    For the 4th-order numerical diffusion term, -nu*d^4(var_zm)/dz^4 (which
    !    is actually -d[nu*d^3(var_zm)/dz^3]/dz with a constant coefficient, 
    !    nu), the flux is:
    !
    !    F = +nu*d^3(var_zm)/dz^3.
    !
    !    In order to have zero-flux boundary conditions, the third derivative of
    !    var_zm, d^3(var_zm)/dz^3, needs to equal 0 at both the lower boundary
    !    and the upper boundary.
    !
    !    Fourth-order numerical diffusion is used in conjunction with
    !    second-order eddy diffusion, +d[(K_zt+nu)*d(var_zm)/dz]/dz, where the
    !    coefficient of eddy diffusivity, (K_zt+nu), varies in the vertical.
    !    Both 4th-order numerical diffusion and 2nd-order eddy diffusion use the
    !    same boundary condition type at all times, which in this case is
    !    zero-flux boundary conditions.  For 2nd-order eddy diffusion, the flux
    !    is:  F = -(K_zt+nu)*d(var_zm)/dz.  In order to have zero-flux boundary
    !    conditions, the derivative of var_zm, d(var_zm)/dz, needs to equal 0 at
    !    both the lower boundary and the upper boundary.
    !
    !    Thus, the boundary conditions used for 4th-order numerical diffusion
    !    are:  d^3(var_zm)/dz^3 = 0 and d(var_zm)/dz = 0 at both the upper
    !    boundary and the lower boundary, resulting in four boundary conditions,
    !    which is the number of boundary conditions needed for a 4th-order term.
    !
    !    In order to discretize the lower boundary condition, consider a new
    !    level outside the model (momentum level 0) just below the lower
    !    boundary level (momentum level 1).  The value of var_zm at the level
    !    just outside the model is defined to be the same as the value of var_zm
    !    at the lower boundary level.  Therefore, the value of d(var_zm)/dz
    !    between the level just outside the model and the lower boundary level
    !    is 0, satisfying one of the boundary conditions.  The boundary
    !    condition d^3(var_zm)/dz^3 = 0 is also set at this level.  The rest of
    !    the levels involved are discretized normally, as listed above.
    !
    !    Since the normal discretization includes two levels on either side of
    !    the central level, the lower boundary begins to effect the
    !    discretization at momentum level 2.
    !
    !    =var_zm(4)=============================================== m(4)
    !
    !    ------d(var_zm)/dz--------------------------------------- t(4)
    !
    !    =var_zm(3)====d^2(var_zm)/dz^2=========================== m(3)
    !
    !    ------d(var_zm)/dz--------d^3(var_zm)/dz^3--------------- t(3)
    !
    !    =var_zm(2)====d^2(var_zm)/dz^2========d^4(var_zm)/dz^4=== m(2)
    !
    !    ------d(var_zm)/dz--------d^3(var_zm)/dz^3--------------- t(2)
    !
    !    =var_zm(1)====d^2(var_zm)/dz^2=========================== m(1) Boundary
    !
    !    ------[d(var_zm)/dz = 0]--------------------------------- t(1)
    !
    !    =[var_zm(0) = var_zm(1)]=====(level outside model)======= m(0)
    !
    !    The discretization of -nu*d^4(var_zm)/dz^4 at momentum level (k=2) is
    !    written out as follows:
    !
    !    -nu*dzm(k)*[ dzt(k+1)*{ dzm(k+1)*( dzt(k+2)*(var_zm(k+2)-var_zm(k+1))
    !                                      -dzt(k+1)*(var_zm(k+1)-var_zm(k)) )
    !                           -dzm(k)*( dzt(k+1)*(var_zm(k+1)-var_zm(k))
    !                                    -dzt(k)*(var_zm(k)-var_zm(k-1)) ) }
    !                -dzt(k)*{ dzm(k)*( dzt(k+1)*(var_zm(k+1)-var_zm(k))
    !                                  -dzt(k)*(var_zm(k)-var_zm(k-1)) )
    !                         -dzm(k-1)*dzt(k)*(var_zm(k)-var_zm(k-1)) } ].
    !
    !    Again, the term is treated completely implicitly, so the leading "-"
    !    sign changes to a "+" sign when the term is brought over to the
    !    left-hand side, and var_zm is considered to be at timestep (t+1).
    !
    !    The result is dependent only on values of var_zm found at momentum
    !    levels 1, 2, 3, and 4.  Thus, it only affects 4 diagonals on the
    !    left-hand side matrix.
    !
    !    The lower boundary also effects the discretization at momentum
    !    level 1.
    !
    !    =var_zm(3)=============================================== m(3)
    !
    !    ------d(var_zm)/dz--------------------------------------- t(3)
    !
    !    =var_zm(2)====d^2(var_zm)/dz^2=========================== m(2)
    !
    !    ------d(var_zm)/dz--------d^3(var_zm)/dz^3--------------- t(2)
    !
    !    =var_zm(1)====d^2(var_zm)/dz^2========d^4(var_zm)/dz^4=== m(1) Boundary
    !
    !    ------[d(var_zm)/dz = 0]--[d^3(var_zm)/dz^3 = 0]--------- t(1)
    !
    !    =[var_zm(0) = var_zm(1)]=====(level outside model)======= m(0)
    !
    !    The discretization of -nu*d^4(var_zm)/dz^4 at momentum level (k=1) is
    !    written out as follows:
    !
    !    -nu*dzm(k)*[dzt(k+1)*{ dzm(k+1)*( dzt(k+2)*(var_zm(k+2)-var_zm(k+1))
    !                                     -dzt(k+1)*(var_zm(k+1)-var_zm(k)) )
    !                          -dzm(k)*dzt(k+1)*(var_zm(k+1)-var_zm(k)) } ].
    !
    !    Again, the term is treated completely implicitly, so the leading "-"
    !    sign changes to a "+" sign when the term is brought over to the
    !    left-hand side, and var_zm is considered to be at timestep (t+1).
    !
    !    The result is dependent only on values of var_zm found at momentum
    !    levels 1, 2, and 3.  Thus, it only affects 3 diagonals on the left-hand
    !    side matrix.
    !
    !    The same method can be used to discretize the upper boundary by
    !    considering a new level outside the model just above the upper boundary
    !    level.
    !
    ! 2) Fixed-point boundary conditions.
    !    Many equations in the model use fixed-point boundary conditions rather
    !    than zero-flux boundary conditions.  This means that the value of
    !    var_zm stays the same over the course of the timestep at the lower
    !    boundary, as well as at the upper boundary.
    !
    !    For a 4th-order term, four boundary conditions are needed.  Two
    !    boundary conditions are applied at each boundary.  For the case of
    !    fixed-point boundary conditions, one of those two conditions is setting
    !    var_zm = A, where A is a constant value.  One more condition is needed.
    !    Setting the values of d(var_zm)/dz and d^3(var_zm)/dz^3 are inherently
    !    used for zero-flux (or perhaps fixed-flux) boundary conditions.
    !    Fixed-point and zero-flux boundary conditions inherently should not be
    !    invoked at the same time.  The only remaining choice for a second
    !    boundary condition for the fixed-point case is setting
    !    d^2(var_zm)/dz^2.  As it turns out, setting d^2(var_zm)/dz^2 = 0 is the
    !    appropriate condition to use because it prevents values of var_zm at
    !    levels outside the model from being involved in the discretization of
    !    -nu*d^4(var_zm)/dz^4 at momentum level 2.  Setting d^3(var_zm)/dz^3 = 0
    !    does not accomplish the same thing for the discretization of
    !    -nu*d^4(var_zm)/dz^4 at momentum level 2.  Also, as stated above,
    !    fourth-order numerical diffusion is used in conjunction with
    !    second-order eddy diffusion, +d[(K_zt+nu)*d(var_zm)/dz]/dz, where the
    !    coefficient of eddy diffusivity, (K_zt+nu), varies in the vertical.
    !    Both 4th-order numerical diffusion and 2nd-order eddy diffusion use the
    !    same boundary condition type at all times, which in this case is
    !    fixed-point boundary conditions.  For 2nd-order eddy diffusion,
    !    fixed-point boundary conditions set var_zm = A, and do not set
    !    d(var_zm)/dz.  Thus, d(var_zm)/dz cannot be set for fixed-point
    !    boundary conditions.  As previously stated, the only other boundary
    !    condition that can be invoked for a fixed-point boundary case is
    !    d^2(var_zm)/dz^2 = 0.
    !
    !    Since the normal discretization includes two levels on either side of
    !    the central level, the lower boundary begins to effect the
    !    discretization at momentum level 2.
    !
    !    =var_zm(4)=============================================== m(4)
    !
    !    ------d(var_zm)/dz--------------------------------------- t(4)
    !
    !    =var_zm(3)====d^2(var_zm)/dz^2=========================== m(3)
    !
    !    ------d(var_zm)/dz--------d^3(var_zm)/dz^3--------------- t(3)
    !
    !    =var_zm(2)====d^2(var_zm)/dz^2========d^4(var_zm)/dz^4=== m(2)
    !
    !    ------d(var_zm)/dz--------d^3(var_zm)/dz^3--------------- t(2)
    !
    !    =var_zm(1)====[d^2(var_zm)/dz^2 = 0]===================== m(1) Boundary
    !
    !    ------d(var_zm)/dz--------------------------------------- t(1)
    !
    !    =var_zm(0)===================(level outside model)======= m(0)
    !
    !    The discretization of -nu*d^4(var_zm)/dz^4 at momentum level (k=2) is
    !    written out as follows:
    !
    !    -nu*dzm(k)*[ dzt(k+1)*{ dzm(k+1)*( dzt(k+2)*(var_zm(k+2)-var_zm(k+1))
    !                                      -dzt(k+1)*(var_zm(k+1)-var_zm(k)) )
    !                           -dzm(k)*( dzt(k+1)*(var_zm(k+1)-var_zm(k))
    !                                    -dzt(k)*(var_zm(k)-var_zm(k-1)) ) }
    !                -dzt(k)*{ dzm(k)*( dzt(k+1)*(var_zm(k+1)-var_zm(k))
    !                                  -dzt(k)*(var_zm(k)-var_zm(k-1)) ) } ].
    !
    !    Again, the term is treated completely implicitly, so the leading "-"
    !    sign changes to a "+" sign when the term is brought over to the
    !    left-hand side, and var_zm is considered to be at timestep (t+1).
    !
    !    The result is dependent only on values of var_zm found at momentum
    !    levels 1, 2, 3, and 4.  Thus, it only affects 4 diagonals on the
    !    left-hand side matrix.
    !
    !    The same method can be used to discretize -nu*d^4(var_zm)/dz^4 at the
    !    second-highest momentum level (k=top-1) by setting d^2(var_zm)/dz^2 = 0
    !    at the highest momentum level.
    !
    !    The discretization at momentum level (k=1) is written to simply set the
    !    value var_zm(1) = A.  Likewise, the discretization at momentum level
    !    (k=top) is written to simply set the value var_zm(top) = B.  In order
    !    to discretize the boundary conditions at the lowest and highest
    !    vertical levels for equations requiring fixed-point boundary
    !    conditions, either:
    !    a) in the parent subroutine or function (that calls this function),
    !       loop over all vertical levels from the second-lowest to the
    !       second-highest, ignoring the lowest and highest levels.  Then set
    !       the values at the lowest and highest levels in the parent
    !       subroutine; or
    !    b) in the parent subroutine or function, loop over all vertical levels
    !       and then overwrite the results at the lowest and highest levels.
    !
    !    Either way, at the lowest and highest levels, an array with a value
    !    of 1 at the main diagonal on the left-hand side and with values of 0 at
    !    all other diagonals on the left-hand side will preserve the right-hand
    !    side value at that level, thus satisfying the fixed-point boundary
    !    conditions.
    !
    !
    ! Conservation Properties:
    !
    ! When zero-flux boundary conditions are used, this technique of
    ! discretizing the 4th-order numerical diffusion term leads to conservative
    ! differencing.  When conservative differencing is in place, the column
    ! totals for each column in the left-hand side matrix (for the 4th-order
    ! numerical diffusion term) should be equal to 0.  This ensures that the
    ! total amount of the quantity var_zm over the entire vertical domain is
    ! being conserved, meaning that nothing is lost due to diffusional effects.
    !
    ! To see that this conservation law is satisfied, compute the 4th-order
    ! numerical diffusion of var_zm and integrate vertically.  In discretized
    ! matrix notation (where "i" stands for the matrix column and "j" stands for
    ! the matrix row):
    !
    !  0 = Sum_j Sum_i ( 1/dzm )_i ( nu*dzm*dzt*dzm*dzt )_ij (var_zm)_j.
    !
    ! The left-hand side matrix, ( nu*dzm*dzt*dzm*dzt )_ij, is partially written
    ! below.  The sum over i in the above equation removes the first dzm(k)
    ! everywhere from the matrix below.  The sum over j leaves the column totals
    ! that are desired.
    !
    ! Left-hand side matrix contributions from 4th-order numerical diffusion
    ! (or hyper-diffusion) term; first five vertical levels:
    !
    !         column 1    ||    column 2    ||    column 3    ||    column 4    ||    column 5
    !    ------------------------------------------------------------------------------------------>
    !   | +nu             -nu               +nu
    !   | *dzm(k)         *dzm(k)           *dzm(k)
    !   |  *[ dzt(k+1)     *[ dzt(k+1)       *dzt(k+1)
    !   |     *{ dzm(k+1)     *{ dzm(k+1)     *dzm(k+1)
    !k=1|        *dzt(k+1)       *( dzt(k+2)   *dzt(k+2)               0                  0
    !   |       +dzm(k)            +dzt(k+1) )
    !   |        *dzt(k+1) }    +dzm(k)
    !   |   ]                    *dzt(k+1) } ]
    !   |
    !   | -nu             +nu               -nu               +nu
    !   | *dzm(k)         *dzm(k)           *dzm(k)           *dzm(k)
    !   |  *[ dzt(k+1)     *[ dzt(k+1)       *[ dzt(k+1)       *dzt(k+1)
    !   |     *dzm(k)         *{ dzm(k+1)       *{ dzm(k+1)     *dzm(k+1)
    !   |      *dzt(k)           *dzt(k+1)         *( dzt(k+2)   *dzt(k+2)
    !   |    +dzt(k)            +dzm(k)              +dzt(k+1) )
    !k=2|     *{ dzm(k)          *( dzt(k+1)      +dzm(k)                                 0
    !   |        *dzt(k)           +dzt(k) ) }     *dzt(k+1) }
    !   |       +dzm(k-1)    +dzt(k)           +dzt(k)
    !   |        *dzt(k) } ]  *{ dzm(k)         *dzm(k)
    !   |                        *( dzt(k+1)     *dzt(k+1) ]
    !   |                          +dzt(k) )
    !   |                       +dzm(k-1)
    !   |                        *dzt(k) } ]
    !   |
    !   | +nu             -nu               +nu               -nu               +nu
    !   | *dzm(k)         *dzm(k)           *dzm(k)           *dzm(k)           *dzm(k)
    !   |  *dzt(k)         *[ dzt(k+1)       *[ dzt(k+1)       *[ dzt(k+1)       *dzt(k+1)
    !   |   *dzm(k-1)         *dzm(k)           *{ dzm(k+1)       *{ dzm(k+1)     *dzm(k+1)
    !   |    *dzt(k-1)         *dzt(k)             *dzt(k+1)         *( dzt(k+2)   *dzt(k+2)
    !   |                    +dzt(k)              +dzm(k)              +dzt(k+1) )
    !k=3|                     *{ dzm(k)            *( dzt(k+1)      +dzm(k)
    !   |                        *dzt(k)             +dzt(k) ) }     *dzt(k+1) }
    !   |                       +dzm(k-1)      +dzt(k)           +dzt(k)
    !   |                        *( dzt(k)      *{ dzm(k)         *dzm(k)
    !   |                          +dzt(k-1) )     *( dzt(k+1)     *dzt(k+1) ]
    !   |                      } ]                   +dzt(k) )
    !   |                                         +dzm(k-1)
    !   |                                          *dzt(k) } ]
    !   |
    !   |                 +nu               -nu               +nu               -nu
    !   |                 *dzm(k)           *dzm(k)           *dzm(k)           *dzm(k)
    !   |                  *dzt(k)           *[ dzt(k+1)       *[ dzt(k+1)       *[ dzt(k+1)
    !   |                   *dzm(k-1)           *dzm(k)           *{ dzm(k+1)       *{ dzm(k+1)
    !   |                    *dzt(k-1)           *dzt(k)             *dzt(k+1)         *( dzt(k+2)
    !   |                                      +dzt(k)              +dzm(k)              +dzt(k+1) )
    !k=4|        0                              *{ dzm(k)            *( dzt(k+1)      +dzm(k)
    !   |                                          *dzt(k)             +dzt(k) ) }     *dzt(k+1) }
    !   |                                         +dzm(k-1)      +dzt(k)           +dzt(k)
    !   |                                          *( dzt(k)      *{ dzm(k)         *dzm(k)
    !   |                                            +dzt(k-1) )     *( dzt(k+1)     *dzt(k+1) ]
    !   |                                        } ]                   +dzt(k) )
    !   |                                                           +dzm(k-1)
    !   |                                                            *dzt(k) } ]
    !   |
    !   |                                   +nu               -nu               +nu
    !   |                                   *dzm(k)           *dzm(k)           *dzm(k)
    !   |                                    *dzt(k)           *[ dzt(k+1)       *[ dzt(k+1)
    !   |                                     *dzm(k-1)           *dzm(k)           *{ dzm(k+1)
    !   |                                      *dzt(k-1)           *dzt(k)             *dzt(k+1)
    !   |                                                        +dzt(k)              +dzm(k)
    !k=5|        0                  0                             *{ dzm(k)            *( dzt(k+1)
    !   |                                                            *dzt(k)             +dzt(k) ) }
    !   |                                                           +dzm(k-1)      +dzt(k)
    !   |                                                            *( dzt(k)      *{ dzm(k)
    !   |                                                              +dzt(k-1) )     *( dzt(k+1)
    !   |                                                          } ]                   +dzt(k) )
    !   |                                                                             +dzm(k-1)
    !   |                                                                              *dzt(k) } ]
    !  \ /
    !
    ! Note:  The super-super diagonal term from level 4 and both the super
    !        diagonal and super-super diagonal terms from level 5 are not shown
    !        on this diagram.
    !
    ! Note:  The matrix shown is a five-diagonal matrix.  For a nine-diagonal
    !        matrix, there would be an extra row between each of the rows shown
    !        and an extra column between each of the columns shown.  However,
    !        for the purposes of the var_zm 4th-order hyper-diffusion term,
    !        those extra row and column values are all 0, and the conservation
    !        properties of the matrix aren't effected.
    !
    ! For the case of fixed-point boundary conditions, the contributions of the
    ! 4th-order hyper-diffusion term are as follows (only the top 2 levels
    ! differ from the matrix diagram above):
    !
    !         column 1    ||    column 2    ||    column 3    ||    column 4    ||    column 5
    !    ------------------------------------------------------------------------------------------>
    !k=1|        0                 0                 0                 0                  0
    !   |
    !   | -nu             +nu               -nu               +nu
    !   | *dzm(k)         *dzm(k)           *dzm(k)           *dzm(k)
    !   |  *[ dzt(k+1)     *[ dzt(k+1)       *[ dzt(k+1)       *dzt(k+1)
    !   |     *dzm(k)         *{ dzm(k+1)       *{ dzm(k+1)     *dzm(k+1)
    !   |      *dzt(k)           *dzt(k+1)         *( dzt(k+2)   *dzt(k+2)
    !k=2|    +dzt(k)            +dzm(k)              +dzt(k+1) )                          0
    !   |     *dzm(k)            *( dzt(k+1)      +dzm(k)
    !   |      *dzt(k) ]           +dzt(k) ) }     *dzt(k+1) }
    !   |                    +dzt(k)           +dzt(k)
    !   |                     *dzm(k)           *dzm(k)
    !   |                      *( dzt(k+1)       *dzt(k+1) ]
    !   |                        +dzt(k) ) ]
    !  \ /
    !
    ! For the left-hand side matrix as a whole, the matrix entries at level 1
    ! (k=1) read:  1  0  0  0  0.  For the case of fixed-point boundary
    ! conditions, conservative differencing is not in play.  The total amount of
    ! var_zm over the entire vertical domain is not being conserved, as amounts
    ! of var_zm may be fluxed out through the upper boundary or lower boundary
    ! through the effects of diffusion.
    !
    ! Brian Griffin.  September 28, 2008.

    ! References:
    ! None
    !-----------------------------------------------------------------------

    implicit none

    ! Constant parameters

    ! Input Variables

    ! Return Variable



    return

  end function hyper_dfsn_4th_ord_zm_lhs

!===============================================================================

end module hyper_diffusion_4th_ord
