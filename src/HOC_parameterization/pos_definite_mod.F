        module pos_definite_mod

        implicit none

        public :: pos_definite_adj, band_mult

        private

        contains
!-----------------------------------------------------------------------
        subroutine pos_definite_adj
     .            ( dt, flux_lhs, solution, field_n, field_pd, flux_pd )
!       Description:
!       References:
!       Comment on ``A Positive Definite Advection Scheme Obtained by
!       Nonlinear Renormalization of the Advective Fluxes''
!-----------------------------------------------------------------------
        use grid_class
        use constants, only : eps
        use stats_prec, only: time_prec

        implicit none

        ! Constant parameters
        real, parameter :: beta = 1.0 ! Coefficient of 'Y' (i.e. field_np1)

        intrinsic :: eoshift, kind, any, min, max

        real(kind=time_prec), intent(in) ::
     .  dt ! Timestep    [s]

        real, dimension(3,gr%nnzp), intent(inout) :: 
     .  flux_lhs ! The matrix 'A' in sgbmv

        real, dimension(gr%nnzp), intent(in) :: 
     .  field_n ! The field (e.g. rtm) at n, prior to n+1

        real, dimension(gr%nnzp), intent(out) :: 
     .  flux_pd, ! Budget of the change in the flux term due to the scheme
     .  field_pd ! Budget of the change in the mean term due to the scheme

        ! Output Variables
        real, target, dimension(2*gr%nnzp), intent(inout) :: solution

        ! Local Variables

        real :: alpha ! Coefficient for matrix 'A'

        integer :: 
     .  ndim,  ! Dimension of solution vector
     .  k, km1, kp1 ! Loop indices

        real, pointer, dimension(:) ::
     .  field_np1,  ! Field at n+1 (e.g. rtm in [kg/kg])
     .  flux_np1    ! Flux applied to field

        real, dimension(gr%nnzp) ::
     .  flux_plus, flux_minus,! [F_i+1/2]^+ [F_i+1/2]^- in Smolarkiewicz 
     .  fout,                 ! (A4) F_i{}^OUT, or the sum flux_plus+flux_minus
     .  flux_lim,             ! Correction applied to flux at n+1
     .  field_nonlim          ! Temporary variable for calculation

!-----------------------------------------------------------------------

        ndim = 2*gr%nnzp

        ! Odd indices contain the mean field & even indices contain the flux
        field_np1 => solution(1:ndim:2)
        flux_np1  => solution(2:ndim:2)

        ! If all the values are positive or the values at the previous 
        ! timestep were negative, then just return
        if ( .not. any( field_np1 < 0. ) .or. any( field_n < 0. ) ) then
          flux_pd  = 0.
          field_pd = 0.
          return
        end if

        !print *, "Correcting flux"

        do k = 1, gr%nnzp, 1
          km1 = max( k-1, 1 )
          kp1 = min( k+1, gr%nnzp )
        
          flux_plus(k)  =  max( 0., flux_np1(k) )
          flux_minus(k) = -min( 0., flux_np1(k) )

          ! We place a limiter of eps to prevent a divide by zero
          fout(k) = max( flux_plus(k) + flux_minus(km1), eps )
        end do 

        do k = 2, gr%nnzp, 1
          kp1 = min( k+1, gr%nnzp )

          ! Eqn 10 from Smolarkowicz (1989)

          flux_lim(k)
     .    = max( min( flux_np1(k),
     .                ( flux_plus(k)/fout(k) ) * field_n(k)
     .              ),
     .           -( ( flux_minus(k)/fout(kp1) ) * field_n(kp1) )
     .         )
        end do

        flux_lim(1) = flux_np1(1)

        flux_pd = flux_lim - flux_np1

        field_nonlim = field_np1

        alpha = -dt
        ! Appy change to field at n+1
        call band_mult( 'N', gr%nnzp, gr%nnzp, 1, 1, 1, 1,
     .                  alpha, beta, flux_lhs, flux_pd, 
     .                  field_np1 ) 

        ! Total change in field due to this calculation
        field_pd = field_nonlim - field_np1

        ! Replace the non-limited flux with the limited flux
        flux_np1 = flux_lim

        return
        end subroutine pos_definite_adj
!-----------------------------------------------------------------------

        subroutine band_mult( trans, ndim, mdim, nsup, nsub, yinc, xinc,
     .                        alpha, beta, lhs, xvec, yvec )
!       Description:
!       Wrapper subroutine for banded matrix by vector multiplication in
!       the level 2 BLAS library.

!       References:
!       <http://www.netlib.org/blas/>
!-----------------------------------------------------------------------

        implicit none

        ! External
        ! Level 2 BLAS to multiply a vector by a band diagonal matrix
        external :: sgbmv, dgbmv

        character(len=1) :: 
     .  trans  ! Whether to use the transposition of the lhs matrix

        integer, intent(in) ::
     .  ndim, mdim,  ! Dimensions of the matrix when not compact
     .  nsup, nsub,  ! Super and Sub diagonals
     .  yinc, xinc   ! Increments of y and x vector

        real, intent(in) ::
     .  alpha,       ! Coefficient of the matrix lhs
     .  beta         ! Coefficient of Y vector

        real, dimension(nsup+nsub+1, ndim), intent(in) ::
     .  lhs          ! The matrix 'A' in the blas subroutine

        real, dimension(ndim), intent(in) ::
     .  xvec ! The vector X

        real, dimension(ndim), intent(inout) ::
     .  yvec ! The vector Y


!-----------------------------------------------------------------------
!       *** BLAS 2 routine ***
!       SUBROUTINE DGBMV(TRANS,M,N,KL,KU,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
!-----------------------------------------------------------------------
        ! Multiply so that Y := alpha*A*X + beta*Y

        if ( kind( lhs(1,1) ) == 4 ) then
          call sgbmv( trans, ndim, mdim, nsub, nsup, 
     .                alpha, lhs, nsup+nsub+1, 
     .                xvec, xinc, beta, yvec, yinc )

        else if ( kind( lhs(1,1) ) == 8 ) then
          call dgbmv( trans, ndim, mdim, nsub, nsup, 
     .                alpha, lhs, nsup+nsub+1, 
     .                xvec, xinc, beta, yvec, yinc )
        end if

        return
        end subroutine band_mult
!-----------------------------------------------------------------------

        end module pos_definite_mod
