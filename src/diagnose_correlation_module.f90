! $Id$
module diagnose_correlations_module 

  implicit none 

  public :: diagnose_KK_corr, &
            diagnose_LH_corr

  private :: diagnose_corr

  contains 



!-----------------------------------------------------------------------
  subroutine diagnose_KK_corr( Ncm, rrainm, Nrm, & ! intent(in)
                               Ncp2_on_Ncm2, rrp2_on_rrm2, Nrp2_on_Nrm2, &
                               corr_sw, corr_rrw, corr_Nrw, corr_Ncw, &
                               corr_rrNr, corr_srr, corr_sNr, corr_sNc ) ! intent(inout)

! Description:
!   This subroutine diagnoses the correlation matrix in order to feed it 
!   into KK microphysics.   

! References:
!   Larson et al. (2011), J. of Geophysical Research, Vol. 116, D00T02
!   (see CLUBB-Trac:ticket:514)
!-----------------------------------------------------------------------


		use clubb_precision, only: &
		    core_rknd ! Variable(s)

	 	use grid_class, only: &
				gr ! Variable(s)

		! Local Constants
		integer, parameter :: &
			n_variables = 5

		! Input Variables
		real( kind = core_rknd ), intent(in) :: &
			Ncm, &
			rrainm, &
			Nrm, &
			Ncp2_on_Ncm2, &
			rrp2_on_rrm2, &
			Nrp2_on_Nrm2, &
      corr_sw, &
      corr_rrw, &
      corr_Nrw, &
      corr_Ncw 	


		! Input/Output Variables
		real( kind = core_rknd ), intent(inout) :: &
			corr_rrNr, &
			corr_srr,  &
			corr_sNr,	 &
			corr_sNc


		! Local Variables
		real( kind = core_rknd ), dimension(n_variables, n_variables) :: &
			corr_matrix_approx

		real( kind = core_rknd ) :: &
			s_mellorp2_on_s_mellorm2

		integer :: &
			ii_w = 1,		   &
			ii_s = 2,			 &
			ii_rrain = 3,  &
			ii_Nr = 4,		 &
			ii_Nc = 5	

		real( kind = core_rknd ), dimension(n_variables) :: & 
		  xp2_on_xm2 ! ratios of x_variance over x_mean^2

		integer :: i, j ! Loop Iterators
		

		!-------------------- Begin code --------------------

    ! S_i is set to 1 for s_mellor, because s_mellorm could be 0
    s_mellorp2_on_s_mellorm2 = 1

		! set up xp2_on_xm2
		xp2_on_xm2(ii_w) = 1
		xp2_on_xm2(ii_s) = s_mellorp2_on_s_mellorm2
		xp2_on_xm2(ii_rrain) = rrp2_on_rrm2
		xp2_on_xm2(ii_Nr) = Nrp2_on_Nrm2
		xp2_on_xm2(ii_Nc) = Ncp2_on_Ncm2
	

		! initialize the correlation matrix with 0
		do i=2, n_variables
		  do j=1, n_variables	
				corr_matrix_approx(i,j) = 0	
			end do
		end do

		! set diagonal of the corraltion matrix to 1
		do i = 1, n_variables
			corr_matrix_approx(i,i) = 1
		end do
		
		! set the first row to the corresponding prescribed correlations
		corr_matrix_approx(1,ii_s) = corr_sw
		corr_matrix_approx(1,ii_rrain) = corr_rrw
		corr_matrix_approx(1,ii_Nr) = corr_Nrw
		corr_matrix_approx(1,ii_Nc) = corr_Ncw

    call diagnose_corr( n_variables, xp2_on_xm2, & !intent(in)
                            corr_matrix_approx ) ! intent(inout)    

    corr_rrNr = corr_matrix_approx(ii_rrain, ii_Nr)
    corr_srr = corr_matrix_approx(ii_s, ii_rrain)
    corr_sNr = corr_matrix_approx(ii_s, ii_Nr)
    corr_sNc = corr_matrix_approx(ii_s, ii_Nc)
    
	end subroutine diagnose_KK_corr

!-----------------------------------------------------------------------
  subroutine diagnose_LH_corr( xp2_on_xm2, & !intent(in)
                            corr_matrix_approx ) ! intent(inout)
! Description:
!   This subroutine diagnoses the correlation matrix in order to feed it 
!   into SILHS microphysics.   

! References:
!   Larson et al. (2011), J. of Geophysical Research, Vol. 116, D00T02
!   (see CLUBB Trac ticket#514)
!-----------------------------------------------------------------------

    use clubb_precision, only: &
		    core_rknd ! Variable(s)

		! Local Constants
		integer, parameter :: &
			n_variables = 5

    ! Input/Output variables
    real( kind = core_rknd ), dimension(n_variables, n_variables), intent(inout) :: &
			corr_matrix_approx

    ! Local Variables
		real( kind = core_rknd ), dimension(n_variables) :: & 
		  xp2_on_xm2 ! ratios of x_variance over x_mean^2

  

	end subroutine diagnose_LH_corr



!-----------------------------------------------------------------------
  subroutine diagnose_corr( n_variables, xp2_on_xm2, & !intent(in)
                            corr_matrix_approx ) ! intent(inout)
! Description:
!   This subroutine diagnoses the correlation matrix for each timestep.   

! References:
!   Larson et al. (2011), J. of Geophysical Research, Vol. 116, D00T02
!   (see CLUBB Trac ticket#514)
!-----------------------------------------------------------------------

		use error_code, only:  & 
			clubb_at_least_debug_level ! Function

		use clubb_precision, only: &
      core_rknd ! Variable(s)

		use parameters_tunable, only:  & 
			alpha_corr ! Constant(s)

    intrinsic :: &
      sqrt, abs, sign

    ! Input Variables
    integer :: &
      n_variables  ! number of variables in the correlation matrix
    
    real( kind = core_rknd ), dimension(n_variables), intent(in) :: & 
      xp2_on_xm2 ! ratios of x_variance over x_mean^2

    ! Input/Output Variables
    real( kind = core_rknd ), dimension(n_variables,n_variables), intent(inout) :: &
      corr_matrix_approx ! correlation matrix


    ! Local Variables
    integer :: i, j ! Loop iterator
    real( kind = core_rknd ) :: f_ij
		real( kind = core_rknd ), dimension(n_variables) :: &
			s_1j, & ! s_1j = sqrt(1-c_1j^2)
			sqrt_xp2_on_xm2 ! = sqrt(xp2_on_xm2)

	
	  !-------------------- Begin code --------------------


		! calculate all square roots
		do i = 1, n_variables
			s_1j(i) = sqrt(1-corr_matrix_approx(1,i)**2)
		end do
		
		do i = 1, n_variables
			sqrt_xp2_on_xm2(i) = sqrt(xp2_on_xm2(i))
		end do

    ! Diagnose the missing correlations (upper triangle)
    do i = 2, (n_variables-1)
      do j = (i+1), n_variables
        
        !print *, "-- i = ", i, ", j = ", j
        !print *, "---- c_1", i, " = ", corr_matrix_approx(1,i)
        !print *, "---- c_1", j, " = ", corr_matrix_approx(1,j)
        !print *, "---- s_1", i, " = ", s_1j(i)
        !print *, "---- s_1", j, " = ", s_1j(j)

				! formular (16) in the ref. paper 
        f_ij = alpha_corr * sqrt_xp2_on_xm2(i) * sqrt_xp2_on_xm2(j) &
        * sign(1.0_core_rknd,corr_matrix_approx(1,i)*corr_matrix_approx(1,j)) 

        !print *, "---- S_", i, " = ", sqrt_xp2_on_xm2(i)
        !print *, "---- S_", j, " = ", sqrt_xp2_on_xm2(j)
        !print *, "---- f_", i, j, " = ", f_ij

				! make sure -0.99<=f_ij<=0.99
				if ( clubb_at_least_debug_level( 2 ) ) then		
					if ( f_ij < -0.99 ) then
						f_ij = -0.99
					else if ( f_ij > 0.99 ) then
						f_ij = 0.99
					end if
				end if

        ! formular (15) in the ref. paper
        corr_matrix_approx(i,j) = corr_matrix_approx(1,i) * corr_matrix_approx(1,j) &
				+ f_ij * s_1j(i) * s_1j(j)

        !print *, "---- c_", i, j, " = ", corr_matrix_approx(i,j)

      end do ! do j
    end do ! do i
    
  end subroutine diagnose_corr 

end module diagnose_correlations_module 
