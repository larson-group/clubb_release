! $Id$
!===============================================================================
module new_pdf

  ! Description:
  ! The portion of CLUBB's multivariate, two-component PDF that is the
  ! trivariate, two-component normal PDF of vertical velocity (w), total water
  ! mixing ratio (rt), and liquid water potential temperature (thl).

  ! References:
  ! Griffin and Larson (2018)
  !-------------------------------------------------------------------------

  implicit none

  contains

  !=============================================================================
  subroutine new_pdf_driver

    ! Description:
    ! Selects which variable is used to set the mixture fraction for the PDF
    ! ("the setter") and which variables are handled after that mixture fraction
    ! has been set ("the responders").  Traditionally, w has been used to set
    ! the PDF.  However, here, the variable with the greatest magnitude of
    ! skewness is used to set the PDF.

    ! References:
    !-----------------------------------------------------------------------

    implicit none

    ! Input Variable(s)

    ! Output Variable(s)

    ! Local Variable(s)





    return

  end subroutine new_pdf_driver

  !=============================================================================
  !
  ! DESCRIPTION OF METHOD FOR THE VARIABLE THAT SETS THE MIXTURE FRACTION
  ! =====================================================================
  !
  ! In order to find equations for the five PDF parameters for the setting
  ! variable:  mu_w_1, mu_w_2, sigma_w_1, sigma_w_2, and mixt_frac, five
  ! equations are needed.  These five equations are the equations for <w>,
  ! <w'^2>, and <w'^3> as found by integrating over the PDF.  Additionally,
  ! two more equations, which involve tunable parameters F_w and zeta_w, and
  ! which are used to help control the spread of the PDF component means and
  ! the size of the PDF component variances compared to each other, are used
  ! in this equation set.  The five equations are:
  !
  ! <w> = mixt_frac * mu_w_1 + ( 1 - mixt_frac ) * mu_w_2;
  !
  ! <w'^2> = mixt_frac * ( ( mu_w_1 - <w> )^2 + sigma_w_1^2 )
  !          + ( 1 - mixt_frac ) * ( ( mu_w_2 - <w> )^2 + sigma_w_2^2 );
  !
  ! <w'^3> = mixt_frac * ( mu_w_1 - <w> )
  !                    * ( ( mu_w_1 - <w> )^2 + 3 * sigma_w_1^2 )
  !          + ( 1 - mixt_frac ) * ( mu_w_2 - <w> )
  !                              * ( ( mu_w_2 - <w> )^2 + 3 * sigma_w_2^2 );
  !
  ! mu_w_1 - <w> = sqrt(F_w) * ( sqrt( 1 - mixt_frac ) / sqrt( mixt_frac ) )
  !                * sqrt( <w'^2> );
  !
  ! where 0 <= F_w <= 1; and
  !
  ! 1 + zeta_w = ( mixt_frac * sigma_w_1^2 )
  !              / ( ( 1 - mixt_frac ) * sigma_w_2^2 );
  !
  ! where zeta_w > -1.

  !=============================================================================

end module new_pdf
