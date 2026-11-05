"""JAX implementation of ADG1 PDF driver (adg1_adg2_3d_luhar_pdf.F90).

Contains code to calculate the PDF parameters using the Analytic Double
Gaussian 1 (ADG1) PDF closure, the Analytic Double Gaussian 2 (ADG2) PDF
closure, or the 3D Luhar PDF closure.  The ADG1 PDF is the original PDF used
by CLUBB in interactive runs.

Porting deviations:
- Fortran output and inout arguments are returned as tuples or dictionaries.
- The Fortran driver routines are split into smaller pure functions so each
  closure component is directly testable and differentiable.
- ``err_info`` debug/fatal mutation is not threaded through these low-level
  pure functions.  Callers are expected to enforce tolerance guards before
  entering the jitted closure path.
- ``calc_comp_corrs_binormal`` lives in its Fortran home
  ``pdf_utilities.F90`` and is imported/re-exported here because the Fortran
  ADG driver obtains it via ``use pdf_utilities``.
- Comments tied to commented-out or inactive Fortran branches are omitted when
  the corresponding JAX code is absent.

Implements:
  ADG1_w_closure             - mixture fraction and w PDF component params
  ADG1_ADG2_responder_params - PDF component params for rt/thl/u/v/sclr
  ADG1_pdf_driver            - top-level ADG1 PDF parameter driver
  calc_comp_corrs_binormal   - PDF component correlations (pdf_utilities.F90;
                                   sqrt_clipped-hardened grad variant of pdf_utilities.calc_comp_corrs_binormal)
  calc_Luhar_params / close_Luhar_pdf / max_cubic_root / backsolve_Luhar_params
  Luhar_3D_pdf_driver / ADG2_pdf_driver

The PDF moment-integral closures (calc_{wp2xp,wpxp2,wp2xp2,wp4,wpxpyp}_pdf and calc_w_up_in_cloud)
live in pdf_closure_module.py, mirroring their Fortran home pdf_closure_module.F90.
"""

import jax
import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.advance_helper_module import sqrt_clipped

from clubb_jax.src.CLUBB_core.constants_clubb import rt_tol, thl_tol, w_tol_sqd, zero_threshold, eps



def ADG1_w_closure(wm, wp2, Skw, sigma_sqd_w, sqrt_wp2, mixt_frac_max_mag):
    """Calculates the mixture fraction, the PDF component means of w, and the PDF
    component variances of w for the Analytic Double Gaussian 1 (ADG1)
    closure.  It sets the widths of both w Gaussians to be the same, and
    furthermore, relates them to the overall variance of w (<w'^2>) by a
    parameter, sigma_sqd_w.  The equation is:

    sigma_w_1^2 = sigma_w_2^2 = sigma_sqd_w * <w'^2>;

    where sigma_w_1^2 is the variance of w in the 1st PDF component,
    sigma_w_2^2 is the variance of w in the 2nd PDF component, and parameter
    sigma_sqd_w must have a value within the range 0 <= sigma_sqd_w <= 1.

    References:
    Golaz, J-C., V. E. Larson, and W. R. Cotton, 2002a: A PDF-based model for
    boundary layer clouds. Part I: Method and model description. J. Atmos.
    Sci., 59, 3540–3551.

    Vincent E. Larson and Jean-Christophe Golaz, 2005: Using Probability
    Density Functions to Derive Consistent Closure Relationships among
    Higher-Order Moments. Mon. Wea. Rev., 133, 1023–1042.

    Returns:
        (w_1, w_2, w_1_n, w_2_n, varnce_w_1, varnce_w_2, mixt_frac)

    Fortran declarations (ADG1_w_closure in src/CLUBB_core/adg1_adg2_3d_luhar_pdf.F90):

        !--------------------------- Input Variables ---------------------------
        integer, intent(in) :: &
          ngrdcol,  & ! Number of grid columns
          nz          ! Number of vertical level

        real( kind = core_rknd ), dimension(ngrdcol,nz),  intent(in) :: &
          wm,                & ! Mean of w (overall)                   [m/s]
          wp2,               & ! Variance of w (overall)               [m^2/s^2]
          Skw,               & ! Skewness of w                         [-]
          sigma_sqd_w,       & ! Widths of each w Gaussian             [-]
          sqrt_wp2             ! Square root of the variance of w      [m/s]

        real( kind = core_rknd ), dimension(ngrdcol),  intent(in) :: &
          mixt_frac_max_mag ! Maximum allowable mag. of mixt_frac   [-]

        !--------------------------- Input/Output Variables ---------------------------
        type(err_info_type), intent(inout) :: &
          err_info

        !--------------------------- Output Variables ---------------------------
        real( kind = core_rknd ), dimension(ngrdcol,nz), intent(out) :: &
          w_1,        & ! Mean of w (1st PDF component)                [m/s]
          w_2,        & ! Mean of w (2nd PDF component)                [m/s]
          w_1_n,      & ! Normalized mean of w (1st PDF component)     [-]
          w_2_n,      & ! Normalized mean of w (2nd PDF component)     [-]
          varnce_w_1, & ! Variance of w (1st PDF component)            [m^2/s^2]
          varnce_w_2, & ! Variance of w (2nd PDF component)            [m^2/s^2]
          mixt_frac     ! Mixture fraction                             [-]
    """
    # TODO: Fortran's debug tolerance error/header depends on err_info and its
    # fatal-return contract, absent from this pure closure; skipped in print-only work.
    _SKW_TOL = 1.0e-5

    # Width (standard deviation) parameters are non-zero

    # The variable "mixt_frac" is the weight of the 1st PDF component.  The
    # weight of the 2nd PDF component is "1-mixt_frac".  If there isn't any
    # skewness of w (Sk_w = 0 because w'^3 = 0), mixt_frac = 0.5, and both
    # PDF components are equally weighted.  If there is positive skewness of
    # w (Sk_w > 0 because w'^3 > 0), 0 < mixt_frac < 0.5, and the 2nd PDF
    # component has greater weight than does the 1st PDF component.  If there
    # is negative skewness of w (Sk_w < 0 because w'^3 < 0),
    # 0.5 < mixt_frac < 1, and the 1st PDF component has greater weight than
    # does the 2nd PDF component.
    denom_sq = 4.0 * (1.0 - sigma_sqd_w) ** 3 + Skw ** 2
    mf_formula = 0.5 * (1.0 - Skw / jnp.sqrt(denom_sq))
    mixt_frac = jnp.where(jnp.abs(Skw) <= _SKW_TOL, 0.5, mf_formula)

    # Clip mixt_frac, and 1 - mixt_frac, to avoid dividing by a small number.
    # Formula for mixt_frac_max_mag =
    # 1 - ( 1/2 * ( 1 - Skw_max
    #                   / sqrt( 4*( 1 - sigma_sqd_w )^3 + Skw_max^2 ) ) ),
    # where sigma_sqd_w is fixed at 0.4 for this calculation.
    limit = jnp.asarray(mixt_frac_max_mag)[..., None]
    mixt_frac = jnp.clip(mixt_frac, 1.0 - limit, limit)

    # The normalized mean of w for Gaussian "plume" 1 is w_1_n.  It's value
    # will always be greater than 0.  As an example, a value of 1.0 would
    # indicate that the actual mean of w for Gaussian "plume" 1 is found 1.0
    # standard deviation above the overall mean for w.
    one_minus_mf = 1.0 - mixt_frac
    sigma_factor = 1.0 - sigma_sqd_w
    # sqrt_clipped (REFACTOR B5): sigma_factor = 1-sigma_sqd_w -> 0 in well-mixed/surface layers, so the bare
    # sqrt has an inf reverse-mode gradient there (the bomex thlm surface nan). Forward-identical (arg >= 0).
    w_1_n = sqrt_clipped(one_minus_mf / mixt_frac * sigma_factor)
    # The normalized mean of w for Gaussian "plume" 2 is w_2_n.  It's value
    # will always be less than 0.  As an example, a value of -0.5 would
    # indicate that the actual mean of w for Gaussian "plume" 2 is found 0.5
    # standard deviations below the overall mean for w.
    w_2_n = -sqrt_clipped(mixt_frac / one_minus_mf * sigma_factor)

    # The mean of w for Gaussian "plume" 1 is w_1.
    w_1 = wm + sqrt_wp2 * w_1_n
    # The mean of w for Gaussian "plume" 2 is w_2.
    w_2 = wm + sqrt_wp2 * w_2_n

    # The variance of w for Gaussian "plume" 1 for varnce_w_1.
    varnce_w_1 = sigma_sqd_w * wp2
    # The variance of w for Gaussian "plume" 2 for varnce_w_2.
    # The variance in both Gaussian "plumes" is defined to be the same.
    varnce_w_2 = sigma_sqd_w * wp2

    return w_1, w_2, w_1_n, w_2_n, varnce_w_1, varnce_w_2, mixt_frac


def ADG1_ADG2_responder_params(xm, xp2, wp2, sqrt_wp2, wpxp,
                                   w_1_n, w_2_n, mixt_frac,
                                   sigma_sqd_w, beta):
    """Calculates the PDF component means and PDF component variances for thl,
    rt, or sclr when either the ADG1 PDF or ADG2 PDF are used.

    The normalized variance for thl, rt, and sclr for "plume" 1 is:

    { 1 - [1/(1-sigma_sqd_w)]*[ <w'x'>^2 / (<w'^2> * <x'^2>) ] / mixt_frac }
    * { (1/3)*beta + mixt_frac*( 1 - (2/3)*beta ) };

    where "x" stands for thl, rt, or sclr; "mixt_frac" is the weight of
    Gaussian "plume" 1, and 0 <= beta <= 3.

    The factor { (1/3)*beta + mixt_frac*( 1 - (2/3)*beta ) } does not depend
    on which varable "x" stands for.  The factor is multiplied by 2 and
    defined as width_factor_1.

    The factor
    { 1 - [1/(1-sigma_sqd_w)]*[ (w'x')^2 / (w'^2 * x'^2) ] / mixt_frac }
    depends on which variable "x" stands for.  It is multiplied by 1/2 and
    defined as alpha_x, where "x" stands for thl, rt, or sclr.

    Returns:
        (x_1, x_2, varnce_x_1, varnce_x_2, alpha_x)

    Fortran declarations (ADG1_ADG2_responder_params in src/CLUBB_core/adg1_adg2_3d_luhar_pdf.F90):

        integer, intent(in) :: &
          ngrdcol,  & ! Number of grid columns
          nz          ! Number of vertical level

        !-------------------------- Input Variables --------------------------
        real ( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: &
          xm,          & ! Mean of x (overall)                    [units vary]
          xp2,         & ! Variance of x (overall)                [(units vary)^2]
          wp2,         & ! Variance of w (overall)                [m^2/s^2]
          sqrt_wp2,    & ! Square root of the variance of w       [m/s]
          wpxp,        & ! Covariance of w and x                  [m/s (units vary)]
          w_1_n,       & ! Normalized mean of w (1st PDF comp.)   [-]
          w_2_n,       & ! Normalized mean of w (2nd PDF comp.)   [-]
          mixt_frac,   & ! Mixture fraction                       [-]
          sigma_sqd_w    ! Width of individual w plumes           [-]

        real ( kind = core_rknd ), dimension(ngrdcol), intent(in) :: &
          beta           ! CLUBB tunable parameter beta           [-]

        real ( kind = core_rknd ), intent(in) :: &
          x_tol          ! Tolerance value for x                  [units vary]

        !-------------------------- InOut Variables --------------------------
        type(err_info_type), intent(inout) :: &
          err_info

        !-------------------------- Output Variables --------------------------
        real ( kind = core_rknd ), dimension(ngrdcol,nz), intent(out) :: &
          x_1,        & ! Mean of x (1st PDF component)             [units vary]
          x_2,        & ! Mean of x (2nd PDF component)             [units vary]
          varnce_x_1, & ! Variance of x (1st PDF component)         [(units vary)^2]
          varnce_x_2, & ! Variance of x (2nd PDF component)         [(units vary)^2]
          alpha_x       ! Factor relating to normalized variance for x           [-]
    """
    # TODO: Fortran's responder tolerance diagnostic requires x_tol and err_info,
    # neither carried by this reduced interface; restoring the guard is separate.
    x_1 = xm - wpxp / (sqrt_wp2 * w_2_n)
    x_2 = xm - wpxp / (sqrt_wp2 * w_1_n)

    alpha_x = 0.5 * (1.0 - wpxp ** 2 / ((1.0 - sigma_sqd_w) * wp2 * xp2))
    alpha_x = jnp.clip(alpha_x, zero_threshold, 1.0)

    # variables for a generalization of Chris Golaz' closure
    # varies width of plumes in theta_l, rt
    if hasattr(beta, 'shape') and len(beta.shape) == 1:
        beta_bc = beta[:, None]  # (ngrdcol, 1)
    else:
        beta_bc = beta
    two_thirds_beta = (2.0 / 3.0) * beta_bc
    width_factor_1 = two_thirds_beta + 2.0 * mixt_frac * (1.0 - two_thirds_beta)

    # Vince Larson multiplied original expressions by width_factor_1,2
    #   to generalize scalar skewnesses.  05 Nov 03
    varnce_x_1 = width_factor_1 * xp2 * alpha_x / mixt_frac
    varnce_x_2 = (2.0 - width_factor_1) * xp2 * alpha_x / (1.0 - mixt_frac)

    return x_1, x_2, varnce_x_1, varnce_x_2, alpha_x


def ADG1_pdf_driver(wm, rtm, thlm, um, vm,
                        wp2, rtp2, thlp2, up2, vp2,
                        Skw, wprtp, wpthlp, upwp, vpwp,
                        sqrt_wp2, sigma_sqd_w,
                        beta, mixt_frac_max_mag):
    """Calculates the PDF parameters using the Analytic Double Gaussian 1 (ADG1)
    PDF closure.

    Returns:
        dict with keys w_1, w_2, w_1_n, w_2_n, varnce_w_1, varnce_w_2,
          mixt_frac, rt_1, rt_2, varnce_rt_1, varnce_rt_2, alpha_rt, thl_1,
          thl_2, varnce_thl_1, varnce_thl_2, alpha_thl, u_1, u_2, varnce_u_1,
          varnce_u_2, alpha_u, v_1, v_2, varnce_v_1, varnce_v_2, alpha_v

    Fortran declarations (ADG1_pdf_driver in src/CLUBB_core/adg1_adg2_3d_luhar_pdf.F90):

        !--------------------------- Input Variables ---------------------------
        integer, intent(in) :: &
          nz,           & ! Number of vertical levels
          ngrdcol,      & ! Number of grid columns
          sclr_dim        ! Number of passive scalars

        real( kind = core_rknd ), intent(in), dimension(sclr_dim) :: &
          sclr_tol          ! Threshold(s) on the passive scalars  [units vary]

        real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) ::  &
          wm,          & ! Mean of w-wind comp. (vert. vel.)    [m/s]
          rtm,         & ! Mean of total water mixing ratio     [kg/kg]
          thlm,        & ! Mean of liquid water potential temp. [K]
          um,          & ! Mean of eastward wind                [m/s]
          vm,          & ! Mean of northward wind               [m/s]
          wp2,         & ! Variance of w (overall)              [m^2/s^2]
          rtp2,        & ! Variance of r_t (overall)            [(kg/kg)^2]
          thlp2,       & ! Variance of th_l (overall)           [K^2]
          up2,         & ! Variance of eastward wind (overall)  [m^2/s^2]
          vp2,         & ! Variance of northward wind (overall) [m^2/s^2]
          Skw,         & ! Skewness of w                        [-]
          wprtp,       & ! Covariance of w and r_t              [(kg/kg)(m/s)]
          wpthlp,      & ! Covariance of w and th_l             [K(m/s)]
          upwp,        & ! Covariance of u and w                [m^2/s^2]
          vpwp,        & ! Covariance of v and w                [m^2/s^2]
          sqrt_wp2,    & ! Square root of variance of w         [m/s]
          sigma_sqd_w    ! Width of individual w plumes         [-]

        real( kind = core_rknd ), dimension(ngrdcol), intent(in) ::  &
          beta                    ! CLUBB tunable parameter beta         [-]

        real( kind = core_rknd ), dimension(ngrdcol), intent(in) :: &
          mixt_frac_max_mag ! Maximum allowable mag. of mixt_frac  [-]

        real( kind = core_rknd ), dimension(ngrdcol, nz, sclr_dim), intent(in) ::  &
          sclrm,   & ! Mean of passive scalar (overall)           [units vary]
          sclrp2,  & ! Variance of passive scalar (overall)       [(units vary)^2]
          wpsclrp    ! Covariance of w and passive scalar         [m/s (units vary)]

        logical, intent(in) :: &
          l_scalar_calc    ! Flag to perform calculations for passive scalars

        !--------------------------- Input/Output Variables ---------------------------
        type(err_info_type), intent(inout) :: &
          err_info

        !--------------------------- Output Variables ---------------------------
        real( kind = core_rknd ), dimension(ngrdcol,nz), intent(out) ::  &
          w_1,          & ! Mean of w (1st PDF component)                      [m/s]
          w_2,          & ! Mean of w (2nd PDF component)                      [m/s]
          rt_1,         & ! Mean of r_t (1st PDF component)                  [kg/kg]
          rt_2,         & ! Mean of r_t (2nd PDF component)                  [kg/kg]
          thl_1,        & ! Mean of th_l (1st PDF component)                     [K]
          thl_2,        & ! Mean of th_l (2nd PDF component)                     [K]
          u_1,          & ! Mean of u (eastward wind) (1st PDF component)      [m/s]
          u_2,          & ! Mean of u (eastward wind) (2nd PDF component)      [m/s]
          v_1,          & ! Mean of v (northward wind) (1st PDF component)     [m/s]
          v_2,          & ! Mean of v (northward wind) (2nd PDF component)     [m/s]
          varnce_w_1,   & ! Variance of w (1st PDF component)              [m^2/s^2]
          varnce_w_2,   & ! Variance of w (2nd PDF component)              [m^2/s^2]
          varnce_rt_1,  & ! Variance of r_t (1st PDF component)          [kg^2/kg^2]
          varnce_rt_2,  & ! Variance of r_t (2nd PDF component)          [kg^2/kg^2]
          varnce_thl_1, & ! Variance of th_l (1st PDF component)               [K^2]
          varnce_thl_2, & ! Variance of th_l (2nd PDF component)               [K^2]
          varnce_u_1,   & ! Variance of u wind (1st PDF component)         [m^2/s^2]
          varnce_u_2,   & ! Variance of u wind (2nd PDF component)         [m^2/s^2]
          varnce_v_1,   & ! Variance of v wind (1st PDF component)         [m^2/s^2]
          varnce_v_2,   & ! Variance of v wind (2nd PDF component)         [m^2/s^2]
          mixt_frac,    & ! Mixture fraction (weight of 1st PDF component)       [-]
          alpha_thl,    & ! Factor relating to normalized variance for th_l      [-]
          alpha_rt,     & ! Factor relating to normalized variance for r_t       [-]
          alpha_u,      & ! Factor relating to normalized variance for u wind    [-]
          alpha_v         ! Factor relating to normalized variance for v wind    [-]

        real( kind = core_rknd ), dimension(ngrdcol, nz, sclr_dim), intent(out) ::  &
          sclr_1,        & ! Mean of passive scalar (1st PDF component) [units vary]
          sclr_2,        & ! Mean of passive scalar (2nd PDF component) [units vary]
          varnce_sclr_1, & ! Variance of pass. sclr (1st PDF comp.) [(units vary)^2]
          varnce_sclr_2, & ! Variance of pass. sclr (2nd PDF comp.) [(units vary)^2]
          alpha_sclr       ! Factor relating to normalized variance for sclr     [-]
    """
    # w closure
    (w_1, w_2, w_1_n, w_2_n,
     varnce_w_1, varnce_w_2, mixt_frac) = ADG1_w_closure(
        wm, wp2, Skw, sigma_sqd_w, sqrt_wp2, mixt_frac_max_mag)

    # Responder parameters for rt
    (rt_1, rt_2, varnce_rt_1, varnce_rt_2, alpha_rt) = ADG1_ADG2_responder_params(
        rtm, rtp2, wp2, sqrt_wp2, wprtp, w_1_n, w_2_n, mixt_frac, sigma_sqd_w, beta)

    # Responder parameters for thl
    (thl_1, thl_2, varnce_thl_1, varnce_thl_2, alpha_thl) = ADG1_ADG2_responder_params(
        thlm, thlp2, wp2, sqrt_wp2, wpthlp, w_1_n, w_2_n, mixt_frac, sigma_sqd_w, beta)

    # Responder parameters for u
    (u_1, u_2, varnce_u_1, varnce_u_2, alpha_u) = ADG1_ADG2_responder_params(
        um, up2, wp2, sqrt_wp2, upwp, w_1_n, w_2_n, mixt_frac, sigma_sqd_w, beta)

    # Responder parameters for v
    (v_1, v_2, varnce_v_1, varnce_v_2, alpha_v) = ADG1_ADG2_responder_params(
        vm, vp2, wp2, sqrt_wp2, vpwp, w_1_n, w_2_n, mixt_frac, sigma_sqd_w, beta)

    return {
        'w_1': w_1, 'w_2': w_2,
        'w_1_n': w_1_n, 'w_2_n': w_2_n,
        'varnce_w_1': varnce_w_1, 'varnce_w_2': varnce_w_2,
        'mixt_frac': mixt_frac,
        'rt_1': rt_1, 'rt_2': rt_2,
        'varnce_rt_1': varnce_rt_1, 'varnce_rt_2': varnce_rt_2, 'alpha_rt': alpha_rt,
        'thl_1': thl_1, 'thl_2': thl_2,
        'varnce_thl_1': varnce_thl_1, 'varnce_thl_2': varnce_thl_2, 'alpha_thl': alpha_thl,
        'u_1': u_1, 'u_2': u_2,
        'varnce_u_1': varnce_u_1, 'varnce_u_2': varnce_u_2, 'alpha_u': alpha_u,
        'v_1': v_1, 'v_2': v_2,
        'varnce_v_1': varnce_v_1, 'varnce_v_2': varnce_v_2, 'alpha_v': alpha_v,
    }


# calc_comp_corrs_binormal lives in its Fortran home pdf_utilities.F90; imported here
# (the ADG1 PDF driver `use pdf_utilities` for it) and re-exported under its bare Fortran name.
from clubb_jax.src.CLUBB_core.pdf_utilities import (
    calc_comp_corrs_binormal,
)


def calc_Luhar_params(Skx, wpxp, xp2, x_tol_sqd):
    """For the Luhar closure, this subroutine takes Skx (and Skw) as input and
    outputs the mixture fraction, big_m, and small_m. This code was written
    using the equations and nomenclature of Larson et al. (2002) Appendix
    section e.

    The relationship between skewness of x (Skx), mixture fraction (a), and
    Luhar's small m (m) is given by:

    Skx^2 = ( m^2 * ( m^2 + 3 )^2 / ( m^2 + 1 )^3 )
            * ( 1 - 2*a )^2 / ( a * ( 1 - a ) ).

    Luhar's large M (M) is used to more easily express the factor involving
    the m's:

    M = ( m^2 + 1 )^3 / ( m^2 * ( m^2 + 3 )^2 ).

    The equation involving skewness of x becomes:

    Skx^2 = ( 1 / M ) * ( 1 - 2*a )^2 / ( a * ( 1 - a ) );

    or:

    M * Skx^2 = ( 1 - 2*a )^2 / ( a * ( 1 - a ) ).

    This equation can be rewritten as:

    ( a * ( 1 - a ) ) * M * Skx^2 = ( 1 - 2*a )^2;

    as well as:

    ( a - a^2 ) * M * Skx^2 = 1 - 4*a + 4*a^2;

    and eventually as:

    ( 4 + M * Skx^2 ) * a^2 - ( 4 + M * Skx^2 ) * a + 1 = 0.

    Solving the quadratic equation for a:

    a = (1/2) * ( 1 +- Skx * sqrt( 1 / ( 4/M + Skx^2 ) ) ).

    Since by definition, mu_w_1 >= mu_w_2, a < 0.5 when Skw > 0, the equation
    for mixture fraction is:

    a = (1/2) * ( 1 - Skx * sqrt( 1 / ( 4/M + Skx^2 ) ) ).

    For 3-D Luhar, the variable (w, rt, or theta-l) with the greatest
    magnitude of skewness is used to calculate mixture fraction.  Since it is
    desirable to still have a < 0.5 when Skw > 0 and a > 0.5 when Skw < 0, the
    sign function is used.  The value of Skx is replaced by:

    Skx|_adj = sgn( <w'x'> ) * Skx;

    where

    sgn( <w'x'> ) = | 1 when <w'x'> >= 0
                    | -1 when <w'x'> < 0.

    Since Skx|_adj^2 = ( sgn( <w'x'> ) * Skx )^2 = ( sgn( <w'x'> ) )^2 * Skx^2
    = Skx^2, the equation for mixture fraction is:

    a = (1/2) * ( 1 - sgn( <w'x'> ) * Skx * sqrt( 1 / ( 4/M + Skx^2 ) ) ).

    When using the ADG2 closure or when using the 3-D Luhar closure when the
    variable with the greatest magnitude of skewness is w, Skw = Skx and
    sgn( <w'^2> ) is always equal to 1, reducing the equation to its previous
    form.

    References:
    Vincent E. Larson, Jean-Christophe Golaz, and William R. Cotton, 2002:
    Small-Scale and Mesoscale Variability in Cloudy Boundary Layers: Joint
    Probability Density Functions. J. Atmos. Sci., 59, 3519–3539.

    Returns:
        (mixt_frac, big_m, small_m)

    Fortran declarations (calc_Luhar_params in src/CLUBB_core/adg1_adg2_3d_luhar_pdf.F90):

        !--------------------------- Input Variables ---------------------------
        integer, intent(in) :: &
          nz    ! Number of vertical levels

        real( kind = core_rknd ), dimension(nz), intent(in) :: &
          Skx,  & ! Skewness of x ( <x'^3> / <x'2>^(3/2) )      [-]
          wpxp, & ! Covariance of w and x                       [m/s (x units)]
          xp2     ! Variance of x                               [(x units)^2]

        real( kind = core_rknd ), intent(in) :: &
          x_tol_sqd    ! Tolerance value of x squared           [(x units)^2]

        !--------------------------- Output Variables ---------------------------
        real( kind = core_rknd ), dimension(nz), intent(out) :: &
          mixt_frac, & ! Mixture fraction                       [-]
          big_m,     & ! Luhar's M                              [-]
          small_m      ! Luhar's m                              [-]
    """

    # Width (standard deviation) parameters are non-zero

    # Calculate Luhar's m (small m).
    # If Skx is very small, then small_m will tend to zero which risks
    # divide-by-zero.  To ameliorate this problem, we enforce abs( x_1_n )
    # and abs( x_2_n ) > 0.05.
    # Note:  Luhar's small_m (m) is the only tunable parameter in the Luhar
    #        closure, so this equation can be changed.  However, the value
    #        of m should go toward 0 as Skx goes toward 0 so that the double
    #        Gaussian reduces to a single Gaussian when the distribution is
    #        unskewed.
    # Guarded cube root (cbrt'(0)=inf otherwise; masked by the max floor but would poison grad as 0*inf).
    absS = jnp.abs(Skx)
    cbrt = jnp.where(absS > 0.0, jnp.cbrt(jnp.where(absS > 0.0, absS, 1.0)), 0.0)
    small_m = jnp.maximum((2.0 / 3.0) * cbrt, 0.05)
    # Calculate m^2.
    m_sqd = small_m ** 2
    # Calculate Luhar's M (big M).
    big_m = (1.0 + m_sqd) ** 3 / ((3.0 + m_sqd) ** 2 * m_sqd)

    # Calculate sgn( <w'x'> ).
    sgn = jnp.where(wpxp >= 0.0, 1.0, -1.0)
    # Calculate mixture fraction.
    mixt_frac_vary = 0.5 * (1.0 - sgn * Skx * jnp.sqrt(1.0 / ((4.0 / big_m) + Skx ** 2)))

    vary = xp2 > x_tol_sqd
    mixt_frac = jnp.where(vary, mixt_frac_vary, 0.5)
    big_m = jnp.where(vary, big_m, 0.0)
    small_m = jnp.where(vary, small_m, 0.0)
    return mixt_frac, big_m, small_m


def close_Luhar_pdf(xm, xp2, mixt_frac, small_m, wpxp, x_tol_sqd):
    """For the Luhar closure, this subroutine takes Skx, xm, xp2, and mixt_frac,
    big_m, and small_m (calculated in calc_Luhar_params) as input and outputs
    the PDF component means and variances of a variable x in the joint-PDF
    according to Luhar et al. (1996).  This code was written using the
    equations and nomenclature of Larson et al. (2002) Appendix section e.

    References:
    Vincent E. Larson, Jean-Christophe Golaz, and William R. Cotton, 2002:
    Small-Scale and Mesoscale Variability in Cloudy Boundary Layers: Joint
    Probability Density Functions. J. Atmos. Sci., 59, 3519–3539.

    Returns:
        (sigma_sqd_x_1, sigma_sqd_x_2, varnce_x_1, varnce_x_2, x_1_n, x_2_n,
          x_1, x_2)

    Fortran declarations (close_Luhar_pdf in src/CLUBB_core/adg1_adg2_3d_luhar_pdf.F90):

        !--------------------------- Input Variables ---------------------------
        integer, intent(in) :: &
          nz    ! Number of vertical levels

        real( kind = core_rknd ), dimension(nz), intent(in) :: &
          xm,        & ! Mean (overall) of x, <x>                    [(x units)]
          xp2,       & ! Variance (overall) of x, <x'^2>             [(x units)^2]
          mixt_frac, & ! Mixture fraction                            [-]
          small_m,   & ! Luhar's small m                             [-]
          wpxp         ! Covariance of w and x, <w'x'>               [m/s (x units)]

        real( kind = core_rknd ), intent(in) :: &
          x_tol_sqd    ! Tolerance value of x squared                [(x units)^2]

        !--------------------------- Output Variables ---------------------------
        real( kind = core_rknd ), dimension(nz), intent(out) :: &
          sigma_sqd_x_1, & ! Normalized width parameter of x (1st PDF component) [-]
          sigma_sqd_x_2, & ! Normalized width parameter of x (1st PDF component) [-]
          varnce_x_1,    & ! Variance of x (1st PDF component)         [(x units)^2]
          varnce_x_2,    & ! Variance of x (2nd PDF component)         [(x units)^2]
          x_1_n,         & ! Normalized mean of x (1st PDF component)            [-]
          x_2_n,         & ! Normalized mean of x (2nd PDF component)            [-]
          x_1,           & ! Mean of x (1st PDF component)               [(x units)]
          x_2              ! Mean of x (2nd PDF component)               [(x units)]
    """
    omf = 1.0 - mixt_frac
    msq = 1.0 + small_m ** 2
    sgn = jnp.where(wpxp >= 0.0, 1.0, -1.0)

    # Width (standard deviation) parameters are non-zero

    # Calculate sgn( <w'x'> ).

    # Calculate the square root of the overall variance of x.
    vary = xp2 > x_tol_sqd
    # Normalized width parameter of x in the 1st PDF component.
    sigma_sqd_x_1 = jnp.where(vary, omf / (mixt_frac * msq), 1.0 / msq)
    # Normalized width parameter of x in the 2nd PDF component.
    sigma_sqd_x_2 = jnp.where(vary, mixt_frac / (omf * msq), 1.0 / msq)
    # The variance of x in the 1st PDF component.
    varnce_x_1 = jnp.where(vary, sigma_sqd_x_1 * xp2, 0.0)
    # The variance of x in the 2nd PDF component.
    varnce_x_2 = jnp.where(vary, sigma_sqd_x_2 * xp2, 0.0)
    # Normalized mean of x in the 1st PDF component.
    x_1_n = sgn * small_m * sqrt_clipped(sigma_sqd_x_1)
    # Normalized mean of x in the 2nd PDF component.
    x_2_n = -sgn * small_m * sqrt_clipped(sigma_sqd_x_2)
    sqrt_xp2 = sqrt_clipped(xp2)
    # The mean of x in the 1st PDF component.
    x_1 = jnp.where(vary, xm + sqrt_xp2 * x_1_n, xm)
    # The mean of x in the 2nd PDF component.
    x_2 = jnp.where(vary, xm + sqrt_xp2 * x_2_n, xm)
    return sigma_sqd_x_1, sigma_sqd_x_2, varnce_x_1, varnce_x_2, x_1_n, x_2_n, x_1, x_2


def max_cubic_root(a_coef, b_coef, c_coef, d_coef):
    """Calculates the largest root that results from solving a cubic equation of
    the form a*x^3 + b*x^2 + c*x + d = 0.

    This is done to backsolve for m^2 for the 3-D Luhar closure, given the
    values of mixt_frac and Skx.

    Returns:
        max_root

    Fortran declarations (max_cubic_root in src/CLUBB_core/adg1_adg2_3d_luhar_pdf.F90):

        !--------------------------- Input Variables ---------------------------
        real( kind = core_rknd ), dimension(1), intent(in) :: &
          a_coef, & ! Coefficient a (of x^3) in a*x^3 + b*x^2 + c^x + d = 0    [-]
          b_coef, & ! Coefficient b (of x^2) in a*x^3 + b*x^2 + c^x + d = 0    [-]
          c_coef, & ! Coefficient c (of x) in a*x^3 + b*x^2 + c^x + d = 0      [-]
          d_coef    ! Coefficient d in a*x^3 + b*x^2 + c^x + d = 0             [-]

        !--------------------------- Output Variables ---------------------------
        real( kind = core_rknd ) :: &
          max_root    ! Maximum root that solves the cubic equation            [-]
    """
    from clubb_jax.src.CLUBB_core.calc_roots import cubic_solve, quadratic_solve
    from clubb_jax.src.CLUBB_core.constants_clubb import eps

    a_thresh = 0.001 * jnp.maximum(jnp.maximum(jnp.abs(b_coef), jnp.abs(c_coef)), jnp.abs(d_coef))
    b_thresh = 0.001 * jnp.maximum(jnp.abs(c_coef), jnp.abs(d_coef))

    cr = cubic_solve(a_coef, b_coef, c_coef, d_coef)                      # (..., 3) complex
    r0, r1, r2 = jnp.real(cr[..., 0]), jnp.real(cr[..., 1]), jnp.real(cr[..., 2])
    is_3real = (jnp.abs(jnp.imag(cr[..., 1])) < eps) & (jnp.abs(jnp.imag(cr[..., 2])) < eps)
    cubic_max = jnp.where(is_3real, jnp.maximum(jnp.maximum(r0, r1), r2), r0)

    qr = quadratic_solve(b_coef, c_coef, d_coef)                     # (..., 2) complex
    quad_max = jnp.maximum(jnp.real(qr[..., 0]), jnp.real(qr[..., 1]))

    c_safe = jnp.where(c_coef != 0.0, c_coef, 1.0)
    linear = -d_coef / c_safe

    max_root = jnp.where(jnp.abs(a_coef) > a_thresh, cubic_max,
                         jnp.where(jnp.abs(b_coef) > b_thresh, quad_max, linear))
    return max_root


def backsolve_Luhar_params(Sk_max, Skx, big_m_max, mixt_frac):
    """This subroutine calculates Luhar's big_m and small_m for the variate 'x'
    consistent with the mixture fraction of the variate with the largest
    skewness.

    The relationship between skewness of x (Skx), mixture fraction (a), and
    Luhar's small m (m) is given by:

    Skx^2 = ( m^2 * ( m^2 + 3 )^2 / ( m^2 + 1 )^3 )
            * ( 1 - 2*a )^2 / ( a * ( 1 - a ) ).

    Moving the factor involving mixture fraction to the right-hand side:

    ( ( a * ( 1 - a ) ) / ( 1 - 2*a )^2 ) * Skx^2
    = m^2 * ( m^2 + 3 )^2 / ( m^2 + 1 )^3.

    This can be rewritten as:

    ( ( a * ( 1 - a ) ) / ( 1 - 2*a )^2 ) * Skx^2
    = ( m^6 + 6*m^4 + 9*m^2 ) / ( m^6 + 3*m^4 + 3*m^2 + 1 ).

    Setting alpha = ( ( a * ( 1 - a ) ) / ( 1 - 2*a )^2 ) * Skx^2, the
    equation can be rewritten as:

    ( m^6 + 3*m^4 + 3*m^2 + 1 ) * alpha = m^6 + 6*m^4 + 9*m^2.

    This can be rearranged and rewritten as:

    ( alpha - 1 ) * m^6 + ( 3 * alpha - 6 ) * m^4
    + ( 3 * alpha - 9 ) * m^2 + alpha = 0.

    This can be rewritten again as:

    ( alpha - 1 ) * (m^2)^3 + ( 3 * alpha - 6 ) * (m^2)^2
    + ( 3 * alpha - 9 ) * (m^2) + alpha = 0.

    The goal is to solve for m^2, and then take the square root of m^2 to
    solve for m.  This can be accomplished by using the cubic formula (with
    the l_use_cubic_backsolve option), or else by a quadratic approximation.

    Only the cubic-formula branch is ported: in Fortran,
    l_use_cubic_backsolve is a parameter set to .true., so the quadratic
    approximation is never used.

    Returns:
        (big_m_x, small_m_x)

    Fortran declarations (backsolve_Luhar_params in src/CLUBB_core/adg1_adg2_3d_luhar_pdf.F90):

        !--------------------------- Input Variables ---------------------------
        real( kind = core_rknd ), intent(in) :: &
          Sk_max,    & ! Maximum skewness                                      [-]
          Skx,       & ! Skewness of the variate solving small_m and big_m for [-]
          big_m_max, & ! Luhar's big_m of the variate with maximum skewness    [-]
          mixt_frac    ! Mixture fraction                                      [-]

        !--------------------------- Output Variables ---------------------------
        real( kind = core_rknd ), intent(out) :: &
          big_m_x,   & ! Luhar's big_m for the variate being solved for    [-]
          small_m_x    ! Luhar's small_m for the variate being solved for  [-]
    """
    from clubb_jax.src.CLUBB_core.constants_clubb import eps

    near_half = jnp.abs(mixt_frac - 0.5) < 0.001
    skx_zero = jnp.abs(Skx) < eps

    m_nh = 0.05
    bigm_nh = (1.0 + m_nh ** 2) ** 3 / ((3.0 + m_nh ** 2) ** 2 * m_nh ** 2)

    denom = jnp.where(near_half, 1.0, (1.0 - 2.0 * mixt_frac) ** 2)        # guard the 1/(1-2mf)^2 at mf=0.5
    alpha = mixt_frac * (1.0 - mixt_frac) / denom * Skx ** 2
    alpha_safe = jnp.where(alpha > 0.0, alpha, 1.0)
    bigm_else = 1.0 / alpha_safe
    m2 = max_cubic_root(alpha - 1.0, 3.0 * alpha - 6.0, 3.0 * alpha - 9.0, alpha)
    m_else = jnp.sqrt(jnp.maximum(m2, 0.05 ** 2))

    small_m_x = jnp.where(near_half, m_nh, jnp.where(skx_zero, 0.0, m_else))
    big_m_x = jnp.where(near_half, bigm_nh, jnp.where(skx_zero, jnp.inf, bigm_else))
    return big_m_x, small_m_x


def Luhar_3D_pdf_driver(wm, rtm, thlm, wp2, rtp2, thlp2, Skw, Skrt, Skthl, wprtp, wpthlp):
    """Calculates the PDF parameters using the 3D Luhar PDF closure.

    Returns:
        (w1, w2, rt1, rt2, thl1, thl2, vw1, vw2, vrt1, vrt2, vthl1, vthl2,
          mixt_frac)

    Fortran declarations (Luhar_3D_pdf_driver in src/CLUBB_core/adg1_adg2_3d_luhar_pdf.F90):

        !--------------------------- Input Variables ---------------------------
        integer, intent(in) :: &
          nz

        real( kind = core_rknd ), dimension(nz), intent(in) ::  &
          wm,     & ! Mean of w-wind component (vertical velocity)  [m/s]
          rtm,    & ! Mean of total water mixing ratio              [kg/kg]
          thlm,   & ! Mean of liquid water potential temperature    [K]
          wp2,    & ! Variance of w (overall)                       [m^2/s^2]
          rtp2,   & ! Variance of r_t (overall)                     [(kg/kg)^2]
          thlp2,  & ! Variance of th_l (overall)                    [K^2]
          Skw,    & ! Skewness of w                                 [-]
          Skrt,   & ! Skewness of rt                                [-]
          Skthl,  & ! Skewness of th_l                              [-]
          wprtp,  & ! Covariance of w and r_t                       [(kg/kg)(m/s)]
          wpthlp    ! Covariance of w and th_l                      [K(m/s)]

        !--------------------------- Output Variables ---------------------------
        real( kind = core_rknd ), dimension(nz), intent(out) ::  &
          w_1,          & ! Mean of w (1st PDF component)                      [m/s]
          w_2,          & ! Mean of w (2nd PDF component)                      [m/s]
          rt_1,         & ! Mean of r_t (1st PDF component)                  [kg/kg]
          rt_2,         & ! Mean of r_t (2nd PDF component)                  [kg/kg]
          thl_1,        & ! Mean of th_l (1st PDF component)                     [K]
          thl_2,        & ! Mean of th_l (2nd PDF component)                     [K]
          varnce_w_1,   & ! Variance of w (1st PDF component)              [m^2/s^2]
          varnce_w_2,   & ! Variance of w (2nd PDF component)              [m^2/s^2]
          varnce_rt_1,  & ! Variance of r_t (1st PDF component)          [kg^2/kg^2]
          varnce_rt_2,  & ! Variance of r_t (2nd PDF component)          [kg^2/kg^2]
          varnce_thl_1, & ! Variance of th_l (1st PDF component)               [K^2]
          varnce_thl_2, & ! Variance of th_l (2nd PDF component)               [K^2]
          mixt_frac       ! Mixture fraction (weight of 1st PDF component)       [-]
    """
    w_tsq, rt_tsq, thl_tsq = w_tol_sqd, rt_tol ** 2, thl_tol ** 2

    # calc_Luhar_params for each variable as a potential setter (w uses wp2 as the covariance -> sgn = +1).
    mf_w, bigm_w, m_w = calc_Luhar_params(Skw, wp2, wp2, w_tsq)
    mf_rt, bigm_rt, m_rt = calc_Luhar_params(Skrt, wprtp, rtp2, rt_tsq)
    mf_thl, bigm_thl, m_thl = calc_Luhar_params(Skthl, wpthlp, thlp2, thl_tsq)

    aw, art, athl = jnp.abs(Skw), jnp.abs(Skrt), jnp.abs(Skthl)
    w_set = (aw >= athl) & (aw >= art)
    thl_set = (~w_set) & (athl > aw) & (athl >= art)
    rt_set = (~w_set) & (~thl_set)

    mixt_frac = jnp.where(w_set, mf_w, jnp.where(thl_set, mf_thl, mf_rt))
    setter_Sk = jnp.where(w_set, Skw, jnp.where(thl_set, Skthl, Skrt))
    setter_bigm = jnp.where(w_set, bigm_w, jnp.where(thl_set, bigm_thl, bigm_rt))

    # Each variable: its own m if it is the setter, else backsolve from the setter.
    _, bm_w_resp = backsolve_Luhar_params(setter_Sk, Skw, setter_bigm, mixt_frac)
    _, bm_rt_resp = backsolve_Luhar_params(setter_Sk, Skrt, setter_bigm, mixt_frac)
    _, bm_thl_resp = backsolve_Luhar_params(setter_Sk, Skthl, setter_bigm, mixt_frac)
    m_w_f = jnp.where(w_set, m_w, bm_w_resp)
    m_rt_f = jnp.where(rt_set, m_rt, bm_rt_resp)
    m_thl_f = jnp.where(thl_set, m_thl, bm_thl_resp)

    _, _, vw1, vw2, _, _, w1, w2 = close_Luhar_pdf(wm, wp2, mixt_frac, m_w_f, wp2, w_tsq)
    _, _, vrt1, vrt2, _, _, rt1, rt2 = close_Luhar_pdf(rtm, rtp2, mixt_frac, m_rt_f, wprtp, rt_tsq)
    _, _, vthl1, vthl2, _, _, thl1, thl2 = close_Luhar_pdf(thlm, thlp2, mixt_frac, m_thl_f, wpthlp, thl_tsq)

    return (w1, w2, rt1, rt2, thl1, thl2, vw1, vw2, vrt1, vrt2, vthl1, vthl2, mixt_frac)


def ADG2_pdf_driver(nz, ngrdcol, sclr_dim, sclr_tol,
                    wm, rtm, thlm, wp2, rtp2, thlp2,
                    Skw, wprtp, wpthlp, sqrt_wp2, beta,
                    sclrm, sclrp2, wpsclrp, l_scalar_calc):
    """Calculates the PDF parameters using the Analytic Double Gaussian 2 (ADG2)
    PDF closure.

    Returns:
        (w1, w2, rt1, rt2, thl1, thl2, vw1, vw2, vrt1, vrt2, vthl1, vthl2,
          mixt_frac, alpha_rt, alpha_thl, sigma_sqd_w, sclr1, sclr2,
          varnce_sclr1, varnce_sclr2, alpha_sclr)

    Fortran declarations (ADG2_pdf_driver in src/CLUBB_core/adg1_adg2_3d_luhar_pdf.F90):

        !--------------------------- Input Variables ---------------------------
        integer, intent(in) :: &
          nz,           & ! Number of vertical levels
          ngrdcol,      & ! Number of grid columns
          sclr_dim        ! Number of passive scalars

        real( kind = core_rknd ), intent(in), dimension(sclr_dim) :: &
          sclr_tol          ! Threshold(s) on the passive scalars  [units vary]

        real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) ::  &
          wm,       & ! Mean of w-wind component (vertical velocity)  [m/s]
          rtm,      & ! Mean of total water mixing ratio              [kg/kg]
          thlm,     & ! Mean of liquid water potential temperature    [K]
          wp2,      & ! Variance of w (overall)                       [m^2/s^2]
          rtp2,     & ! Variance of r_t (overall)                     [(kg/kg)^2]
          thlp2,    & ! Variance of th_l (overall)                    [K^2]
          Skw,      & ! Skewness of w                                 [-]
          wprtp,    & ! Covariance of w and r_t                       [(kg/kg)(m/s)]
          wpthlp,   & ! Covariance of w and th_l                      [K(m/s)]
          sqrt_wp2    ! Square root of variance of w                  [m/s]

        real ( kind = core_rknd ), dimension(ngrdcol), intent(in) :: &
          beta       ! CLUBB tunable parameter beta               [-]

        real( kind = core_rknd ), dimension(ngrdcol,nz, sclr_dim), intent(in) ::  &
          sclrm,   & ! Mean of passive scalar (overall)           [units vary]
          sclrp2,  & ! Variance of passive scalar (overall)       [(units vary)^2]
          wpsclrp    ! Covariance of w and passive scalar         [m/s (units vary)]

        logical, intent(in) :: &
          l_scalar_calc    ! Flag to perform calculations for passive scalars

        !--------------------------- Input/Output Variables ---------------------------
        type(err_info_type), intent(inout) :: &
          err_info

        !--------------------------- Output Variables ---------------------------
        real( kind = core_rknd ), dimension(ngrdcol,nz), intent(out) ::  &
          w_1,          & ! Mean of w (1st PDF component)                      [m/s]
          w_2,          & ! Mean of w (2nd PDF component)                      [m/s]
          rt_1,         & ! Mean of r_t (1st PDF component)                  [kg/kg]
          rt_2,         & ! Mean of r_t (2nd PDF component)                  [kg/kg]
          thl_1,        & ! Mean of th_l (1st PDF component)                     [K]
          thl_2,        & ! Mean of th_l (2nd PDF component)                     [K]
          varnce_w_1,   & ! Variance of w (1st PDF component)              [m^2/s^2]
          varnce_w_2,   & ! Variance of w (2nd PDF component)              [m^2/s^2]
          varnce_rt_1,  & ! Variance of r_t (1st PDF component)          [kg^2/kg^2]
          varnce_rt_2,  & ! Variance of r_t (2nd PDF component)          [kg^2/kg^2]
          varnce_thl_1, & ! Variance of th_l (1st PDF component)               [K^2]
          varnce_thl_2, & ! Variance of th_l (2nd PDF component)               [K^2]
          mixt_frac,    & ! Mixture fraction (weight of 1st PDF component)       [-]
          alpha_thl,    & ! Factor relating to normalized variance for th_l      [-]
          alpha_rt,     & ! Factor relating to normalized variance for r_t       [-]
          sigma_sqd_w     ! Width of individual w plumes (ADG1 parameter)        [-]

        real( kind = core_rknd ), dimension(ngrdcol,nz, sclr_dim), intent(out) ::  &
          sclr_1,        & ! Mean of passive scalar (1st PDF component) [units vary]
          sclr_2,        & ! Mean of passive scalar (2nd PDF component) [units vary]
          varnce_sclr_1, & ! Variance of pass. sclr (1st PDF comp.) [(units vary)^2]
          varnce_sclr_2, & ! Variance of pass. sclr (2nd PDF comp.) [(units vary)^2]
          alpha_sclr       ! Factor relating to normalized variance for sclr     [-]
    """
    del nz, ngrdcol
    w_tsq = w_tol_sqd

    mixt_frac, _, small_m_w = calc_Luhar_params(Skw, wp2, wp2, w_tsq)
    _, _, vw1, vw2, w_1_n, w_2_n, w1, w2 = close_Luhar_pdf(wm, wp2, mixt_frac, small_m_w, wp2, w_tsq)
    sigma_sqd_w = jnp.minimum(1.0 / (1.0 + small_m_w ** 2), 0.99)

    rt1, rt2, vrt1, vrt2, alpha_rt = ADG1_ADG2_responder_params(
        rtm, rtp2, wp2, sqrt_wp2, wprtp, w_1_n, w_2_n, mixt_frac, sigma_sqd_w, beta)
    thl1, thl2, vthl1, vthl2, alpha_thl = ADG1_ADG2_responder_params(
        thlm, thlp2, wp2, sqrt_wp2, wpthlp, w_1_n, w_2_n, mixt_frac, sigma_sqd_w, beta)

    if l_scalar_calc:
        sclr1, sclr2, varnce_sclr1, varnce_sclr2, alpha_sclr = jax.vmap(
            lambda sclrm_s, sclrp2_s, wpsclrp_s: ADG1_ADG2_responder_params(
                sclrm_s,
                sclrp2_s,
                wp2,
                sqrt_wp2,
                wpsclrp_s,
                w_1_n,
                w_2_n,
                mixt_frac,
                sigma_sqd_w,
                beta,
            ),
            in_axes=(2, 2, 2),
            out_axes=-1,
        )(sclrm, sclrp2, wpsclrp)
    else:
        sclr1 = jnp.zeros_like(sclrm)
        sclr2 = jnp.zeros_like(sclrm)
        varnce_sclr1 = jnp.zeros_like(sclrm)
        varnce_sclr2 = jnp.zeros_like(sclrm)
        alpha_sclr = jnp.zeros_like(sclrm)

    return (w1, w2, rt1, rt2, thl1, thl2, vw1, vw2, vrt1, vrt2, vthl1, vthl2,
            mixt_frac, alpha_rt, alpha_thl, sigma_sqd_w, sclr1, sclr2,
            varnce_sclr1, varnce_sclr2, alpha_sclr)
