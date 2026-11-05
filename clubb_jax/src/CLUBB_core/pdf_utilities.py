"""JAX port of pdf_utilities.F90 - lognormal <-> normal moment/correlation conversions.

CLUBB's precipitating hydrometeors (r_r, N_r) and the simplified cloud-nuclei
concentration N_cn have an assumed LOGNORMAL marginal in each PDF component. The
upscaled-KK rate functions (KK_upscaled_means) consume the moments of ln(hm) — the
NORMAL-space mean, standard deviation, and correlations — but the PDF setup carries
the LINEAR-space mean/variance. These elemental functions convert between the two:

  mean_L2N   : linear mean   -> mean of ln x          (Garvey 2000, App. B)
  stdev_L2N  : linear var    -> std dev of ln x        (Garvey 2000, App. B)
  corr_NL2NN : corr(x, y)    -> corr(x, ln y)   (x normal, y lognormal; Garvey Eq. B-1)
  corr_LL2NN : corr(x, y)    -> corr(ln x, ln y) (x,y lognormal; Garvey Eq. C-3)
  corr_NN2NL : corr(x, ln y) -> corr(x, y)      (inverse of corr_NL2NN)
  corr_NN2LL : corr(ln x,ln y)->corr(x, y)      (inverse of corr_LL2NN)
  calc_corr_{chi,eta}_x : corr(chi/eta, x) from corr(rt, x), corr(thl, x) (chi=crt rt-cthl thl,
                          eta=crt rt+cthl thl decomposition)
  calc_corr_{rt,thl}_x  : the inverses — corr(rt/thl, x) from corr(chi, x), corr(eta, x)

Porting deviations:
- The Fortran `_dp` duplicate routines are not separately ported; JAX runs with
  x64 enabled and uses the core-rknd behavior where it differs.
- Fortran routines with explicit `nz, ngrdcol` loops broadcast over the shapes
  supplied by callers.
- Correlation clipping uses `jnp.clip` where the Fortran uses if/elseif blocks.
- `sqrt_clipped` is a JAX-only helper for finite inactive-branch gradients; it is
  forward-equivalent to sqrt(max(value,0)) for the non-negative variance inputs
  expected here.

Oracle: pdf_utilities.F90 (the core_rknd variants - mean_L2N keeps its max(.,tiny)
guard; the _dp variants drop it). Inputs use the ratio sigma2_on_mu2 = sigma_x^2/mu_x^2.

mean_L2N, stdev_L2N, corr_NN2NL, corr_NN2LL, calc_corr_chi_x, calc_corr_eta_x are bit-to-bit
verifiable against the f2py API; corr_NL2NN/corr_LL2NN (not exposed) are checked vs Monte-Carlo
and calc_corr_rt_x/calc_corr_thl_x via the chi/eta<->rt/thl round-trip.
See clubb_jax/auto_llm_generated_tests/pytests/unit_tests/test_pdf_utilities.py. All functions are jnp and differentiable.
"""
import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.clubb_precision import core_rknd
from clubb_jax.src.CLUBB_core.advance_helper_module import sqrt_clipped

# constants_clubb.F90 (pdf_utilities.F90 `use constants_clubb`)
from clubb_jax.src.CLUBB_core.constants_clubb import (
    max_mag_correlation as MAX_MAG_CORRELATION,
    chi_tol as _CHI_TOL, eta_tol as _ETA_TOL, rt_tol as _RT_TOL, thl_tol as _THL_TOL,
    eps, min_max_smth_mag,
)
_TINY = jnp.finfo(core_rknd).tiny  # Fortran tiny(mu_x)



def mean_L2N(mu_x, sigma2_on_mu2):
    """For a lognormally-distributed variable x, this function finds the mean of
    ln x (mu_x_n) for the ith component of the PDF, given the mean of x (mu_x)
    and the variance of x (sigma_sqd_x) for the ith component of the PDF. The
    value ln x is distributed normally when x is distributed lognormally.

    References:
    Garvey, P. R., 2000: Probability methods for cost uncertainty analysis.
      Marcel Dekker, 401 pp.
    -- App. B.

    Returns:
        mu_x_n

    Fortran declarations (mean_L2N in src/CLUBB_core/pdf_utilities.F90):

        !--------------------------- Input Variables ---------------------------
        real( kind = core_rknd ), intent(in) ::  &
          mu_x,          & ! Mean of x (ith PDF component)                     [-]
          sigma2_on_mu2    ! Ratio:  sigma_x^2 / mu_x^2 (ith PDF component)    [-]

        !--------------------------- Output Variables ---------------------------
        real( kind = core_rknd ) ::  &
          mu_x_n  ! Mean of ln x (ith PDF component)           [-]
    """
    # Find the mean of ln x for the ith component of the PDF.
    # The max( mu_x / sqrt( 1 + sigma_x^2 / mu_x^2 ), tiny( mu_x ) ) statement
    # is used to prevent taking ln 0, which will produce a result of -infinity.
    # This would happen when mu_x is 0.  However, this code should not be
    # entered when mu_x has a value of 0.
    mu_x_n = jnp.log(jnp.maximum(mu_x / jnp.sqrt(1.0 + sigma2_on_mu2), _TINY))
    return mu_x_n


def stdev_L2N(sigma2_on_mu2):
    """For a lognormally-distributed variable x, this function finds the standard
    deviation of ln x (sigma_x_n) for the ith component of the PDF, given the
    mean of x (mu_x) and the variance of x (sigma_sqd_x) for the ith component
    of the PDF.  The value ln x is distributed normally when x is distributed
    lognormally.

    References:
    Garvey, P. R., 2000: Probability methods for cost uncertainty analysis.
      Marcel Dekker, 401 pp.
    -- App. B.

    Returns:
        sigma_x_n

    Fortran declarations (stdev_L2N in src/CLUBB_core/pdf_utilities.F90):

        !--------------------------- Input Variables ---------------------------
        real( kind = core_rknd ), intent(in) ::  &
          sigma2_on_mu2    ! Ratio:  sigma_x^2 / mu_x^2 (ith PDF component)    [-]

        !--------------------------- Output Variables ---------------------------
        real( kind = core_rknd ) ::  &
          sigma_x_n  ! Standard deviation of ln x (ith PDF component)   [-]
    """
    # Find the standard deviation of ln x for the ith component of the PDF.
    sigma_x_n = jnp.sqrt(jnp.log(1.0 + sigma2_on_mu2))
    return sigma_x_n


def corr_NL2NN(corr_x_y, sigma_y_n, y_sigma2_on_mu2):
    """For a normally-distributed variable x and a lognormally-distributed
    variable y, this function finds the correlation of x and ln y (corr_x_y_n)
    for the ith component of the PDF, given the correlation of x and y
    (corr_x_y) and the standard deviation of ln y (sigma_y_n) for the ith
    component of the PDF.  The value ln y is distributed normally when y is
    distributed lognormally.

    References:
    Garvey, P. R., 2000: Probability methods for cost uncertainty analysis.
      Marcel Dekker, 401 pp.
    -- Eq. B-1.

    Returns:
        corr_x_y_n

    Fortran declarations (corr_NL2NN in src/CLUBB_core/pdf_utilities.F90):

        !--------------------------- Input Variables ---------------------------
        real( kind = core_rknd ), intent(in) :: &
          corr_x_y,        & ! Correlation of x and y (ith PDF component)       [-]
          sigma_y_n,       & ! Standard deviation of ln y (ith PDF component)   [-]
          y_sigma2_on_mu2    ! Ratio:  sigma_y^2 / mu_y^2 (ith PDF component)   [-]

        !--------------------------- Output Variables ---------------------------
        real( kind = core_rknd ) ::  &
          corr_x_y_n  ! Correlation of x and ln y (ith PDF component) [-]
    """
    # Find the correlation of x and ln y for the ith component of the PDF.
    # When sigma_y = 0 and mu_y > 0, y_sigma2_on_mu2 = 0.  This results in
    # sigma_y_n = 0.  The resulting corr_x_y_n is undefined.  However, the
    # divide-by-zero problem needs to be addressed in the code.
    sy_safe = jnp.where(sigma_y_n > 0.0, sigma_y_n, 1.0)
    corr_n = jnp.where(sigma_y_n > 0.0,
                       corr_x_y * sqrt_clipped(y_sigma2_on_mu2) / sy_safe,
                       corr_x_y)
    # Clip the magnitude of the correlation of x and ln y in the ith PDF
    # component, just in case the correlation (ith PDF component) of x and y and
    # the standard deviation (ith PDF component) of ln y are inconsistent,
    # resulting in an unrealizable value for corr_x_y_n.
    corr_x_y_n = jnp.clip(corr_n, -MAX_MAG_CORRELATION, MAX_MAG_CORRELATION)
    return corr_x_y_n


def corr_LL2NN(corr_x_y, sigma_x_n, sigma_y_n, x_sigma2_on_mu2, y_sigma2_on_mu2):
    """For lognormally-distributed variables x and y, this function finds the
    correlation of ln x and ln y (corr_x_y_n) for the ith component of the
    PDF, given the correlation of x and y (corr_x_y), the standard deviation
    of ln x (sigma_x_n), and the standard deviation of ln y (sigma_y_n) for
    the ith component of the PDF.  The value of ln x (or ln y) is distributed
    normally when x (or y) is distributed lognormally.

    References:
    Garvey, P. R., 2000: Probability methods for cost uncertainty analysis.
      Marcel Dekker, 401 pp.
    -- Eq. C-3.

    Returns:
        corr_x_y_n

    Fortran declarations (corr_LL2NN in src/CLUBB_core/pdf_utilities.F90):

        !--------------------------- Input Variables ---------------------------
        real( kind = core_rknd ), intent(in) ::  &
          corr_x_y,        & ! Correlation of x and y (ith PDF component)       [-]
          sigma_x_n,       & ! Standard deviation of ln x (ith PDF component)   [-]
          sigma_y_n,       & ! Standard deviation of ln y (ith PDF component)   [-]
          x_sigma2_on_mu2, & ! Ratio:  sigma_x^2 / mu_x^2 (ith PDF component)   [-]
          y_sigma2_on_mu2    ! Ratio:  sigma_y^2 / mu_y^2 (ith PDF component)   [-]

        !--------------------------- Output Variables ---------------------------
        real( kind = core_rknd ) ::  &
          corr_x_y_n  ! Correlation of ln x and ln y (ith PDF component)  [-]
    """
    # Find the correlation of ln x and ln y for the ith component of the PDF.
    # If sigma_x_n = 0 or sigma_y_n = 0, the value of corr_x_y_n is undefined,
    # so set it to corr_x_y.
    both_pos = (sigma_x_n > 0.0) & (sigma_y_n > 0.0)
    denom = jnp.where(both_pos, sigma_x_n * sigma_y_n, 1.0)
    log_arg = 1.0 + corr_x_y * sqrt_clipped(x_sigma2_on_mu2 * y_sigma2_on_mu2)
    corr_n = jnp.where(both_pos, jnp.log(log_arg) / denom, corr_x_y)
    # Clip the magnitude of the correlation of ln x and ln y.
    corr_x_y_n = jnp.clip(corr_n, -MAX_MAG_CORRELATION, MAX_MAG_CORRELATION)
    return corr_x_y_n


def corr_NN2NL(corr_x_y_n, sigma_y_n, y_sigma2_on_mu2):
    """For a normally-distributed variable x and a lognormally-distributed
    variable y, this function finds the correlation of x and y (corr_x_y) for
    the ith component of the PDF, given the correlation of x and ln y
    (corr_x_y_n) and the standard deviation of ln y (sigma_y_n) for the ith
    component of the PDF.  The value ln y is distributed normally when y is
    distributed lognormally.

    References:
    Garvey, P. R., 2000: Probability methods for cost uncertainty analysis.
      Marcel Dekker, 401 pp.
    -- Eq. B-1.

    Returns:
        corr_x_y

    Fortran declarations (corr_NN2NL in src/CLUBB_core/pdf_utilities.F90):

        !--------------------------- Input Variables ---------------------------
        integer, intent(in) :: &
          nz,      & ! Number of vertical levels
          ngrdcol    ! Number of grid columns

        real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: &
          corr_x_y_n,      & ! Correlation of x and ln y (ith PDF component)    [-]
          sigma_y_n,       & ! Standard deviation of ln y (ith PDF component)   [-]
          y_sigma2_on_mu2    ! Ratio:  sigma_y^2 / mu_y^2 (ith PDF component)   [-]

        !--------------------------- Output Variables ---------------------------
        real( kind = core_rknd ), dimension(ngrdcol,nz), intent(out) ::  &
          corr_x_y  ! Correlation of x and y (ith PDF component) [-]
    """
    # Find the correlation of x and y for the ith component of the PDF.
    corr = jnp.where(sigma_y_n > 0.0,
                     corr_x_y_n * sigma_y_n / jnp.sqrt(jnp.maximum(y_sigma2_on_mu2, _TINY)),
                     corr_x_y_n)
    # Clip the magnitude of the correlation of x and y.
    corr_x_y = jnp.clip(corr, -MAX_MAG_CORRELATION, MAX_MAG_CORRELATION)
    return corr_x_y


def corr_NN2LL(corr_x_y_n, sigma_x_n, sigma_y_n, x_sigma2_on_mu2, y_sigma2_on_mu2):
    """For lognormally-distributed variables x and y, this function finds the
    correlation of x and y (corr_x_y) for the ith component of the PDF, given
    the correlation of ln x and ln y (corr_x_y_n), the standard deviation of
    ln x (sigma_x_n), and the standard deviation of ln y (sigma_y_n) for
    the ith component of the PDF.  The value of ln x (or ln y) is distributed
    normally when x (or y) is distributed lognormally.

    References:
    Garvey, P. R., 2000: Probability methods for cost uncertainty analysis.
      Marcel Dekker, 401 pp.
    -- Eq. C-3.

    Returns:
        corr_x_y

    Fortran declarations (corr_NN2LL in src/CLUBB_core/pdf_utilities.F90):

        !--------------------------- Input Variables ---------------------------
        integer, intent(in) :: &
          nz,      & ! Number of vertical levels
          ngrdcol    ! Number of grid columns

        real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) ::  &
          corr_x_y_n,      & ! Correlation of ln x and ln y (ith PDF component) [-]
          sigma_x_n,       & ! Standard deviation of ln x (ith PDF component)   [-]
          sigma_y_n,       & ! Standard deviation of ln y (ith PDF component)   [-]
          x_sigma2_on_mu2, & ! Ratio:  sigma_x^2 / mu_x^2 (ith PDF component)   [-]
          y_sigma2_on_mu2    ! Ratio:  sigma_y^2 / mu_y^2 (ith PDF component)   [-]

        !--------------------------- Output Variables ---------------------------
        real( kind = core_rknd ), dimension(ngrdcol,nz), intent(out) ::  &
          corr_x_y  ! Correlation of x and y (ith PDF component)  [-]
    """
    both_pos = (sigma_x_n > 0.0) & (sigma_y_n > 0.0)
    denom = jnp.where(both_pos, sqrt_clipped(x_sigma2_on_mu2 * y_sigma2_on_mu2), 1.0)
    corr = jnp.where(both_pos,
                     (jnp.exp(sigma_x_n * sigma_y_n * corr_x_y_n) - 1.0) / denom,
                     corr_x_y_n)
    corr_x_y = jnp.clip(corr, -MAX_MAG_CORRELATION, MAX_MAG_CORRELATION)
    return corr_x_y


def calc_corr_chi_x(crt_i, cthl_i, sigma_rt_i, sigma_thl_i, sigma_chi_i,
                    corr_rt_x_i, corr_thl_x_i):
    """This function calculates the correlation of extended liquid water mixing
    ratio, chi (old s), and a generic variable x, within the ith component of
    the PDF.  The variable chi can be split into mean and turbulent
    components, such that:

    chi = <chi> + chi';

    where < > denotes a mean field an ' denotes a turbulent component.

    The linearized equation for chi' is given in Larson et al. (2001), where
    within the ith component of the PDF:

    chi_(i)' = Coef_rt(i) * r_t(i)' - Coef_thl(i) * th_l(i)'.

    The equation for chi' can be multiplied by x'.  The equation becomes:

    chi'x'_(i) = Coef_rt(i) * r_t'x'_(i) - Coef_thl(i) * th_l'x'_(i).

    Averaging both sides, the covariance <chi'x'> is given by the equation:

    <chi'x'_(i)> = Coef_rt(i) * <r_t'x'_(i)> - Coef_thl(i) * <th_l'x'_(i)>.

    This equation can be rewritten as:

    sigma_chi(i) * sigma_x(i) * corr_chi_x(i)
      = Coef_rt(i) * sigma_rt(i) * sigma_x(i) * corr_rt_x(i)
        - Coef_thl(i) * sigma_thl(i) * sigma_x(i) * corr_thl_x(i).

    This equation can be solved for corr_chi_x(i):

    corr_chi_x(i)
    = Coef_rt(i) * ( sigma_rt(i) / sigma_chi(i) ) * corr_rt_x(i)
      - Coef_thl(i) * ( sigma_thl(i) / sigma_chi(i) ) * corr_thl_x(i).

    The correlation of chi and x within the ith component of the PDF is
    calculated.

    References:
    Eq. (13) and Eq. (14) of Larson, V. E., R. Wood, P. R. Field, J.-C. Golaz,
    T. H. Vonder Haar, W. R. Cotton, 2001:  Systematic Biases in the
    Microphysics and Thermodynamics of Numerical Models That Ignore
    Subgrid-Scale Variability. J. Atmos. Sci., 58, 1117--1128,
    doi:https://doi.org/10.1175/1520-0469(2001)058%3C1117:SBITMA%3E2.0.CO;2.

    Eq. (A29) of Griffin, B. M., 2016:  Improving the Subgrid-Scale
    Representation of Hydrometeors and Microphysical Feedback Effects Using a
    Multivariate PDF.  Doctoral dissertation, University of
    Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
    http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.

    Returns:
        corr_chi_x_i

    Fortran declarations (calc_corr_chi_x in src/CLUBB_core/pdf_utilities.F90):

        !--------------------------- Input Variables ---------------------------
        real( kind = core_rknd ), intent(in) :: &
          crt_i,        & ! Coefficient of r_t for chi (old s) (ith PDF comp.)   [-]
          cthl_i,       & ! Coefficient of th_l for chi (ith PDF comp.)  [(kg/kg)/K]
          sigma_rt_i,   & ! Standard deviation of r_t (ith PDF component)    [kg/kg]
          sigma_thl_i,  & ! Standard deviation of th_l (ith PDF component)       [K]
          sigma_chi_i,  & ! Standard deviation of chi (ith PDF component)    [kg/kg]
          corr_rt_x_i,  & ! Correlation of r_t and x (ith PDF component)         [-]
          corr_thl_x_i    ! Correlation of th_l and x (ith PDF component)        [-]

        !--------------------------- Output Variables ---------------------------
        real( kind = core_rknd ) :: &
          corr_chi_x_i  ! Correlation of chi and x (ith PDF component)   [-]
    """
    sig_safe = jnp.where(sigma_chi_i > _CHI_TOL, sigma_chi_i, 1.0)
    corr = jnp.where(
        sigma_chi_i > _CHI_TOL,
        crt_i * (sigma_rt_i / sig_safe) * corr_rt_x_i
        - cthl_i * (sigma_thl_i / sig_safe) * corr_thl_x_i,
        0.0)
    corr_chi_x_i = jnp.clip(corr, -MAX_MAG_CORRELATION, MAX_MAG_CORRELATION)
    return corr_chi_x_i


def calc_corr_eta_x(crt_i, cthl_i, sigma_rt_i, sigma_thl_i, sigma_eta_i,
                    corr_rt_x_i, corr_thl_x_i):
    """This function calculates the correlation of the variable that is
    orthogonal to extended liquid water mixing ratio in a PDF transformation,
    eta (old t), and a generic variable x, within the ith component of
    the PDF.

    References:
    Eq. (A30) of Griffin, B. M., 2016:  Improving the Subgrid-Scale
    Representation of Hydrometeors and Microphysical Feedback Effects Using a
    Multivariate PDF.  Doctoral dissertation, University of
    Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
    http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.

    Returns:
        corr_eta_x_i

    Fortran declarations (calc_corr_eta_x in src/CLUBB_core/pdf_utilities.F90):

        !--------------------------- Input Variables ---------------------------
        real( kind = core_rknd ), intent(in) :: &
          crt_i,        & ! Coefficient of r_t for chi (old s) (ith PDF comp.)   [-]
          cthl_i,       & ! Coefficient of th_l for chi (ith PDF comp.)  [(kg/kg)/K]
          sigma_rt_i,   & ! Standard deviation of r_t (ith PDF component)    [kg/kg]
          sigma_thl_i,  & ! Standard deviation of th_l (ith PDF component)       [K]
          sigma_eta_i,  & ! Standard deviation of eta (ith PDF component)    [kg/kg]
          corr_rt_x_i,  & ! Correlation of r_t and x (ith PDF component)         [-]
          corr_thl_x_i    ! Correlation of th_l and x (ith PDF component)        [-]
    """
    sig_safe = jnp.where(sigma_eta_i > _ETA_TOL, sigma_eta_i, 1.0)
    corr = jnp.where(
        sigma_eta_i > _ETA_TOL,
        crt_i * (sigma_rt_i / sig_safe) * corr_rt_x_i
        + cthl_i * (sigma_thl_i / sig_safe) * corr_thl_x_i,
        0.0)
    corr_eta_x_i = jnp.clip(corr, -MAX_MAG_CORRELATION, MAX_MAG_CORRELATION)
    return corr_eta_x_i


def compute_mean_binormal(mu_x_1, mu_x_2, mixt_frac):
    """Computes the overall grid-box mean of a binormal distribution from the
    mean of each component

    Returns:
        xm

    Fortran declarations (compute_mean_binormal in src/CLUBB_core/pdf_utilities.F90):

        !--------------------------- Input Variables ---------------------------
        real( kind = core_rknd ), intent(in) :: &
          mu_x_1,    & ! First PDF component mean of 'x'                       [?]
          mu_x_2,    & ! Second PDF component mean of 'x'                      [?]
          mixt_frac    ! Weight of the first PDF component                     [-]

        !--------------------------- Output Variables ---------------------------
        real( kind = core_rknd ) :: &
          xm           ! Mean of 'x' (overall)                                 [?]
    """
    # Calculate the overall mean of x.
    xm = mixt_frac * mu_x_1 + (1.0 - mixt_frac) * mu_x_2
    return xm


def compute_variance_binormal(xm, mu_x_1, mu_x_2, stdev_x_1, stdev_x_2, mixt_frac):
    """Computes the overall grid-box variance of a binormal distribution from the
    variance of each component.

    Returns:
        xp2

    Fortran declarations (compute_variance_binormal in src/CLUBB_core/pdf_utilities.F90):

        !--------------------------- Input Variables ---------------------------
        real( kind = core_rknd ), intent(in) :: &
          xm,        & ! Overall mean of 'x'                                   [?]
          mu_x_1,    & ! First PDF component mean of 'x'                       [?]
          mu_x_2,    & ! Second PDF component mean of 'x'                      [?]
          stdev_x_1, & ! Standard deviation of 'x' in the first PDF component  [?]
          stdev_x_2, & ! Standard deviation of 'x' in the second PDF component [?]
          mixt_frac    ! Weight of the first PDF component                     [-]

        !--------------------------- Output Variables ---------------------------
        real( kind = core_rknd ) :: &
          xp2          ! Variance of 'x' (overall)                             [?^2]
    """
    # Calculate the overall variance of x.
    xp2 = (mixt_frac * ((mu_x_1 - xm) ** 2 + stdev_x_1 ** 2)
           + (1.0 - mixt_frac) * ((mu_x_2 - xm) ** 2 + stdev_x_2 ** 2))
    return xp2


def calc_corr_rt_x(crt_i, sigma_rt_i, sigma_chi_i, sigma_eta_i,
                   corr_chi_x_i, corr_eta_x_i):
    """This function calculates the correlation of rt and x based on the
    correlation of chi and x and the correlation of eta and x.

    Returns:
        corr_rt_x_i

    Fortran declarations (calc_corr_rt_x in src/CLUBB_core/pdf_utilities.F90):

        !--------------------------- Input Variables ---------------------------
        real( kind = core_rknd ), intent(in) :: &
          crt_i,        & ! Coef. of r_t in chi/eta eqns. (ith PDF component)    [-]
          sigma_rt_i,   & ! Standard deviation of r_t (ith PDF component)    [kg/kg]
          sigma_chi_i,  & ! Standard deviation of chi (ith PDF component)    [kg/kg]
          sigma_eta_i,  & ! Standard deviation of eta (ith PDF component)    [kg/kg]
          corr_chi_x_i, & ! Correlation of chi and x (ith PDF component)         [-]
          corr_eta_x_i    ! Correlation of eta and x (ith PDF component)         [-]

        !--------------------------- Output Variables ---------------------------
        real( kind = core_rknd ) :: &
          corr_rt_x_i   ! Correlation of rt and x (ith PDF component)            [-]
    """
    sig_safe = jnp.where(sigma_rt_i > _RT_TOL, sigma_rt_i, 1.0)
    corr = jnp.where(
        sigma_rt_i > _RT_TOL,
        (sigma_eta_i * corr_eta_x_i + sigma_chi_i * corr_chi_x_i) / (2.0 * crt_i * sig_safe),
        0.0)
    corr_rt_x_i = jnp.clip(corr, -MAX_MAG_CORRELATION, MAX_MAG_CORRELATION)
    return corr_rt_x_i


def calc_corr_thl_x(cthl_i, sigma_thl_i, sigma_chi_i, sigma_eta_i,
                    corr_chi_x_i, corr_eta_x_i):
    """This function calculates the correlation of thl and x based on the
    correlation of chi and x and the correlation of eta and x.

    Returns:
        corr_thl_x_i

    Fortran declarations (calc_corr_thl_x in src/CLUBB_core/pdf_utilities.F90):

        !--------------------------- Input Variables ---------------------------
        real( kind = core_rknd ), intent(in) :: &
          cthl_i,       & ! Coef. of thl:  chi/eta eqns. (ith PDF comp.) [(kg/kg)/K]
          sigma_thl_i,  & ! Standard deviation of thl (ith PDF component)        [K]
          sigma_chi_i,  & ! Standard deviation of chi (ith PDF component)    [kg/kg]
          sigma_eta_i,  & ! Standard deviation of eta (ith PDF component)    [kg/kg]
          corr_chi_x_i, & ! Correlation of chi and x (ith PDF component)         [-]
          corr_eta_x_i    ! Correlation of eta and x (ith PDF component)         [-]

        !--------------------------- Output Variables ---------------------------
        real( kind = core_rknd ) :: &
          corr_thl_x_i    ! Correlation of thl and x (ith PDF component)         [-]
    """
    sig_safe = jnp.where(sigma_thl_i > _THL_TOL, sigma_thl_i, 1.0)
    corr = jnp.where(
        sigma_thl_i > _THL_TOL,
        (sigma_eta_i * corr_eta_x_i - sigma_chi_i * corr_chi_x_i) / (2.0 * cthl_i * sig_safe),
        0.0)
    corr_thl_x_i = jnp.clip(corr, -MAX_MAG_CORRELATION, MAX_MAG_CORRELATION)
    return corr_thl_x_i


def smooth_corr_quotient(numerator, denominator, denom_thresh):
    """Smoothly bounded correlation quotient num/denom.

    Two smooth_max operations keep the result a valid correlation: the
    denominator is first raised to at least |num|/max_mag_correlation and then
    to at least denom_thresh.

    Fortran declarations (smooth_corr_quotient in src/CLUBB_core/pdf_utilities.F90):

        ! ---------------- Input Variables ----------------
        integer, intent(in) :: &
          ngrdcol, & ! Number of grid columns
          nz         ! Number of vertical levels

        real ( kind = core_rknd ), intent(in) :: &
          denom_thresh

        real ( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: &
          numerator, &
          denominator

        ! ---------------- Output Variable ----------------
        real ( kind = core_rknd ), dimension(ngrdcol,nz), intent(out) :: &
          quotient
    """
    coef = jnp.minimum(min_max_smth_mag, denom_thresh)

    def _smax(a, b):
        return 0.5 * ((a + b) + jnp.sqrt((a - b) ** 2 + coef ** 2))

    tmp = _smax(jnp.abs(numerator) / MAX_MAG_CORRELATION, denominator)
    tmp = _smax(tmp, denom_thresh)
    quotient = numerator / tmp
    return quotient


def calc_comp_corrs_binormal(xpyp, xm, ym, mu_x_1, mu_x_2, mu_y_1, mu_y_2,
                             sigma_x_1_sqd, sigma_x_2_sqd, sigma_y_1_sqd, sigma_y_2_sqd, mixt_frac):
    """Calculates the PDF component correlations of variables x and y, where
    x and y are both distributed as two-component normals (or binormals).
    The PDF component correlations are set equal to each other.

    The overall covariance of x and y, <x'y'>, can be expressed in terms of
    PDF parameters by integrating over the PDF:

    <x'y'> = INT(-inf:inf) INT(-inf:inf) ( x - <x> ) ( y - <y> ) P(x,y) dy dx;

    where <x> is the overall mean of x, <y> is the overall mean of y, and
    P(x,y) is the equation for the two-component normal PDF of x and y.

    The integral is evaluated, and the equation for <x'y'> is:

    <x'y'> = mixt_frac * ( ( mu_x_1 - <x> ) * ( mu_y_1 - <y> )
                           + corr_x_y_1 * sigma_x_1 * sigma_y_1 )
             + ( 1 - mixt_frac ) * ( ( mu_x_2 - <x> ) * ( mu_y_2 - <y> )
                                     + corr_x_y_2 * sigma_x_2 * sigma_y_2 );

    where mu_x_1 is the mean of x in the 1st PDF component, mu_x_2 is the mean
    of x in the 2nd PDF component, mu_y_1 is the mean of y in the 1st PDF
    component, mu_y_2 is the mean of y in the 2nd PDF component, sigma_x_1 is
    the standard deviation of x in the 1st PDF component, sigma_x_2 is the
    standard deviation of x in the 2nd PDF component, sigma_y_1 is the
    standard deviation of y in the 1st PDF component, sigma_y_2 is the
    standard deviation of y in the 2nd PDF component, corr_x_y_1 is the
    correlation of x and y in the 1st PDF component, corr_x_y_2 is the
    correlation of x and y in the 2nd PDF component, and mixt_frac is the
    mixture fraction (weight of the 1st PDF component).

    This equation can be rewritten as:

    <x'y'> = mixt_frac * ( mu_x_1 - <x> ) * ( mu_y_1 - <y> )
             + mixt_frac * corr_x_y_1 * sigma_x_1 * sigma_y_1
             + ( 1 - mixt_frac ) * ( mu_x_2 - <x> ) * ( mu_y_2 - <y> )
             + ( 1 - mixt_frac ) * corr_x_y_2 * sigma_x_2 * sigma_y_2.

    Setting the two PDF component correlations equal to each other
    (corr_x_y_1 = corr_x_y_2), the equation can be solved for the PDF
    component correlations:

    corr_x_y_1 = corr_x_y_2
    = ( <x'y'> - mixt_frac * ( mu_x_1 - <x> ) * ( mu_y_1 - <y> )
               - ( 1 - mixt_frac ) * ( mu_x_2 - <x> ) * ( mu_y_2 - <y> ) )
      / ( mixt_frac * sigma_x_1 * sigma_y_1
          + ( 1 - mixt_frac ) * sigma_x_2 * sigma_y_2 );

    where -1 <= corr_x_y_1 = corr_x_y_2 <= 1.

    When sigma_x_1 * sigma_y_1 = 0 and sigma_x_2 * sigma_y_2 = 0, at least one
    of x or y are constant within each PDF component, and both PDF component
    correlations are undefined.

    Returns:
        (corr, corr)

    Fortran declarations (calc_comp_corrs_binormal in src/CLUBB_core/pdf_utilities.F90):

        ! ---------------- Input Variables ----------------
        integer, intent(in) :: &
          nz,      & ! Number of vertical levels
          ngrdcol    ! Number of grid columns

        real ( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: &
          xpyp,          & ! Covariance of x and y (overall)    [(x units)(y units)]
          xm,            & ! Mean of x (overall)                [x units]
          ym,            & ! Mean of y (overall)                [y units]
          mu_x_1,        & ! Mean of x (1st PDF component)      [x units]
          mu_x_2,        & ! Mean of x (2nd PDF component)      [x units]
          mu_y_1,        & ! Mean of y (1st PDF component)      [y units]
          mu_y_2,        & ! Mean of y (2nd PDF component)      [y units]
          sigma_x_1_sqd, & ! Variance of x (1st PDF component)  [(x units)^2]
          sigma_x_2_sqd, & ! Variance of x (2nd PDF component)  [(x units)^2]
          sigma_y_1_sqd, & ! Variance of y (1st PDF component)  [(y units)^2]
          sigma_y_2_sqd, & ! Variance of y (2nd PDF component)  [(y units)^2]
          mixt_frac        ! Mixture fraction                   [-]

        ! ---------------- Output Variables ----------------
        real ( kind = core_rknd ), dimension(ngrdcol,nz), intent(out) :: &
          corr_x_y_1, & ! Correlation of x and y (1st PDF component)    [-]
          corr_x_y_2    ! Correlation of x and y (2nd PDF component)    [-]
    """
    # Calculate the numerator of the component correlation equation.
    numerator = (xpyp - mixt_frac * (mu_x_1 - xm) * (mu_y_1 - ym)
                 - (1.0 - mixt_frac) * (mu_x_2 - xm) * (mu_y_2 - ym))
    # Calculate the denominator of the component correlation equation.
    denominator = (mixt_frac * sqrt_clipped(sigma_x_1_sqd * sigma_y_1_sqd)
                   + (1.0 - mixt_frac) * sqrt_clipped(sigma_x_2_sqd * sigma_y_2_sqd))
    # Calculate corr_x_y_1 and corr_x_y_2.
    corr = smooth_corr_quotient(numerator, denominator, eps)
    return corr, corr
