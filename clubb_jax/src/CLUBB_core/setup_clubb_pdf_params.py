"""JAX port of setup_clubb_pdf_params.F90.

Porting deviations:
- Fortran output and inout arguments are returned as values so the routines can
  be used in functional/JIT-friendly call graphs.
- Fortran derived types are represented by JAX-friendly Python dataclasses and
  pytrees; optional ``hydromet_pdf_params`` is represented with ``None`` when
  absent.
- Fortran loops over grid columns, levels, hydrometeors, and PDF variables are
  vectorized with JAX arrays and small Python loops over static dimensions.
- Fortran immediate returns after fatal errors are replaced
  by functional updates to ``err_info``; the caller handles the fatal state
  after leaving the JAX region.
- Stats calls update JAX stats banks instead of calling the
  Fortran netCDF stats API directly.

Hydrometeor PDF component moments:

calc_comp_mu_sigma_hm — the in-precip component means/stdevs (mu_hm_1/2,
sigma_hm_1/2) of a precipitating hydrometeor, which are the rate-function inputs
(KK accr/evap consume mu_rr_i/sigma_rr_i etc.). Oracle: setup_clubb_pdf_params.F90:1653.

The two in-precip component means are solved so the overall mean <hm> and overall
variance <hm'^2> are preserved (Griffin 2015):
  <hm> = a f_p_1 mu_1 + (1-a) f_p_2 mu_2
  <hm'^2> = a f_p_1 (1+omicron Rmax (1+zeta)) mu_1^2
            + (1-a) f_p_2 (1+omicron Rmax) mu_2^2 - <hm>^2
with R = omicron*Rmax the ratio sigma_2^2/mu_2^2, and sigma_1^2/mu_1^2 = R(1+zeta).
This gives a quadratic A mu_1^2 + B mu_1 + C = 0; the root is chosen by sign(mu_thl_1 -
mu_thl_2) so the component with more cloud also has the larger in-precip mean. Minimum
bounds (mu_hm_min_coef) and an "emergency" R-recompute handle small/degenerate cases.
Branches: precip in both components / comp 1 only / comp 2 only / neither.

compute_mean_stdev / norm_transform_mean_stdev (Iter131) — the orchestration that stacks
the per-PDF-variable component means/stdevs (chi, eta, w, Ncn, then the precipitating
hydrometeors) into the (ngrdcol, nzt, pdf_dim) arrays the rate functions index, and
transforms the lognormal variables (Ncn + hydrometeors) to normal (log) space. These are
the linear/normal-space inputs assembled by setup_pdf_parameters_api before the KK rate
calls (oracle setup_clubb_pdf_params.F90:818 + :2942). The PDF-variable index layout is
chi, eta, w, Ncn, <hydrometeors in hydromet-array order>, matching the iiPDF indices
(corr_varnce_module.F90:682). For KK that is [chi, eta, w, Ncn, rr, Nr], pdf_dim = 6.
"""
import jax
import sys
import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.advance_helper_module import calc_xpwp
from clubb_jax.src.CLUBB_core.constants_clubb import (
    cloud_frac_min,
    max_mag_correlation,
    Ncn_tol,
    rc_tol,
    w_tol,
    zero_threshold,
)
from clubb_jax.src.CLUBB_core.clip_explicit import clip_covar, clip_wphydrometp
from clubb_jax.src.CLUBB_core.error_code import CLUBB_FATAL_ERROR, clubb_at_least_debug_level
from clubb_jax.src.CLUBB_core.grid_class import zm2zt, zt2zm
from clubb_jax.src.CLUBB_core.hydromet_pdf_parameter_module import (
    HydrometPdfParameter,
    MAX_HYDROMET_DIM,
    PrecipitationFractions,
)
from clubb_jax.src.CLUBB_core.index_mapping import hydromet2pdf_idx, pdf2hydromet_idx
from clubb_jax.src.CLUBB_core.Nc_Ncn_eqns import Nc_in_cloud_to_Ncnm
from clubb_jax.src.CLUBB_core.parameter_indices import ic_K_hm, iomicron, izeta_vrnce_rat
from clubb_jax.src.CLUBB_core.advance_helper_module import sqrt_clipped
from clubb_jax.src.CLUBB_core.pdf_utilities import (
    compute_mean_binormal,
    compute_variance_binormal,
    corr_NN2NL,
    corr_NN2LL,
    mean_L2N,
    stdev_L2N,
)
from clubb_jax.src.CLUBB_core.matrix_operations import Cholesky_factor
from clubb_jax.src.CLUBB_core.precipitation_fraction import precip_fraction

_MU_HM_MIN_COEF = 0.01   # setup_clubb_pdf_params.F90:2176

# Note:  the ADG2 PDF follows the same standards as the ADG1 PDF.
l_follow_ADG1_PDF_standards = True

# PDF-variable index layout (0-based), matching the iiPDF indices for the KK PDF.
IIPDF_CHI = 0
IIPDF_ETA = 1
IIPDF_W = 2
IIPDF_NCN = 3
# Precipitating hydrometeors follow Ncn in hydromet-array order (rr, Nr for KK).


def setup_pdf_parameters_api(
    gr, nzm, nzt, ngrdcol, pdf_dim,            # In
    hydromet_dim,                              # In
    Nc_in_cloud, cloud_frac, Kh_zm,            # In
    ice_supersat_frac, hydromet, wphydrometp,  # In
    corr_array_n_cloud, corr_array_n_below,    # In
    hm_metadata,                               # In
    pdf_params,                                # In
    clubb_params,                              # In
    iiPDF_type,                                # In
    l_use_precip_frac,                         # In
    l_diagnose_correlations,                   # In
    l_calc_w_corr,                             # In
    l_const_Nc_in_cloud,                       # In
    l_fix_w_chi_eta_correlations,              # In
    err_info,                                  # InOut
    precip_fracs,                              # InOut
    hydromet_pdf_params,                       # Out
    stats,                                     # InOut
):
    """Set up hydrometeor PDF parameters.

    Fortran out arguments are returned in the same order after ``err_info``: ``hydrometp2, mu_x_1_n,
    mu_x_2_n, sigma_x_1_n, sigma_x_2_n, corr_array_1_n, corr_array_2_n, corr_cholesky_mtx_1,
    corr_cholesky_mtx_2, precip_fracs, hydromet_pdf_params, stats``.

    Arguments:
        gr: Grid coordinates, interpolation weights and vertical metrics.
        nzm: Number of momentum levels.
        nzt: Number of thermodynamic levels.
        ngrdcol: Number of grid columns.
        pdf_dim: Number of variables in the multivariate PDF.
        hydromet_dim: Number of precipitating hydrometeor fields.
        Nc_in_cloud: Mean (in-cloud) cloud droplet conc. [num/kg]
        cloud_frac: Cloud fraction [-]
        Kh_zm: Eddy diffusivity coef. on momentum levels [m^2/s]
        ice_supersat_frac: Ice supersaturation fraction [-]
        hydromet: Mean of hydrometeor, hm (overall) (t-levs.) [units]
        wphydrometp: Covariance < w'h_m' > (momentum levels) [(m/s)units]
        corr_array_n_cloud: Prescribed normal space corr. array in cloud [-]
        corr_array_n_below: Prescribed normal space corr. array below cl. [-]
        hm_metadata: Hydrometeor/PDF names, zero-based species indices and tolerances.
        pdf_params: PDF parameters [units vary]
        clubb_params: Column-dependent tunable CLUBB parameters.
        iiPDF_type: Selected option for the two-component normal (double Gaussian) PDF type to
            use for the w, rt, and theta-l (or w, chi, and eta) portion of CLUBB's
            multivariate, two-component PDF.
        l_use_precip_frac: Flag to use precipitation fraction in KK microphysics. The
            precipitation fraction is automatically set to 1 when this flag is turned off.
        l_diagnose_correlations: Diagnose correlations instead of using fixed ones
        l_calc_w_corr: Calculate the correlations between w and the hydrometeors
        l_const_Nc_in_cloud: Use a constant cloud droplet conc. within cloud (K&K)
        l_fix_w_chi_eta_correlations: Use a fixed correlation for s and t Mellor(chi/eta)
        err_info: Per-column error state; return any updated fatal status.
        precip_fracs: Precipitation fractions [-]
        hydromet_pdf_params: Hydrometeor PDF parameters [units vary]
        stats: Immutable statistics state; return its updated value.

    Fortran declarations (setup_pdf_parameters_api in src/CLUBB_core/setup_clubb_pdf_params.F90):

        !------------------------ Input Variables ------------------------
        integer, intent(in) :: &
          nzm,          & ! Number of model momentum vertical grid levels
          nzt,          & ! Number of model thermodynamic vertical grid levels
          pdf_dim,      & ! Number of variables in the correlation array
          ngrdcol,      & ! Number of grid columns
          hydromet_dim    ! Number of hydrometeor species

        type (grid), intent(in) :: gr

        real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
          Nc_in_cloud,       & ! Mean (in-cloud) cloud droplet conc.       [num/kg]
          cloud_frac,        & ! Cloud fraction                            [-]
          ice_supersat_frac    ! Ice supersaturation fraction              [-]

        real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(in) :: &
          Kh_zm                ! Eddy diffusivity coef. on momentum levels [m^2/s]

        real( kind = core_rknd ), dimension(ngrdcol,nzt,hydromet_dim), intent(in) :: &
          hydromet       ! Mean of hydrometeor, hm (overall) (t-levs.) [units]

        real( kind = core_rknd ), dimension(ngrdcol,nzm,hydromet_dim), intent(in) :: &
          wphydrometp    ! Covariance < w'h_m' > (momentum levels)     [(m/s)units]

        real( kind = core_rknd ), dimension(pdf_dim,pdf_dim), intent(in) :: &
          corr_array_n_cloud, & ! Prescribed normal space corr. array in cloud  [-]
          corr_array_n_below    ! Prescribed normal space corr. array below cl. [-]

        type (hm_metadata_type), intent(in) :: &
          hm_metadata

        type(pdf_parameter), intent(in) :: &
          pdf_params    ! PDF parameters                               [units vary]

        real( kind = core_rknd ), dimension(ngrdcol,nparams), intent(in) :: &
          clubb_params    ! Array of CLUBB's tunable parameters    [units vary]

        integer, intent(in) :: &
          iiPDF_type    ! Selected option for the two-component normal (double
                        ! Gaussian) PDF type to use for the w, rt, and theta-l (or
                        ! w, chi, and eta) portion of CLUBB's multivariate,
                        ! two-component PDF.

        logical, intent(in) :: &
          l_use_precip_frac,            & ! Flag to use precipitation fraction in KK microphysics. The
                                          ! precipitation fraction is automatically set to 1 when this
                                          ! flag is turned off.
          l_diagnose_correlations,      & ! Diagnose correlations instead of using fixed ones
          l_calc_w_corr,                & ! Calculate the correlations between w and the hydrometeors
          l_const_Nc_in_cloud,          & ! Use a constant cloud droplet conc. within cloud (K&K)
          l_fix_w_chi_eta_correlations    ! Use a fixed correlation for s and t Mellor(chi/eta)

        type(stats_type), intent(inout) :: &
          stats

        type(err_info_type), intent(inout) :: &
          err_info      ! err_info struct containing err_code and err_header
        !------------------------ Output Variables ------------------------
        real( kind = core_rknd ), dimension(ngrdcol,nzm,hydromet_dim), intent(out) :: &
          hydrometp2    ! Variance of a hydrometeor (overall) (m-levs.)   [units^2]

        real( kind = core_rknd ), dimension(ngrdcol,nzt,pdf_dim), intent(out) :: &
          mu_x_1_n,    & ! Mean array (normal space): PDF vars. (comp. 1) [un. vary]
          mu_x_2_n,    & ! Mean array (normal space): PDF vars. (comp. 2) [un. vary]
          sigma_x_1_n, & ! Std. dev. array (normal space): PDF vars (comp. 1) [u.v.]
          sigma_x_2_n    ! Std. dev. array (normal space): PDF vars (comp. 2) [u.v.]

        real( kind = core_rknd ), dimension(ngrdcol,nzt,pdf_dim,pdf_dim), intent(out) :: &
          corr_array_1_n, & ! Corr. array (normal space) of PDF vars. (comp. 1)  [-]
          corr_array_2_n    ! Corr. array (normal space) of PDF vars. (comp. 2)  [-]

        real( kind = core_rknd ), dimension(ngrdcol,nzt,pdf_dim,pdf_dim), intent(out) :: &
          corr_cholesky_mtx_1, & ! Transposed corr. cholesky matrix, 1st comp. [-]
          corr_cholesky_mtx_2    ! Transposed corr. cholesky matrix, 2nd comp. [-]

        ! This is only an output, but it contains allocated arrays, so we need to treat it as inout
        type(precipitation_fractions), intent(inout) :: &
          precip_fracs           ! Precipitation fractions      [-]

        type(hydromet_pdf_parameter), dimension(ngrdcol,nzt), optional, intent(out) :: &
          hydromet_pdf_params    ! Hydrometeor PDF parameters        [units vary]
    """
    del precip_fracs

    # ---- Begin Code ----

    # Assertion check
    # Check that all hydrometeors are positive otherwise exit the program
    if clubb_at_least_debug_level(0):
        negative_hydromet = jnp.any(hydromet < zero_threshold, axis=(1, 2))
        def print_negative_hydromet(values, info):
            import numpy as np
            for j, k, i in np.argwhere((values < zero_threshold).transpose(2, 1, 0)):
                print(info.column_error_header_host(int(i)), file=sys.stderr)
                print(" at beginning of setup_pdf_parameters_api: ", hm_metadata.hydromet_list[int(j)] + " = ",
                      values[i,k,j], " < ", zero_threshold, " at k = ", k+1, file=sys.stderr)
        jax.lax.cond(
            jnp.any(negative_hydromet),
            lambda _: jax.debug.callback(print_negative_hydromet, hydromet, err_info, ordered=True),
            lambda _: None, operand=None,
        )
        err_info = err_info.set_fatal(mask=negative_hydromet)
        # JAX adaptation: Fortran returns immediately after setting a fatal
        # code. A jitted path cannot branch on that host boolean, so the fatal
        # code is carried functionally and handled by the caller after the JAX
        # region.

    # Setup some of the PDF parameters
    sigma_w_1 = jnp.sqrt(pdf_params.varnce_w_1)
    sigma_w_2 = jnp.sqrt(pdf_params.varnce_w_2)

    # Compute rcm_pdf for use within SILHS
    rcm_pdf = compute_mean_binormal(
        pdf_params.rc_1,
        pdf_params.rc_2,
        pdf_params.mixt_frac,
    )

    # Note on hydrometeor PDF shape:
    # To use a single lognormal over the entire grid level, turn off the
    # l_use_precip_frac flag and set omicron to 1 and zeta_vrnce_rat to 0 in
    # tunable_parameters.in.
    # To use a single delta-lognormal (single lognormal in-precip.), enable the
    # l_use_precip_frac flag and set omicron to 1 and zeta_vrnce_rat to 0 in
    # tunable_parameters.in.
    # Otherwise, with l_use_precip_frac enabled and omicron and zeta_vrnce_rat
    # values that are not 1 and 0, respectively, the PDF shape is a double
    # delta-lognormal (two independent lognormals in-precip.).

    # Calculate precipitation fraction.
    if l_use_precip_frac:
        (
            err_info,
            precip_frac,
            precip_frac_1,
            precip_frac_2,
            precip_frac_tol,
            stats,
        ) = precip_fraction(
            gr, nzt, ngrdcol, hydromet_dim,                                               # In
            hydromet, cloud_frac, pdf_params.cloud_frac_1,                                # In
            hm_metadata.l_mix_rat_hm, hm_metadata.l_frozen_hm, hm_metadata.hydromet_tol,  # In
            pdf_params.cloud_frac_2, ice_supersat_frac,                                   # In
            pdf_params.ice_supersat_frac_1, pdf_params.ice_supersat_frac_2,               # In
            pdf_params.mixt_frac, clubb_params,                                           # In
            err_info,                                                                     # InOut
            stats,                                                                        # InOut
        )
        if clubb_at_least_debug_level(0):
            jax.lax.cond(
                err_info.any_fatal(),
                lambda _: jax.debug.callback(
                    lambda info: print(
                        info.global_error_header_host(),
                        " in setup_pdf_parameters_api after calling precip_fraction",
                        sep="\n",
                        file=sys.stderr,
                    ),
                    err_info, ordered=True,
                ),
                lambda _: None, operand=None,
            )
        # JAX adaptation: preserve the fatal code from precip_fraction, but do
        # not return early inside the traced region.
    else:
        precip_frac = jnp.ones((ngrdcol, nzt))
        precip_frac_1 = jnp.ones((ngrdcol, nzt))
        precip_frac_2 = jnp.ones((ngrdcol, nzt))
        precip_frac_tol = jnp.full((ngrdcol,), cloud_frac_min)

    precip_fracs = PrecipitationFractions(
        ngrdcol=ngrdcol,
        nzt=nzt,
        precip_frac=precip_frac,
        precip_frac_1=precip_frac_1,
        precip_frac_2=precip_frac_2,
    )

    # Calculate <N_cn> from Nc_in_cloud, whether Nc_in_cloud is predicted or
    # based on a prescribed value, and whether the value is constant or varying
    # over the grid level.
    if not l_const_Nc_in_cloud:
        # Ncn varies at each vertical level.
        const_Ncnp2_on_Ncnm2 = jnp.full_like(
            Nc_in_cloud,
            hm_metadata.Ncnp2_on_Ncnm2,
        )
        stdev_const_Ncnp2_on_Ncnm2 = stdev_L2N(const_Ncnp2_on_Ncnm2)
    else:
        # Ncn is constant at each vertical level.
        const_Ncnp2_on_Ncnm2 = jnp.zeros_like(Nc_in_cloud)
        stdev_const_Ncnp2_on_Ncnm2 = jnp.zeros_like(Nc_in_cloud)

    const_corr_chi_Ncn_n_cloud = jnp.full_like(
        Nc_in_cloud,
        corr_array_n_cloud[hm_metadata.iiPDF_Ncn, hm_metadata.iiPDF_chi],
    )
    const_corr_chi_Ncn = corr_NN2NL(
        const_corr_chi_Ncn_n_cloud,
        stdev_const_Ncnp2_on_Ncnm2,
        const_Ncnp2_on_Ncnm2,
    )

    Ncnm = Nc_in_cloud_to_Ncnm(
        pdf_params.chi_1, pdf_params.chi_2, pdf_params.stdev_chi_1,  # In
        pdf_params.stdev_chi_2, pdf_params.mixt_frac, Nc_in_cloud,   # In
        pdf_params.cloud_frac_1, pdf_params.cloud_frac_2,            # In
        const_Ncnp2_on_Ncnm2,                                        # In
        const_corr_chi_Ncn,                                          # In
    )

    # Calculate the overall variance of a precipitating hydrometeor (hm),
    # <hm'^2>.
    hydromet_tol = hm_metadata.hydromet_tol
    hmp2_ip_on_hmm2_ip = hm_metadata.hmp2_ip_on_hmm2_ip
    precip_frac_safe = jnp.where(precip_frac > 0.0, precip_frac, 1.0)
    hydrometp2_zt = jnp.where(
        hydromet >= hydromet_tol[None, None, :],
        ((hmp2_ip_on_hmm2_ip[None, None, :] + 1.0) / precip_frac_safe[:, :, None] - 1.0)
        * hydromet**2,
        0.0,
    )

    # Interpolate the overall variance of a hydrometeor, <hm'^2>, to its home
    # on momentum grid levels.
    hydrometp2 = jnp.stack(
        [
            zt2zm(
                nzm,
                nzt,
                ngrdcol,
                gr,
                hydrometp2_zt[:, :, j],
                zero_threshold,
            )
            .at[:, gr.k_ub_zm]
            .set(0.0)
            for j in range(hydromet_dim)
        ],
        axis=-1,
    )

    wphydrometp_zt = jnp.zeros((ngrdcol, nzt, hydromet_dim))
    wpNcnp_zt = jnp.zeros((ngrdcol, nzt))

    # Calculate correlations involving w and Ncn by first calculating the
    # overall covariance of w and Ncn using the down-gradient approximation.
    if l_calc_w_corr:
        wm_zt = compute_mean_binormal(
            pdf_params.w_1,
            pdf_params.w_2,
            pdf_params.mixt_frac,
        )
        wp2_zt = compute_variance_binormal(
            wm_zt, pdf_params.w_1, pdf_params.w_2,  # In
            sigma_w_1, sigma_w_2,                   # In
            pdf_params.mixt_frac,                   # In
        )

        wphydrometp_zt = jnp.stack(
            [zm2zt(nzm, nzt, ngrdcol, gr, wphydrometp[:, :, j]) for j in range(hydromet_dim)],
            axis=-1,
        )
        wphydrometp_zt = jnp.where(
            hydromet < hydromet_tol[None, None, :],
            0.0,
            wphydrometp_zt,
        )

        clipped_wphydrometp_zt = []
        for j in range(hydromet_dim):
            wphydrometp_zt_j, _wphydrometp_chnge_j = clip_covar(
                nzt, ngrdcol, clip_wphydrometp, wp2_zt, hydrometp2_zt[:, :, j],  # In
                wphydrometp_zt[:, :, j],                                         # InOut
            )
            clipped_wphydrometp_zt.append(wphydrometp_zt_j)
        wphydrometp_zt = jnp.stack(clipped_wphydrometp_zt, axis=-1)

        if clubb_params.ndim == 1:
            Kh_zm_c_K_hm = -clubb_params[ic_K_hm] * Kh_zm
        else:
            Kh_zm_c_K_hm = -clubb_params[:, ic_K_hm][:, None] * Kh_zm

        wpNcnp_zm = calc_xpwp(nzm, nzt, ngrdcol, gr, Kh_zm_c_K_hm, Ncnm)

        # Boundary conditions; We are assuming zero flux.
        wpNcnp_zm = wpNcnp_zm.at[:, 0].set(0.0)
        wpNcnp_zm = wpNcnp_zm.at[:, nzm - 1].set(0.0)

        # Interpolate the covariances to thermodynamic grid levels.
        wpNcnp_zt = zm2zt(nzm, nzt, ngrdcol, gr, wpNcnp_zm)
        wpNcnp_zt = jnp.where(Ncnm <= Ncn_tol, 0.0, wpNcnp_zt)
    else:
        wm_zt = compute_mean_binormal(
            pdf_params.w_1,
            pdf_params.w_2,
            pdf_params.mixt_frac,
        )

    # Unpack CLUBB parameters.
    if clubb_params.ndim == 1:
        omicron = clubb_params[iomicron]
        zeta_vrnce_rat = clubb_params[izeta_vrnce_rat]
    else:
        omicron = clubb_params[:, iomicron][:, None]
        zeta_vrnce_rat = clubb_params[:, izeta_vrnce_rat][:, None]

    hydromets = [
        (
            hydromet[:, :, j],
            hydrometp2_zt[:, :, j],
            hmp2_ip_on_hmm2_ip[j],
            hydromet_tol[j],
        )
        for j in range(hydromet_dim)
    ]

    # Calculate the means and standard deviations involving PDF variables
    # -- w, chi, eta, N_cn, and any precipitating hydrometeors (hm in-precip)
    # -- for each PDF component.
    (
        mu_x_1,
        mu_x_2,
        sigma_x_1,
        sigma_x_2,
        hm_1,
        hm_2,
        sigma2_on_mu2_ip_1,
        sigma2_on_mu2_ip_2,
    ) = compute_mean_stdev(
        pdf_params.chi_1,
        pdf_params.chi_2,
        pdf_params.stdev_chi_1,
        pdf_params.stdev_chi_2,
        pdf_params.stdev_eta_1,
        pdf_params.stdev_eta_2,
        Ncnm,
        hm_metadata.Ncnp2_on_Ncnm2,
        l_const_Nc_in_cloud,
        hydromets,
        pdf_params.thl_1,
        pdf_params.thl_2,
        pdf_params.mixt_frac,
        precip_frac,
        precip_frac_1,
        precip_frac_2,
        precip_frac_tol,
        omicron,
        zeta_vrnce_rat,
        pdf_params.w_1,
        pdf_params.w_2,
        sigma_w_1,
        sigma_w_2,
    )

    # Transform the component means and standard deviations involving
    # precipitating hydrometeors (hm in-precip) and N_cn -- ln hm and
    # ln N_cn -- to normal space for each PDF component.
    mu_x_1_n, mu_x_2_n, sigma_x_1_n, sigma_x_2_n = norm_transform_mean_stdev(
        mu_x_1,
        mu_x_2,
        sigma_x_1,
        sigma_x_2,
        sigma2_on_mu2_ip_1,
        sigma2_on_mu2_ip_2,
        Ncnm,
        hm_1,
        hm_2,
        hydromet_tol,
        l_const_Nc_in_cloud,
    )

    # Calculate the normal space correlations.
    # The normal space correlations are the same as the true correlations
    # except when at least one of the variables involved is a precipitating
    # hydrometeor or Ncn.  In these cases, the normal space correlation
    # involves the natural logarithm of the precipitating hydrometeors, ln hm
    # (for example, ln r_r and ln N_r), and ln N_cn for each PDF component.
    if l_diagnose_correlations:
        from clubb_jax.src.CLUBB_core.diagnose_correlations_module import (
            calc_cholesky_corr_mtx_approx,
            diagnose_correlations,
        )

        corr_cloud = diagnose_correlations(
            pdf_dim, hm_metadata.iiPDF_w + 1,  # In
            corr_array_n_cloud,                # In
            l_calc_w_corr,                     # In
        )
        corr_below = diagnose_correlations(
            pdf_dim, hm_metadata.iiPDF_w + 1,  # In
            corr_array_n_below,                # In
            l_calc_w_corr,                     # In
        )
        corr_array_1_n = jnp.where(
            (rcm_pdf > rc_tol)[:, :, None, None],
            corr_cloud,
            corr_below,
        )
        corr_array_2_n = corr_array_1_n

        def _diag_cholesky(corr_array_n):
            corr_cholesky_mtx, corr_array_n = calc_cholesky_corr_mtx_approx(
                pdf_dim,
                hm_metadata.iiPDF_w + 1,
                corr_array_n,
            )
            return corr_cholesky_mtx, corr_array_n

        corr_cholesky_mtx_1, corr_array_1_n = jax.vmap(
            jax.vmap(_diag_cholesky, in_axes=0),
            in_axes=0,
        )(corr_array_1_n)
        corr_cholesky_mtx_2, corr_array_2_n = jax.vmap(
            jax.vmap(_diag_cholesky, in_axes=0),
            in_axes=0,
        )(corr_array_2_n)
    else:
        if l_fix_w_chi_eta_correlations and not l_calc_w_corr:
            # When the flags are set this way, the correlation matrices do not
            # vary with any vertical values, and instead are determined entirely
            # by prescribed values. This results in there
            # being only two unique correlation matrices, one for when the grid box is in cloud
            # and one for when it is not. So instead of setting up correlation matrices for all
            # grid boxes then calculating their Cholesky decomps, we can simply set up two correlation
            # matrices, one for in cloud and one for out cloud, calculate the corresponding
            # Cholesky decompositions, then use the value of rc at each grid box to determine whether
            # we assign the in cloud or out of cloud matrices to that grid box.
            (
                corr_array_1_n,
                corr_array_2_n,
                corr_cholesky_mtx_1,
                corr_cholesky_mtx_2,
            ) = calc_corr_norm_and_cholesky_factor(
                corr_array_n_cloud,
                corr_array_n_below,
                pdf_params.rc_1,
                pdf_params.rc_2,
                iiPDF_type,
                hm_metadata.iiPDF_chi,
                hm_metadata.iiPDF_eta,
                hm_metadata.iiPDF_w,
                hm_metadata.iiPDF_Ncn,
                l_follow_ADG1_PDF_standards,
            )
        else:
            # The correlation matrices can vary with vertical values. So we need to set the
            # correlation matrices up for each grid box, then find the Cholesky decomp for each
            # grid box individually. This is very computationally expensive.
            pdf2hydromet = tuple(pdf2hydromet_idx(j, hm_metadata) for j in range(pdf_dim))
            pdf_params_corr = {
                "corr_chi_eta_1": pdf_params.corr_chi_eta_1,
                "corr_chi_eta_2": pdf_params.corr_chi_eta_2,
                "corr_w_chi_1": pdf_params.corr_w_chi_1,
                "corr_w_chi_2": pdf_params.corr_w_chi_2,
            }
            (
                corr_array_1_n,
                corr_array_2_n,
            ) = comp_corr_norm(
                mu_x_1,
                mu_x_2,
                sigma_x_1,
                sigma_x_2,
                sigma_x_1_n,
                sigma_x_2_n,
                wm_zt,
                pdf_params.rc_1,
                pdf_params.rc_2,
                pdf_params.mixt_frac,
                precip_frac_1,
                precip_frac_2,
                wpNcnp_zt,
                wphydrometp_zt,
                corr_array_n_cloud,
                corr_array_n_below,
                hm_metadata.iiPDF_chi,
                hm_metadata.iiPDF_eta,
                hm_metadata.iiPDF_w,
                hm_metadata.iiPDF_Ncn,
                pdf2hydromet,
                hydromet_tol,
                Ncn_tol,
                iiPDF_type,
                l_calc_w_corr,
                l_fix_w_chi_eta_correlations,
                pdf_params_corr,
            )

            def _cholesky(corr_array_n):
                _corr_array_scaling, corr_cholesky_mtx, _l_corr_array_scaling = Cholesky_factor(
                    corr_array_n
                )
                return jnp.tril(corr_cholesky_mtx)

            # Compute choleksy factorization for the correlation matrix of 1st PDF component
            corr_cholesky_mtx_1 = jax.vmap(
                jax.vmap(_cholesky, in_axes=0),
                in_axes=0,
            )(corr_array_1_n)
            # Compute choleksy factorization for the correlation matrix of 2nd PDF component
            corr_cholesky_mtx_2 = jax.vmap(
                jax.vmap(_cholesky, in_axes=0),
                in_axes=0,
            )(corr_array_2_n)

    # hydromet_pdf_params is optional, so if it is not present we simply skip
    # the steps to compute it.
    if hydromet_pdf_params is not None:
        # Calculate the true correlations for each PDF component.
        corr_array_1, corr_array_2 = denorm_transform_corr(
            sigma_x_1_n, sigma_x_2_n,                # In
            sigma2_on_mu2_ip_1, sigma2_on_mu2_ip_2,  # In
            corr_array_1_n,                          # In
            corr_array_2_n,                          # In
            hm_metadata.iiPDF_chi,                   # In
            hm_metadata.iiPDF_eta,                   # In
            hm_metadata.iiPDF_w,                     # In
            hm_metadata.iiPDF_Ncn,                   # In
        )
        # Pack the PDF parameters
        hydromet_pdf_params = pack_hydromet_pdf_params(
            nzt, ngrdcol,                  # In
            hydromet_dim, hm_metadata,     # In
            hm_1, hm_2, pdf_dim, mu_x_1,   # In
            mu_x_2, sigma_x_1, sigma_x_2,  # In
            corr_array_1, corr_array_2,    # In
            hydromet_pdf_params,           # Out
        )
    else:
        corr_array_1 = None
        corr_array_2 = None

    if stats.l_sample:
        for j in range(hydromet_dim):
            hm_type = hm_metadata.hydromet_list[j]
            var_name = f"{hm_type[:2].strip()}p2_zt"
            stats = stats.update(var_name, hydrometp2_zt[:, :, j])

        if corr_array_1 is not None and corr_array_2 is not None:
            stats = pdf_param_hm_stats(
                nzt, ngrdcol, pdf_dim, hydromet_dim, hm_metadata,  # In
                hm_1, hm_2,                                        # In
                mu_x_1, mu_x_2,                                    # In
                sigma_x_1, sigma_x_2,                              # In
                corr_array_1, corr_array_2,                        # In
                stats,                                             # InOut
            )

        # Statistics for normal space PDF parameters involving hydrometeors.
        stats = pdf_param_ln_hm_stats(
            nzt, ngrdcol, pdf_dim, hm_metadata,  # In
            mu_x_1_n, mu_x_2_n,                  # In
            sigma_x_1_n, sigma_x_2_n,            # In
            corr_array_1_n,                      # In
            corr_array_2_n,                      # In
            stats,                               # InOut
        )

        if stats.var_on_stats_list("rtp2_from_chi"):
            rtp2_zt_from_chi = compute_rtp2_from_chi(
                pdf_params.stdev_chi_1, pdf_params.stdev_chi_2,                      # In
                pdf_params.stdev_eta_1, pdf_params.stdev_eta_2,                      # In
                pdf_params.rt_1, pdf_params.rt_2,                                    # In
                pdf_params.crt_1, pdf_params.crt_2,                                  # In
                pdf_params.mixt_frac,                                                # In
                corr_array_1_n[:, :, hm_metadata.iiPDF_chi, hm_metadata.iiPDF_eta],  # In
                corr_array_2_n[:, :, hm_metadata.iiPDF_chi, hm_metadata.iiPDF_eta],  # In
            )
            rtp2_zm_from_chi = zt2zm(
                nzm,
                nzt,
                ngrdcol,
                gr,
                rtp2_zt_from_chi,
                zero_threshold,
            )
            stats = stats.update("rtp2_from_chi", rtp2_zm_from_chi)

        stats = stats.update("precip_frac", precip_frac)
        stats = stats.update("precip_frac_1", precip_frac_1)
        stats = stats.update("precip_frac_2", precip_frac_2)
        stats = stats.update("Ncnm", Ncnm)

    if clubb_at_least_debug_level(2):
        # Preserve symmetry/unit-diagonal checks inside compiled physics.
        from clubb_jax.src.CLUBB_core.constants_clubb import eps

        symmetry_tol = 1.0e-6  # Largest allowed asymmetry of a correlation array [-]
        bad_corr = jnp.zeros((ngrdcol,), dtype=bool)
        for corr_array_n in (corr_array_1_n, corr_array_2_n):
            symmetric = jnp.all(
                jnp.abs(corr_array_n - jnp.swapaxes(corr_array_n, -1, -2)) <= symmetry_tol, axis=(1, 2, 3)
            )
            unit_diagonal = jnp.all(
                jnp.abs(jnp.diagonal(corr_array_n, axis1=-2, axis2=-1) - 1.0) <= eps, axis=(1, 2)
            )
            bad_corr = bad_corr | ~symmetric | ~unit_diagonal
        err_info = err_info.set_fatal(mask=bad_corr)
        def print_corr_errors(first, second, info):
            from clubb_jax.src.CLUBB_core.corr_varnce_module import assert_corr_symmetric
            for i in range(ngrdcol):
                for k in range(nzt):
                    ok_first = assert_corr_symmetric(first[i,k])
                    ok_second = assert_corr_symmetric(second[i,k])
                    if not ok_first or not ok_second or info.err_code[i] == CLUBB_FATAL_ERROR:
                        print(info.column_error_header_host(i), file=sys.stderr)
        jax.lax.cond(
            jnp.any(bad_corr),
            lambda _: jax.debug.callback(print_corr_errors, corr_array_1_n, corr_array_2_n, err_info, ordered=True),
            lambda _: None, operand=None,
        )


    if clubb_at_least_debug_level(0):
        jax.lax.cond(
            err_info.any_fatal(),
            lambda _: jax.debug.callback(
                lambda: print(
                    " in setup_pdf_parameters_api after calling assert_corr_symmetric",
                    file=sys.stderr,
                ), ordered=True,
            ),
            lambda _: None, operand=None,
        )

    return (
        err_info,
        hydrometp2,
        mu_x_1_n,
        mu_x_2_n,
        sigma_x_1_n,
        sigma_x_2_n,
        corr_array_1_n,
        corr_array_2_n,
        corr_cholesky_mtx_1,
        corr_cholesky_mtx_2,
        precip_fracs,
        hydromet_pdf_params,
        stats,
    )


def calc_comp_mu_sigma_hm(hmm, hmp2, hmp2_ip_on_hmm2_ip, mixt_frac,
                          precip_frac, precip_frac_1, precip_frac_2,
                          hm_tol, precip_frac_tol, mu_thl_1, mu_thl_2,
                          omicron, zeta_vrnce_rat_in):
    """Calculates the in-precipitation mean and in-precipitation standard
    deviation in both PDF components for a precipitating hydrometeor.

    When precipitation is found in both PDF components (precip_frac_1 > 0 and
    precip_frac_2 > 0), the method that solves for in-precip. mean and
    in-precip. standard deviation in each PDF component, preserving overall
    mean and overall variance, is used.  When precipitation fraction is found
    in one PDF component but not the other one (precip_frac_1 > 0 and
    precip_frac_2 = 0, or precip_frac_1 = 0 and precip_frac_2 > 0), the
    calculation of component in-precip. mean and in-precip. standard deviation
    is simple.  When precipitation is not found in either component
    (precip_frac_1 = 0 and precip_frac_2 = 0), there isn't any precipitation
    found overall (at that grid level).

    DESCRIPTION OF THE METHOD THAT SOLVES FOR TWO IN-PRECIPITATION COMPONENTS
    =========================================================================

    OVERVIEW

    The goal is to calculate the in-precip. mean of the hydrometeor field in
    each PDF component (mu_hm_1 and mu_hm_2) in a scenario when there is
    precipitation found in both PDF components.  The fields provided are the
    overall mean of the hydrometeor, <hm>, the overall variance of the
    hydrometeor, <hm’^2>, the mixture fraction, a, the overall precipitation
    fraction, f_p, and the precipitation fraction in each PDF component
    (f_p_1 and f_p_2).

    The PDF equation for <hm> is:

    <hm> = a * f_p_1 * mu_hm_1 + ( 1- a ) * f_p_2 * mu_hm_2.

    Likewise, the PDF equation for <hm’^2> is:

    <hm’^2> = a * f_p_1 * ( mu_hm_1^2 + sigma_hm_1^2 )
              + ( 1 - a ) * f_p_2 * ( mu_hm_2^2 + sigma_hm_2^2 )
              - <hm>^2;

    where sigma_hm_1 and sigma_hm_2 are the in-precip. standard deviations of
    the hydrometeor field in each PDF component.  This can be rewritten as:

    <hm’^2>
    = a * f_p_1 * ( 1 + sigma_hm_1^2 / mu_hm_1^2 ) * mu_hm_1^2
      + ( 1 - a ) * f_p_2 * ( 1 + sigma_hm_2^2 / mu_hm_2^2 ) * mu_hm_2^2
      - <hm>^2.

    The ratio of sigma_hm_2^2 to mu_hm_2^2 is denoted R:

    R = sigma_hm_2^2 / mu_hm_2^2.

    In order to allow sigma_hm_1^2 / mu_hm_1^2 to have a different ratio, the
    parameter zeta is introduced, such that:

    R * ( 1 + zeta ) = sigma_hm_1^2 / mu_hm_1^2;

    where zeta > -1.  When -1 < zeta < 0, the ratio sigma_hm_2^2 / mu_hm_2^2
    grows at the expense of sigma_hm_1^2 / mu_hm_1^2, which narrows.  When
    zeta = 0, the ratio sigma_hm_1^2 / mu_hm_1^2 is the same as
    sigma_hm_2^2 / mu_hm_2^2.  When zeta > 0, sigma_hm_1^2 / mu_hm_1^2 grows
    at the expense of sigma_hm_2^2 / mu_hm_2^2, which narrows.  The component
    variances are written as:

    sigma_hm_1^2 = R * ( 1 + zeta ) * mu_hm_1^2; and
    sigma_hm_2^2 = R * mu_hm_2^2,

    and the component standard deviations are simply:

    sigma_hm_1 = sqrt( R * ( 1 + zeta ) ) * mu_hm_1; and
    sigma_hm_2 = sqrt( R ) * mu_hm_2.

    The equation for <hm’^2> can be rewritten as:

    <hm’^2> = a * f_p_1 * ( 1 + R * ( 1 + zeta ) ) * mu_hm_1^2
              + ( 1 - a ) * f_p_2 * ( 1 + R ) * mu_hm_2^2
              - <hm>^2.

    HYDROMETEOR IN-PRECIP. VARIANCE:
    THE SPREAD OF THE MEANS VS. THE STANDARD DEVIATIONS

    Part I:  Minimum and Maximum Values for R

    The in-precip. variance of the hydrometeor is accounted for through a
    combination of the variance of each PDF component and the spread between
    the means of each PDF component.  At one extreme, the standard deviation
    of each component could be set to 0 and the in-prccip. variance could be
    accounted for by spreading the PDF component (in-precip.) means far apart.
    The value of R in this scenario would be its minimum possible value, which
    is 0.  At the other extreme, the means of each component could be set
    equal to each other and the in-precip. variance could be accounted for
    entirely by the PDF component (in-precip.) standard deviations.  The value
    of R in this scenario would be its maximum possible value, which is Rmax.

    In order to calculate the value of Rmax, use the equation set but set
    mu_hm_1 = mu_hm_2 and R = Rmax.  When this happens:

    <hm> = ( a * f_p_1 + ( 1- a ) * f_p_2 ) * mu_hm_i;

    and since f_p = a * f_p_1 + ( 1 - a ) * f_p_2:

    mu_hm_i = <hm> / f_p = <hm|_ip>;

    where <hm|_ip> is the in-precip. mean of the hydrometeor.  The equation
    for hydrometeor variance in this scenario becomes:

    <hm’^2> = <hm|_ip>^2 * ( a * f_p_1 * ( 1 + Rmax * ( 1 + zeta ) )
                             + ( 1 - a ) * f_p_2 * ( 1 + Rmax ) )
              - <hm>^2.

    The general equation for the in-precip. variance of a hydrometeor,
    <hm|_ip’^2>, is given by:

    <hm|_ip’^2> = ( <hm’^2> + <hm>^2 - f_p * <hm|_ip>^2 ) / f_p;

    which can be rewritten as:

    <hm’^2> + <hm>^2 = f_p * ( <hm|_ip’^2> + <hm|_ip>^2 ).

    When the above equation is substituted into the modified PDF equation for
    <hm’^2>, Rmax is solved for and the equation is:

    Rmax = ( f_p / ( a * f_p_1 * ( 1 + zeta ) + ( 1 - a ) * f_p_2 ) )
           * ( <hm|_ip’^2> / <hm|_ip>^2 ).

    Here, in the scenario that zeta = 0, both PDF components have the same
    mean and same variance, which reduces the in-precip. distribution to an
    assumed single lognormal, and the above equation reduces to:

    Rmax = <hm|_ip’^2> / <hm|_ip>^2;

    which is what is expected in that case.

    Part II:  Enter omicron

    A parameter is used to prescribe the ratio of R to its maximum value,
    Rmax.  The prescribed parameter is called omicron, where:

    R = omicron * Rmax;

    where 0 <= omicron <= 1.  When omicron = 0, the standard deviation of each
    PDF component is 0, and mu_hm_1 is spread as far away from mu_hm_2 as it
    needs to be to account for the in-precip. variance.  When omicron = 1,
    mu_hm_1 is equal to mu_hm_2, and the standard deviations of the PDF
    components account for all of the in-precip. variance (and when zeta = 0,
    the PDF shape is a single lognormal in-precip.).  At intermediate values
    of omicron, the means of each PDF component are somewhat spread and each
    PDF component has some width.  The modified parameters are listed below.

    The ratio of sigma_hm_2^2 to mu_hm_2^2 is:

    sigma_hm_2^2 / mu_hm_2^2 = omicron * Rmax;

    and the ratio of sigma_hm_1^2 / mu_hm_1^2 is:

    sigma_hm_1^2 / mu_hm_1^2 = omicron * Rmax * ( 1 + zeta ).

    The component variances are written as:

    sigma_hm_1^2 = omicron * Rmax * ( 1 + zeta ) * mu_hm_1^2; and
    sigma_hm_2^2 = omicron * Rmax * mu_hm_2^2,

    and the component standard deviations are simply:

    sigma_hm_1 = sqrt( omicron * Rmax * ( 1 + zeta ) ) * mu_hm_1; and
    sigma_hm_2 = sqrt( omicron * Rmax ) * mu_hm_2.

    The equation set becomes:

    [1] <hm> = a * f_p_1 * mu_hm_1 + ( 1- a ) * f_p_2 * mu_hm_2; and

    [2] <hm’^2>
        = a * f_p_1 * ( 1 + omicron * Rmax * ( 1 + zeta ) ) * mu_hm_1^2
          + ( 1 - a ) * f_p_2 * ( 1 + omicron * Rmax ) * mu_hm_2^2
          - <hm>^2.

    SOLVING THE EQUATION SET FOR MU_HM_1 AND MU_HM_2.

    The above system of two equations can be solved for mu_hm_1 and mu_hm_2.
    All other quantities in the equation set are known quantities.  The
    equation for <hm> is rewritten to isolate mu_hm_2:

    mu_hm_2 = ( <hm> - a * f_p_1 * mu_hm_1 ) / ( ( 1 - a ) * f_p_2 ).

    The above equation is substituted into the equation for <hm’^2>.  The
    equation for <hm’^2> is rewritten, resulting in:

    [ a * f_p_1 * ( 1 + omicron * Rmax * ( 1 + zeta ) )
      + a^2 * f_p_1^2 * ( 1 + omicron * Rmax ) / ( ( 1 - a ) * f_p_2 ) ]
    * mu_hm_1^2
    + [ - 2 * <hm> * a * f_p_1 * ( 1 + omicron * Rmax )
        / ( ( 1 - a ) * f_p_2 ) ] * mu_hm_1
    + [ - ( <hm’^2>
            + ( 1 - ( 1 + omicron * Rmax ) / ( ( 1 - a ) * f_p_2 ) )
              * <hm>^2 ) ]
    = 0.

    This equation is of the form:

    A * mu_hm_1^2 + B * mu_hm_1 + C = 0;

    so the solution for mu_hm_1 is:

    mu_hm_1 = ( -B +/- sqrt( B^2 - 4*A*C ) ) / (2*A);

    where:

    A = a * f_p_1 * ( 1 + omicron * Rmax * ( 1 + zeta ) )
        + a^2 * f_p_1^2 * ( 1 + omicron * Rmax ) / ( ( 1 - a ) * f_p_2 );

    B = - 2 * <hm> * a * f_p_1 * ( 1 + omicron * Rmax )
        / ( ( 1 - a ) * f_p_2 );

    and

    C = - ( <hm’^2>
            + ( 1 - ( 1 + omicron * Rmax ) / ( ( 1 - a ) * f_p_2 ) )
              * <hm>^2 ).

    The signs of the coefficients:

    1) coefficient A is always positive,
    2) coefficient B is always negative (this means that -B is always
       positive), and
    3) coefficient C can be positive, negative, or zero.

    Since ( 1 - ( 1 + omicron * Rmax ) / ( ( 1 - a ) * f_p_2 ) ) * <hm>^2 is
    always negative and <hm’^2> is always positive, the sign of coefficient C
    depends on which term is greater in magnitude.

    When <hm’^2> is greater, the sign of coefficient C is negative.  This
    means that -4*A*C is positive, which in turn means that
    sqrt( B^2 - 4*A*C ) is greater in magnitude than -B.  If the subtraction
    option of the +/- were to be chosen, the value of mu_hm_1 would be
    negative in this scenerio.  So the natural thing to do would be to always
    choose the addition option.  However, this method requires that mu_hm_1
    equals mu_hm_2 when omicron = 1.  When zeta >= 0, this happens when the
    addition option is chosen, but not when the subtraction option is chosen.
    However, when zeta < 0, this happens when the subtraction option is
    chosen, but not when the addition option is chosen.  So, the equation for
    mu_hm_1 becomes:

    mu_hm_1 = ( -B + sqrt( B^2 - 4*A*C ) ) / (2*A); when zeta >= 0; and
    mu_hm_1 = ( -B - sqrt( B^2 - 4*A*C ) ) / (2*A); when zeta < 0.

    Once this is set, of course:

    mu_hm_2 = ( <hm> - a * f_p_1 * mu_hm_1 ) / ( ( 1 - a ) * f_p_2 ).

    The system has been solved and the in-precip. PDF component means have
    been found!

    NOTES

    Note 1:

    The term B^2 - 4*A*C has been analyzed, and mathematically:

    B^2 - 4*A*C >= 0

    always holds true.  Additionally, the minimum value:

    B^2 - 4*A*C = 0,

    can only occur when omicron = 1 and zeta = 0 (or alternatively to
    zeta = 0, Rmax = 0, but this only occurs when <hm|_ip'^2> / <hm|_ip>^2 has
    a value of 0).

    Numerically, when omicron = 1 and zeta = 0, B^2 - 4*A*C can produce very
    small (on the order of epsilon) negative values.  This is due to numerical
    round off error.  When this happens, the erroneous small, negative value
    of B^2 - 4*A*C is simply reset to the value it's supposed to have, which
    is 0.

    Note 2:

    As the value of <hm|_ip'^2> / <hm|_ip>^2 increases and as the value of
    omicron decreases (narrowing the in-precip standard deviations and
    increasing the spread between the in-precip means), a situtation arises
    where the value of one of the component means will become negative.  This
    is because there is a limit to the amount of in-precip variance that can
    be represented by this kind of distribution.  In order to prevent
    out-of-bounds values of mu_hm_1 or mu_hm_2, lower limits will be
    declared, called mu_hm_1_min and mu_hm_2_min.  The value of the
    hydrometeor in-precip. component mean will be limited from going any
    smaller (or negative) at this value.  From there, the value of the other
    hydrometeor in-precip. component mean is easy to calculate.  Then, both
    values will be entered into the calculation of hydrometeor variance, which
    will be rewritten to solve for R.  Then, both the hydrometeor mean and
    hydrometeor variance will be preserved with a valid distribution.

    In this emergency scenario, the value of R is:

    R = ( <hm'^2> + <hm>^2 - a * f_p_1 * mu_hm_1^2
          - ( 1 - a ) * f_p_2 * mu_hm_2^2 )
        / ( a * f_p_1 * ( 1 + zeta ) * mu_hm_1^2
            + ( 1 - a ) * f_p_2 * mu_hm_2^2 ).

    The minimum values of the in-precip. component means are bounded by:

    mu_hm_1_min >= hm_tol / f_p_1; and
    mu_hm_2_min >= hm_tol / f_p_2.

    These are set this way because hm_1 ( = mu_hm_1 * f_p_1 ) and
    hm_2 ( = mu_hm_2 * f_p_2 ) need to have values of at least hm_tol when
    precipitation is found in both PDF components.

    However, an in-precip. component mean value of hm_tol / f_p_1 or
    hm_tol / f_p_2 often produces a distribution where one component centers
    around values that are too small to be a good match with data taken from
    Large Eddy Simulations (LES).  It is desirable to increase the minimum
    threshold of mu_hm_1 and mu_hm_2.

    As the minimum threshold increases, the value of the in-precip. component
    mean that is from the component that is not being set to the minimum
    threshold decreases.  If the minimum threshold were to be boosted as high
    as <hm> / f_p (in most cases, <hm> / f_p >> hm_tol / f_p_i), both
    components would have a value of <hm> / f_p.  The minimum threshold should
    not be set this high.

    Additionally, the minimum threshold for one in-precip. component mean
    cannot be set so high as to drive the other in-precip. component mean
    below hm_tol / f_p_i.  (This doesn't come into play unless <hm> is close
    to hm_tol.)  The upper limit for the in-precip. mean values are:

    mu_hm_1|_(upper. lim.) = ( <hm> - ( 1 - a ) * f_p_2 * ( hm_tol / f_p_2 ) )
                             / ( a * f_p_1 ); and

    mu_hm_2|_(upper. lim.) = ( <hm> - a * f_p_1 * ( hm_tol / f_p_1 ) )
                             / ( ( 1 - a ) * f_p_2 );

    which reduces to:

    mu_hm_1|_(upper. lim.) = ( <hm> - ( 1 - a ) * hm_tol ) / ( a * f_p_1 );
    and
    mu_hm_2|_(upper. lim.) = ( <hm> - a * hm_tol ) / ( ( 1 - a ) * f_p_2 ).

    An appropriate minimum value for mu_hm_1 can be set by:

    mu_hm_1_min = | min( hm_tol / f_p_1
                  |      + mu_hm_min_coef * ( <hm> / f_p - hm_tol / f_p_1 ),
                  |      ( <hm> - ( 1 - a ) * hm_tol ) / ( a * f_p_1 ) );
                  |    where <hm> / f_p > hm_tol / f_p_1;
                  | hm_tol / f_p_1;
                  |    where <hm> / f_p <= hm_tol / f_p_1;

    and similarly for mu_hm_2:

    mu_hm_2_min = | min( hm_tol / f_p_2
                  |      + mu_hm_min_coef * ( <hm> / f_p - hm_tol / f_p_2 ),
                  |      ( <hm> - a * hm_tol ) / ( ( 1 - a ) * f_p_2 ) );
                  |    where <hm> / f_p > hm_tol / f_p_2;
                  | hm_tol / f_p_2;
                  |    where <hm> / f_p <= hm_tol / f_p_2;

    where mu_hm_min_coef is a coefficient that has a value
    0 <= mu_hm_min_coef < 1.  When the value of mu_hm_min_coef is 0,
    mu_hm_1_min reverts to hm_tol / f_p_1 and mu_hm_2_min reverts to
    hm_tol / f_p_2.  An appropriate value for mu_hm_min_coef should be small,
    such as 0.01 - 0.05.

    Note 3:

    When the value of zeta >= 0, the value of mu_hm_1 tends to be larger than
    the value of mu_hm_2.  Likewise when the value of zeta < 0, the value of
    mu_hm_2 tends to be larger than the value of mu_hm_1.  Since most cloud
    water and cloud fraction tends to be found in PDF component 1, it is
    advantageous to have the larger in-precip. component mean of the
    hydrometeor also found in PDF component 1.  The recommended value of zeta
    is a value greater than or equal to 0.

    Update:

    In order to better represent the increase in <th_l'^2> near the ground
    from the evaporation of rain, the code is modified to tend towards a
    negative correlation of th_l and hm.  When mu_thl_1 <= mu_thl_2,
    mu_hm_1 >= mu_hm_2; otherwise, when mu_thl_1 > mu_thl_2,
    mu_hm_1 <= mu_hm_2.

    In the original derivation, mu_hm_1 >= mu_hm_2 when zeta >= 0, and
    mu_hm_1 < mu_hm_2 when zeta < 0, where zeta is a tunable or adjustable
    parameter.  In order to allow the relationship of mu_hm_1 to mu_hm_2 to
    depend on the relationship of mu_thl_1 to mu_thl_2, the value of zeta must
    also depend on the relationship of mu_thl_1 to mu_thl_2.

    When mu_thl_1 <= mu_thl_2:

    The relationship of mu_hm_1 to mu_hm_2 is mu_hm_1 >= mu_hm_2, so
    zeta >= 0.  The tunable value of zeta is referred to as zeta_in.  When
    zeta_in is already greater than 0 (meaning sigma_hm_1^2 / mu_hm_1^2 is
    greater than sigma_hm_2^2 / mu_hm_2^2), zeta is simply set to zeta_in.  In
    other words, the component with the greater mean value of the hydrometeor
    also has the greater value of component variance to the square of the
    component mean.  However, when zeta_in is less than 0, zeta must be
    adjusted to be greater than 0.  The following equation is used to set the
    value of zeta:

           | zeta_in; when zeta_in >= 0
    zeta = |
           | ( 1 / ( 1 + zeta_in ) ) - 1; when zeta_in < 0.

    Previously, when zeta_in < 0, mu_hm_1 < mu_hm_2, and
    sigma_hm_1^2 / mu_hm_1^2 < sigma_hm_2^2 / mu_hm_2^2.  Now that zeta has to
    be greater than 0, sigma_hm_1^2 / mu_hm_1^2 > sigma_hm_2^2 / mu_hm_2^2.
    The ratio of the greater variance-over-mean-squared to the smaller
    variance-over-mean-squared remains the same when using the equation for
    zeta listed above.

    When mu_thl_1 > mu_thl_2:

    The relationship of mu_hm_1 to mu_hm_2 is mu_hm_1 <= mu_hm_2, so
    zeta <= 0.  When zeta_in is already less than 0 (meaning
    sigma_hm_1^2 / mu_hm_1^2 is less than sigma_hm_2^2 / mu_hm_2^2), zeta is
    simply set to zeta_in.  In other words, the component with the greater
    mean value of the hydrometeor also has the greater value of component
    variance to the square of the component mean.  However, when zeta_in is
    greater than 0, zeta must be adjusted to be less than 0.  The following
    equation is used to set the value of zeta:

           | zeta_in; when zeta_in <= 0
    zeta = |
           | ( 1 / ( 1 + zeta_in ) ) - 1; when zeta_in > 0.

    Previously, when zeta_in > 0, mu_hm_1 > mu_hm_2, and
    sigma_hm_1^2 / mu_hm_1^2 > sigma_hm_2^2 / mu_hm_2^2.  Now that zeta has to
    be less than 0, sigma_hm_1^2 / mu_hm_1^2 < sigma_hm_2^2 / mu_hm_2^2.
    The ratio of the greater variance-over-mean-squared to the smaller
    variance-over-mean-squared remains the same when using the equation for
    zeta listed above.

    Brian Griffin; February 2015.

    Returns:
        (mu_hm_1, mu_hm_2, sigma_hm_1, sigma_hm_2, hm_1, hm_2, s1r, s2r)

    Fortran declarations (calc_comp_mu_sigma_hm in src/CLUBB_core/setup_clubb_pdf_params.F90):

        !--------------------------- Input Variables ---------------------------
        integer, intent(in) :: &
          nzt,     & ! Number of model thermodynamic vertical grid levels
          ngrdcol    ! Number of grid columns

        real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
          hmm,           & ! Hydrometeor mean (overall), <hm>           [hm units]
          hmp2,          & ! Hydrometeor variance (overall), <hm'^2>    [hm units^2]
          mixt_frac,     & ! Mixture fraction                           [-]
          precip_frac,   & ! Precipitation fraction (overall)           [-]
          precip_frac_1, & ! Precipitation fraction (1st PDF component) [-]
          precip_frac_2, & ! Precipitation fraction (2nd PDF component) [-]
          mu_thl_1,      & ! Mean of th_l (1st PDF component)           [K]
          mu_thl_2         ! Mean of th_l (2nd PDF component)           [K]

        real( kind = core_rknd ), intent(in) :: &
          hmp2_ip_on_hmm2_ip, & ! Ratio <hm|_ip'^2> / <hm|_ip>^2                 [-]
          hm_tol                ! Tolerance value of hydrometeor          [hm units]

        real( kind = core_rknd ), dimension(ngrdcol), intent(in) :: &
          precip_frac_tol       ! Min. precip. frac. when hydromet. are present  [-]

        real( kind = core_rknd ), dimension(ngrdcol), intent(in) :: &
          omicron,           & ! Relative width parameter, omicron = R / Rmax    [-]
          zeta_vrnce_rat_in    ! Width parameter for sigma_hm_1^2 / mu_hm_1^2    [-]

        !--------------------------- Output Variables ---------------------------
        real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(out) :: &
          mu_hm_1,    & ! Mean of hm (1st PDF component) in-precip (ip)      [hm un]
          mu_hm_2,    & ! Mean of hm (2nd PDF component) ip                  [hm un]
          sigma_hm_1, & ! Standard deviation of hm (1st PDF component) ip    [hm un]
          sigma_hm_2, & ! Standard deviation of hm (2nd PDF component) ip    [hm un]
          hm_1,       & ! Mean of hm (1st PDF component)                     [hm un]
          hm_2          ! Mean of hm (2nd PDF component)                     [hm un]

        real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(out) :: &
          sigma_hm_1_sqd_on_mu_hm_1_sqd, & ! Ratio sigma_hm_1**2 / mu_hm_1**2    [-]
          sigma_hm_2_sqd_on_mu_hm_2_sqd    ! Ratio sigma_hm_2**2 / mu_hm_2**2    [-]
    """
    pftol = precip_frac_tol[:, None]
    oma = 1.0 - mixt_frac

    # Branch masks (mirror the Fortran if/elseif ladder).
    both = (hmm >= hm_tol) & (precip_frac_1 >= pftol) & (precip_frac_2 >= pftol)
    comp1 = (~both) & (hmm >= hm_tol) & (precip_frac_1 >= pftol)
    comp2 = (~both) & (~comp1) & (hmm >= hm_tol) & (precip_frac_2 >= pftol)

    # Safe denominators so unselected branches never divide by zero / poison gradients.
    fp1_s = jnp.where(precip_frac_1 > 0.0, precip_frac_1, 1.0)
    fp2_s = jnp.where(precip_frac_2 > 0.0, precip_frac_2, 1.0)
    a_s = jnp.where(mixt_frac > 0.0, mixt_frac, 1.0)
    oma_s = jnp.where(oma > 0.0, oma, 1.0)

    # Adjust the value of zeta based on the relationship of mu_thl_1 to
    # mu_thl_2.
    thl_le = mu_thl_1 <= mu_thl_2
    zeta_flip = (1.0 / (1.0 + zeta_vrnce_rat_in)) - 1.0
    zeta = jnp.where(
        zeta_vrnce_rat_in >= 0.0,
        jnp.where(thl_le, zeta_vrnce_rat_in, zeta_flip),
        jnp.where(thl_le, zeta_flip, zeta_vrnce_rat_in),
    )

    # Calculate the value of Rmax.
    # Rmax = ( f_p / ( a * f_p_1 * ( 1 + zeta ) + ( 1 - a ) * f_p_2 ) )
    #        * ( <hm|_ip’^2> / <hm|_ip>^2 ).
    # The parameter zeta is written in the code as zeta_vrnce_rat.
    # Guard the Rmax denominator: at no-precip points fp1=fp2=0 -> 0/0 = nan (forward and grad),
    # which poisons even though the both-precip branch is unselected there.
    den_R = mixt_frac * precip_frac_1 * (1.0 + zeta) + oma * precip_frac_2
    Rmax = (precip_frac / jnp.where(den_R > 0.0, den_R, 1.0)) * hmp2_ip_on_hmm2_ip
    oR = omicron * Rmax
    # Calculate the value of coefficient A.
    # A = a * f_p_1 * ( 1 + omicron * Rmax * ( 1 + zeta ) )
    #     + a^2 * f_p_1^2 * ( 1 + omicron * Rmax ) / ( ( 1 - a ) * f_p_2 ).
    coef_A = (mixt_frac * precip_frac_1 * (1.0 + oR * (1.0 + zeta))
              + mixt_frac ** 2 * precip_frac_1 ** 2 * (1.0 + oR) / (oma_s * fp2_s))
    # Calculate the value of coefficient B.
    # B = - 2 * <hm> * a * f_p_1 * ( 1 + omicron * Rmax )
    #     / ( ( 1 - a ) * f_p_2 ).
    coef_B = -2.0 * hmm * mixt_frac * precip_frac_1 * (1.0 + oR) / (oma_s * fp2_s)
    # Calculate the value of coefficient C.
    # C = - ( <hm’^2>
    #         + ( 1 - ( 1 + omicron * Rmax ) / ( ( 1 - a ) * f_p_2 ) )
    #           * <hm>^2 ).
    coef_C = -(hmp2 + (1.0 - (1.0 + oR) / (oma_s * fp2_s)) * hmm ** 2)
    # Calculate value of B^2 - 4*A*C.
    disc = jnp.maximum(coef_B ** 2 - 4.0 * coef_A * coef_C, 0.0)
    coef_A_s = jnp.where(coef_A != 0.0, coef_A, 1.0)
    # Calculate the mean (in-precip.) of the hydrometeor in the 1st PDF
    # component.
    mu1 = jnp.where(thl_le, (-coef_B + sqrt_clipped(disc)) / (2.0 * coef_A_s),
                    (-coef_B - sqrt_clipped(disc)) / (2.0 * coef_A_s))
    # Calculate the mean (in-precip.) of the hydrometeor in the 2nd PDF
    # component.
    mu2 = (hmm - mixt_frac * precip_frac_1 * mu1) / (oma_s * fp2_s)
    R = oR

    # Calculate the minimum allowable value of the mean (in-precip.) of a
    # hydrometeor in the 1st PDF component.
    hm_ip = hmm / jnp.where(precip_frac > 0.0, precip_frac, 1.0)
    mu1_min = jnp.where(hm_ip > hm_tol / fp1_s,
                        jnp.minimum(hm_tol / fp1_s + _MU_HM_MIN_COEF * (hm_ip - hm_tol / fp1_s),
                                    (hmm - oma * hm_tol) / (a_s * fp1_s)),
                        hm_tol / fp1_s)
    # Calculate the minimum allowable value of the mean (in-precip.) of a
    # hydrometeor in the 2nd PDF component.
    mu2_min = jnp.where(hm_ip > hm_tol / fp2_s,
                        jnp.minimum(hm_tol / fp2_s + _MU_HM_MIN_COEF * (hm_ip - hm_tol / fp2_s),
                                    (hmm - mixt_frac * hm_tol) / (oma_s * fp2_s)),
                        hm_tol / fp2_s)

    def _emergency_R(m1, m2):
        # Calculate new R.
        num = (hmp2 + hmm ** 2 - mixt_frac * precip_frac_1 * m1 ** 2 - oma * precip_frac_2 * m2 ** 2)
        den = (mixt_frac * precip_frac_1 * (1.0 + zeta) * m1 ** 2 + oma * precip_frac_2 * m2 ** 2)
        return jnp.maximum(num / jnp.where(den != 0.0, den, 1.0), 0.0)

    # If mu_hm_1 is less than mu_hm_1_min, then set mu_hm_1 equal to
    # mu_hm_1_min.  Calculate mu_hm_2 to preserve overall mean, and calculate
    # a new value of R to preserve overall variance.
    mu1_e1 = mu1_min
    mu2_e1 = (hmm - mixt_frac * precip_frac_1 * mu1_e1) / (oma_s * fp2_s)
    R_e1 = _emergency_R(mu1_e1, mu2_e1)
    # If mu_hm_2 is less than mu_hm_2_min, then set mu_hm_2 equal to
    # mu_hm_2_min.  Calculate mu_hm_1 to preserve overall mean, and calculate
    # a new value of R to preserve overall variance.
    mu2_e2 = mu2_min
    mu1_e2 = (hmm - oma * precip_frac_2 * mu2_e2) / (a_s * fp1_s)
    R_e2 = _emergency_R(mu1_e2, mu2_e2)

    use_e1 = mu1 < mu1_min
    use_e2 = (~use_e1) & (mu2 < mu2_min)
    mu1_b = jnp.where(use_e1, mu1_e1, jnp.where(use_e2, mu1_e2, mu1))
    mu2_b = jnp.where(use_e1, mu2_e1, jnp.where(use_e2, mu2_e2, mu2))
    R_b = jnp.where(use_e1, R_e1, jnp.where(use_e2, R_e2, R))

    # Calculate the in-precip. standard deviation of the hydrometeor in the 1st
    # and 2nd PDF components.
    sig1_b = sqrt_clipped(R_b * (1.0 + zeta)) * mu1_b
    sig2_b = sqrt_clipped(R_b) * mu2_b
    # Calculate the mean of the hydrometeor in the 1st and 2nd PDF components.
    hm1_b = jnp.maximum(mu1_b * precip_frac_1, hm_tol)
    hm2_b = jnp.maximum(mu2_b * precip_frac_2, hm_tol)

    # Else if there is precipitation found in the 1st PDF component but not
    # the 2nd PDF component.
    # Calculate the mean (in-precip.) of the hydrometeor in the 1st PDF
    # component.
    mu1_c1 = hmm / (a_s * fp1_s)
    # Calculate the in-precip. standard deviation of the hydrometeor in the 1st
    # PDF component.
    sig1_c1 = sqrt_clipped((hmp2 + hmm ** 2 - mixt_frac * precip_frac_1 * mu1_c1 ** 2) / (a_s * fp1_s))
    # Calculate the mean of the hydrometeor in the 1st PDF component.
    hm1_c1 = mu1_c1 * precip_frac_1

    # Else if there is precipitation found in the 2nd PDF component but not
    # the 1st PDF component.
    # Calculate the mean (in-precip.) of the hydrometeor in the 2nd PDF
    # component.
    mu2_c2 = hmm / (oma_s * fp2_s)
    # Calculate the in-precip. standard deviation of the hydrometeor in the 2nd
    # PDF component.
    sig2_c2 = sqrt_clipped((hmp2 + hmm ** 2 - oma * precip_frac_2 * mu2_c2 ** 2) / (oma_s * fp2_s))
    # Calculate the mean of the hydrometeor in the 2nd PDF component.
    hm2_c2 = mu2_c2 * precip_frac_2

    # Else there is no precipitation.
    z = jnp.zeros_like(hmm)
    mu_hm_1 = jnp.where(both, mu1_b, jnp.where(comp1, mu1_c1, z))
    mu_hm_2 = jnp.where(both, mu2_b, jnp.where(comp2, mu2_c2, z))
    sigma_hm_1 = jnp.where(both, sig1_b, jnp.where(comp1, sig1_c1, z))
    sigma_hm_2 = jnp.where(both, sig2_b, jnp.where(comp2, sig2_c2, z))
    hm_1 = jnp.where(both, hm1_b, jnp.where(comp1, hm1_c1, z))
    hm_2 = jnp.where(both, hm2_b, jnp.where(comp2, hm2_c2, z))

    mu1_safe = jnp.where(mu_hm_1 > 0.0, mu_hm_1, 1.0)
    mu2_safe = jnp.where(mu_hm_2 > 0.0, mu_hm_2, 1.0)
    s1r = jnp.where((both | comp1), sigma_hm_1 ** 2 / mu1_safe ** 2, 0.0)
    s2r = jnp.where(both, R_b, jnp.where(comp2, sigma_hm_2 ** 2 / mu2_safe ** 2, 0.0))

    return mu_hm_1, mu_hm_2, sigma_hm_1, sigma_hm_2, hm_1, hm_2, s1r, s2r


def hydrometp2_zt(hmm, precip_frac, ratio):
    """Overall variance <hm'^2> of a precipitating hydrometeor from its prescribed in-precip ratio
    (setup_clubb_pdf_params.F90:449): ((hmp2_ip_on_hmm2_ip + 1)/precip_frac − 1) · hm^2 where hm ≥ tol
    (the caller zeroes the sub-tolerance levels). Pure-jnp → differentiable."""
    pf = jnp.where(precip_frac > 0.0, precip_frac, 1.0)
    hydrometp2_zt = ((ratio + 1.0) / pf - 1.0) * hmm ** 2
    return hydrometp2_zt


def compute_mean_stdev(chi_1, chi_2, stdev_chi_1, stdev_chi_2,
                       stdev_eta_1, stdev_eta_2, Ncnm, Ncnp2_on_Ncnm2,
                       l_const_Nc_in_cloud, hydromets, thl_1, thl_2, mixt_frac,
                       precip_frac, precip_frac_1, precip_frac_2, precip_frac_tol,
                       omicron, zeta_vrnce_rat,
                       w_1=None, w_2=None, stdev_w_1=None, stdev_w_2=None):
    """Calculates the means and standard deviations (for each PDF component) of
    chi, eta, w, Ncn, and the precipitating hydrometeors.  For the
    precipitating hydrometeors, the component means and standard deviations
    are in-precip.

    Fortran declarations (compute_mean_stdev in src/CLUBB_core/setup_clubb_pdf_params.F90):

        !----------------------- Input Variables -----------------------
        integer, intent(in) :: &
          nzt,          & ! Number of model thermodynamic vertical grid levels
          pdf_dim,      & ! Number of PDF variables
          ngrdcol,      & ! Number of grid columns
          hydromet_dim    ! Number of hydrometeor species

        real( kind = core_rknd ), dimension(ngrdcol,nzt,hydromet_dim), intent(in) :: &
          hydromet,       & ! Mean of a hydrometeor (overall)         [hm units]
          hydrometp2_zt     ! Variance of a hydrometeor (overall)     [(hm units)^2]

        type (hm_metadata_type), intent(in) :: &
          hm_metadata

        real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
          Ncnm,          & ! Mean simplified cloud nuclei concentration   [num/kg]
          mixt_frac,     & ! Mixture fraction                             [-]
          precip_frac,   & ! Precipitation fraction (overall)             [-]
          precip_frac_1, & ! Precipitation fraction (1st PDF component)   [-]
          precip_frac_2    ! Precipitation fraction (2nd PDF component)   [-]

        real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
          w_1,         & ! Mean of w (1st PDF component)                  [m/s]
          w_2,         & ! Mean of w (2nd PDF component)                  [m/s]
          stdev_w_1,   & ! Variance of w (1st PDF component)              [m^2/s^2]
          stdev_w_2,   & ! Variance of w (2nd PDF component)              [m^2/s^2]
          chi_1,       & ! Mean of chi (1st PDF component)                [kg/kg]
          chi_2,       & ! Mean of chi (2nd PDF component)                [kg/kg]
          stdev_chi_1, & ! Standard deviation of chi (1st PDF component)  [kg/kg]
          stdev_chi_2, & ! Standard deviation of chi (2nd PDF component)  [kg/kg]
          stdev_eta_1, & ! Standard deviation of eta (1st PDF component)  [kg/kg]
          stdev_eta_2, & ! Standard deviation of eta (2nd PDF component)  [kg/kg]
          thl_1,       & ! Mean of thl (1st PDF component)                [K]
          thl_2          ! Mean of thl (2nd PDF component)                [K]

        real( kind = core_rknd ), dimension(ngrdcol), intent(in) :: &
          precip_frac_tol    ! Minimum precip. frac. when hydromet. are present  [-]

        real( kind = core_rknd ), dimension(ngrdcol), intent(in) :: &
          omicron,        & ! Relative width parameter, omicron = R / Rmax    [-]
          zeta_vrnce_rat    ! Width parameter for sigma_hm_1^2 / mu_hm_1^2    [-]

        logical, intent(in) :: &
          l_const_Nc_in_cloud ! Use a constant cloud droplet conc. within cloud (K&K)

        !----------------------- Output Variables -----------------------
        ! Note:  This code assumes to be these arrays in the same order as the
        ! correlation arrays, etc., which is determined by the iiPDF indices.
        ! The order should be as follows:  chi, eta, w, Ncn, <precip. hydrometeors>
        ! (indices increasing from left to right).
        real( kind = core_rknd ), dimension(ngrdcol,nzt,pdf_dim), intent(out) :: &
          mu_x_1,    & ! Mean array of PDF vars. (1st PDF component)    [units vary]
          mu_x_2,    & ! Mean array of PDF vars. (2nd PDF component)    [units vary]
          sigma_x_1, & ! Standard deviation array of PDF vars (comp. 1) [units vary]
          sigma_x_2    ! Standard deviation array of PDF vars (comp. 2) [units vary]

        real( kind = core_rknd ), dimension(ngrdcol,nzt,hydromet_dim), intent(out) :: &
          hm_1, & ! Mean of a precip. hydrometeor (1st PDF component)   [units vary]
          hm_2    ! Mean of a precip. hydrometeor (2nd PDF component)   [units vary]

        real( kind = core_rknd ), dimension(ngrdcol,nzt,pdf_dim), intent(out) :: &
          sigma_hm_1_sqd_on_mu_hm_1_sqd, & ! Ratio sigma_hm_1^2 / mu_hm_1^2      [-]
          sigma_hm_2_sqd_on_mu_hm_2_sqd    ! Ratio sigma_hm_2^2 / mu_hm_2^2      [-]
    """
    z = jnp.zeros_like(chi_1)
    w_1 = z if w_1 is None else w_1
    w_2 = z if w_2 is None else w_2
    stdev_w_1 = z if stdev_w_1 is None else stdev_w_1
    stdev_w_2 = z if stdev_w_2 is None else stdev_w_2

    # Initialize output variables.

    # Enter the PDF parameters.

    # Vertical velocity, w.
    # Extended liquid water mixing ratio, chi.
    # Coordinate orthogonal to chi, eta.
    # Simplified cloud nuclei concentration, Ncn.
    # Standard deviation of simplified cloud nuclei concentration, Ncn,
    # in PDF component 1.
    if l_const_Nc_in_cloud:
        # Ncn is constant in both PDF components.
        sig_Ncn = z
        s2m2_Ncn = z
    else:
        # Ncn varies in both PDF components.
        sig_Ncn = jnp.sqrt(Ncnp2_on_Ncnm2) * Ncnm
        s2m2_Ncn = jnp.broadcast_to(Ncnp2_on_Ncnm2, chi_1.shape)

    # Note:  This code assumes to be these arrays in the same order as the
    # correlation arrays, etc., which is determined by the iiPDF indices.
    # The order should be as follows:  chi, eta, w, Ncn, <precip. hydrometeors>
    # (indices increasing from left to right).
    mu1_cols = [chi_1, jnp.broadcast_to(z, chi_1.shape), w_1, Ncnm]
    mu2_cols = [chi_2, jnp.broadcast_to(z, chi_1.shape), w_2, Ncnm]
    sig1_cols = [stdev_chi_1,
                 stdev_eta_1, stdev_w_1, sig_Ncn]
    sig2_cols = [stdev_chi_2,
                 stdev_eta_2, stdev_w_2, sig_Ncn]
    s2m2_1_cols = [z, z, z, s2m2_Ncn]
    s2m2_2_cols = [z, z, z, s2m2_Ncn]
    hm1_list, hm2_list = [], []

    for hmm, hmp2_zt, ratio, hm_tol in hydromets:
        mu1, mu2, sig1, sig2, hm1, hm2, s1r, s2r = calc_comp_mu_sigma_hm(
            hmm, hmp2_zt, jnp.broadcast_to(ratio, chi_1.shape), mixt_frac,
            precip_frac, precip_frac_1, precip_frac_2, hm_tol, precip_frac_tol,
            thl_1, thl_2, omicron, zeta_vrnce_rat)
        mu1_cols.append(mu1); mu2_cols.append(mu2)
        sig1_cols.append(sig1); sig2_cols.append(sig2)
        s2m2_1_cols.append(s1r); s2m2_2_cols.append(s2r)
        hm1_list.append(hm1); hm2_list.append(hm2)

    stack = lambda cols: jnp.stack([jnp.broadcast_to(c, chi_1.shape) for c in cols], axis=-1)
    hm_1 = jnp.stack(hm1_list, axis=-1) if hm1_list else jnp.zeros(chi_1.shape + (0,))
    hm_2 = jnp.stack(hm2_list, axis=-1) if hm2_list else jnp.zeros(chi_1.shape + (0,))
    return (stack(mu1_cols), stack(mu2_cols), stack(sig1_cols), stack(sig2_cols),
            hm_1, hm_2, stack(s2m2_1_cols), stack(s2m2_2_cols))


def norm_transform_mean_stdev(mu_x_1, mu_x_2, sigma_x_1, sigma_x_2,
                              s2m2_1, s2m2_2, Ncnm, hm_1, hm_2, hydromet_tols,
                              l_const_Nc_in_cloud):
    """Transforms the means and the standard deviations of PDF variables that
    have assumed lognormal distributions -- which are precipitating
    hydrometeors (in precipitation) and N_cn -- to normal space for each PDF
    component.

    Where the in-precip. mean is below tolerance, Fortran sets the normal-space
    mean to -huge as an "absent" marker; the JAX function uses a finite floor
    (mean_L2N of |mu| floored to 1e-30), which keeps the vanishing weight in
    the integrals that use it while leaving gradients finite.

    Returns:
        (mu_x_1_n, mu_x_2_n, sigma_x_1_n, sigma_x_2_n)

    Fortran declarations (norm_transform_mean_stdev in src/CLUBB_core/setup_clubb_pdf_params.F90):

        !--------------------------- Input Variables ---------------------------
        integer, intent(in) :: &
          nzt,          & ! Number of model thermodynamic vertical grid levels
          pdf_dim,      & ! Number of variables in CLUBB's PDF
          ngrdcol,      & ! Number of grid columns
          hydromet_dim    ! Number of hydrometeor species

        type (hm_metadata_type), intent(in) :: &
          hm_metadata

        real( kind = core_rknd ), dimension(ngrdcol,nzt,hydromet_dim), intent(in) :: &
          hm_1, & ! Mean of a precip. hydrometeor (1st PDF component)   [units vary]
          hm_2    ! Mean of a precip. hydrometeor (2nd PDF component)   [units vary]

        real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
          Ncnm    ! Mean cloud nuclei concentration, < N_cn >               [num/kg]

        real( kind = core_rknd ), dimension(ngrdcol,nzt,pdf_dim), intent(in) :: &
          mu_x_1,    & ! Mean array of PDF vars. (1st PDF component)    [units vary]
          mu_x_2,    & ! Mean array of PDF vars. (2nd PDF component)    [units vary]
          sigma_x_1, & ! Standard deviation array of PDF vars (comp. 1) [units vary]
          sigma_x_2    ! Standard deviation array of PDF vars (comp. 2) [units vary]

        real( kind = core_rknd ), dimension(ngrdcol,nzt,pdf_dim), intent(in) :: &
          sigma2_on_mu2_ip_1, & ! Prescribed ratio array: sigma_hm_1^2/mu_hm_1^2 [-]
          sigma2_on_mu2_ip_2    ! Prescribed ratio array: sigma_hm_2^2/mu_hm_2^2 [-]

        logical, intent(in) :: &
          l_const_Nc_in_cloud ! Use a constant cloud droplet conc. within cloud (K&K)

        !--------------------------- Output Variables ---------------------------
        real( kind = core_rknd ), dimension(ngrdcol,nzt,pdf_dim), intent(out) :: &
          mu_x_1_n,    & ! Mean array (normal space): PDF vars. (comp. 1) [un. vary]
          mu_x_2_n,    & ! Mean array (normal space): PDF vars. (comp. 2) [un. vary]
          sigma_x_1_n, & ! Std. dev. array (normal space): PDF vars (comp. 1) [u.v.]
          sigma_x_2_n    ! Std. dev. array (normal space): PDF vars (comp. 2) [u.v.]
    """
    def _to_n(mu, s2m2):
        mu_floor = 1.0e-30
        mu_safe = jnp.maximum(jnp.abs(mu), mu_floor)
        return mean_L2N(mu_safe, s2m2), stdev_L2N(s2m2)

    mu_x_1_n = mu_x_1
    mu_x_2_n = mu_x_2
    sigma_x_1_n = sigma_x_1
    sigma_x_2_n = sigma_x_2

    # Ncn (single lognormal: components share the moments).
    mu_Ncn_1_n, sig_Ncn_1_n = _to_n(mu_x_1_n[..., IIPDF_NCN], s2m2_1[..., IIPDF_NCN])
    mu_Ncn_2_n, sig_Ncn_2_n = _to_n(mu_x_2_n[..., IIPDF_NCN], s2m2_2[..., IIPDF_NCN])
    if l_const_Nc_in_cloud:
        sig_Ncn_1_n = jnp.zeros_like(sig_Ncn_1_n)
        sig_Ncn_2_n = jnp.zeros_like(sig_Ncn_2_n)
    mu_x_1_n = mu_x_1_n.at[..., IIPDF_NCN].set(mu_Ncn_1_n)
    mu_x_2_n = mu_x_2_n.at[..., IIPDF_NCN].set(mu_Ncn_2_n)
    sigma_x_1_n = sigma_x_1_n.at[..., IIPDF_NCN].set(sig_Ncn_1_n)
    sigma_x_2_n = sigma_x_2_n.at[..., IIPDF_NCN].set(sig_Ncn_2_n)

    # Precipitating hydrometeors (pdf indices after Ncn).
    for j, _hm_tol in enumerate(hydromet_tols):
        ivar = IIPDF_NCN + 1 + j
        m1n, s1n = _to_n(mu_x_1_n[..., ivar], s2m2_1[..., ivar])
        m2n, s2n = _to_n(mu_x_2_n[..., ivar], s2m2_2[..., ivar])
        mu_x_1_n = mu_x_1_n.at[..., ivar].set(m1n)
        mu_x_2_n = mu_x_2_n.at[..., ivar].set(m2n)
        sigma_x_1_n = sigma_x_1_n.at[..., ivar].set(s1n)
        sigma_x_2_n = sigma_x_2_n.at[..., ivar].set(s2n)

    return mu_x_1_n, mu_x_2_n, sigma_x_1_n, sigma_x_2_n


def calc_corr_w_hm_n(wm, wphydrometp, mu_w_1, mu_w_2, mu_hm_1, mu_hm_2,
                     sigma_w_1, sigma_w_2, sigma_hm_1, sigma_hm_2, sigma_hm_1_n, sigma_hm_2_n,
                     mixt_frac, precip_frac_1, precip_frac_2, hm_tol):
    """Calculates the PDF component correlation (in-precip) between vertical
    velocity, w, and the natural logarithm of a hydrometeor, ln hm.  The
    overall covariance of w and hm, <w'hm'> can be written in terms of the PDF
    parameters.  When both w and hm vary in both PDF components, the equation
    is written as:

    <w'hm'> = mixt_frac * precip_frac_1
              * ( mu_w_1 - <w>
                  + corr_w_hm_1_n * sigma_w_1 * sigma_hm_1_n ) * mu_hm_1
              + ( 1 - mixt_frac ) * precip_frac_2
                * ( mu_w_2 - <w>
                    + corr_w_hm_2_n * sigma_w_2 * sigma_hm_2_n ) * mu_hm_2.

    The overall covariance is provided, so the component correlation is solved
    by setting corr_w_hm_1_n = corr_w_hm_2_n ( = corr_w_hm_n ).  The equation
    is:

    corr_w_hm_n
    = ( <w'hm'>
        - mixt_frac * precip_frac_1 * ( mu_w_1 - <w> ) * mu_hm_1
        - ( 1 - mixt_frac ) * precip_frac_2 * ( mu_w_2 - <w> ) * mu_hm_2 )
      / ( mixt_frac * precip_frac_1 * sigma_w_1 * sigma_hm_1_n * mu_hm_1
          + ( 1 - mixt_frac ) * precip_frac_2
            * sigma_w_2 * sigma_hm_2_n * mu_hm_2 );

    again, where corr_w_hm_1_n = corr_w_hm_2_n = corr_w_hm_n.  When either w
    or hm is constant in one PDF component, but both w and hm vary in the
    other PDF component, the equation for <w'hm'> is written as:

    <w'hm'> = mixt_frac * precip_frac_1
              * ( mu_w_1 - <w>
                  + corr_w_hm_1_n * sigma_w_1 * sigma_hm_1_n ) * mu_hm_1
              + ( 1 - mixt_frac ) * precip_frac_2
                * ( mu_w_2 - <w> ) * mu_hm_2.

    In the above equation, either w or hm (or both) is (are) constant in PDF
    component 2, but both w and hm vary in PDF component 1.  When both w and
    hm vary in PDF component 2, but at least one of w or hm is constant in PDF
    component 1, the equation is similar.  The above equation can be rewritten
    to solve for corr_w_hm_1_n, such that:

    corr_w_hm_1_n
    = ( <w'hm'>
        - mixt_frac * precip_frac_1 * ( mu_w_1 - <w> ) * mu_hm_1
        - ( 1 - mixt_frac ) * precip_frac_2 * ( mu_w_2 - <w> ) * mu_hm_2 )
      / ( mixt_frac * precip_frac_1 * sigma_w_1 * sigma_hm_1_n * mu_hm_1 ).

    Since either w or hm is constant in PDF component 2, corr_w_hm_2_n is
    undefined.  When both w and hm vary in PDF component 2, but at least one
    of w or hm is constant in PDF component 1, the equation is similar, but
    is in terms of corr_w_hm_2_n, while corr_w_hm_1_n is undefined.  When
    either w or hm is constant in both PDF components, the equation for
    <w'hm'> is:

    <w'hm'> = mixt_frac * precip_frac_1
              * ( mu_w_1 - <w> ) * mu_hm_1
              + ( 1 - mixt_frac ) * precip_frac_2
                * ( mu_w_2 - <w> ) * mu_hm_2.

    When this is the case, both corr_w_hm_1_n and corr_w_hm_2_n are undefined.

    Returns:
        (corr_1, corr_2)

    Fortran declarations (calc_corr_w_hm_n in src/CLUBB_core/setup_clubb_pdf_params.F90):

        !--------------------------- Input Variables ---------------------------
        real( kind = core_rknd ), intent(in) :: &
          wm,            & ! Mean vertical velocity (overall), <w>             [m/s]
          wphydrometp,   & ! Covariance of w and hm (overall), <w'hm'>  [m/s(hm un)]
          mu_w_1,        & ! Mean of w (1st PDF component)                     [m/s]
          mu_w_2,        & ! Mean of w (2nd PDF component)                     [m/s]
          mu_hm_1,       & ! Mean of hm (1st PDF component) in-precip (ip)   [hm un]
          mu_hm_2,       & ! Mean of hm (2nd PDF component) ip               [hm un]
          sigma_w_1,     & ! Standard deviation of w (1st PDF component)       [m/s]
          sigma_w_2,     & ! Standard deviation of w (2nd PDF component)       [m/s]
          sigma_hm_1,    & ! Standard deviation of hm (1st PDF component) ip [hm un]
          sigma_hm_2,    & ! Standard deviation of hm (2nd PDF component) ip [hm un]
          sigma_hm_1_n,  & ! Standard deviation of ln hm (1st PDF component) ip  [-]
          sigma_hm_2_n,  & ! Standard deviation of ln hm (2nd PDF component) ip  [-]
          mixt_frac,     & ! Mixture fraction                                    [-]
          precip_frac_1, & ! Precipitation fraction (1st PDF component)          [-]
          precip_frac_2    ! Precipitation fraction (2nd PDF component)          [-]

        real( kind = core_rknd ), intent(in) :: &
          hm_tol           ! Hydrometeor tolerance value                     [hm un]

        !--------------------------- Output Variables ---------------------------
        real( kind = core_rknd ), intent(out) :: &
          corr_w_hm_1_n, & ! Correlation of w and ln hm (1st PDF component) ip   [-]
          corr_w_hm_2_n    ! Correlation of w and ln hm (2nd PDF component) ip   [-]
    """
    w1, w2 = mixt_frac * precip_frac_1, (1.0 - mixt_frac) * precip_frac_2
    num = wphydrometp - w1 * (mu_w_1 - wm) * mu_hm_1 - w2 * (mu_w_2 - wm) * mu_hm_2
    den_both = w1 * sigma_w_1 * sigma_hm_1_n * mu_hm_1 + w2 * sigma_w_2 * sigma_hm_2_n * mu_hm_2
    den_1 = w1 * sigma_w_1 * sigma_hm_1_n * mu_hm_1
    den_2 = w2 * sigma_w_2 * sigma_hm_2_n * mu_hm_2

    def _clamp(x):
        return jnp.clip(x, -max_mag_correlation, max_mag_correlation)

    def _safe_div(n, d):
        return n / jnp.where(d != 0.0, d, 1.0)

    c1_vary = (sigma_w_1 > w_tol) & (sigma_hm_1 > hm_tol)
    c2_vary = (sigma_w_2 > w_tol) & (sigma_hm_2 > hm_tol)
    both = c1_vary & c2_vary

    cn = _clamp(_safe_div(num, den_both))
    c1v = _clamp(_safe_div(num, den_1))
    c2v = _clamp(_safe_div(num, den_2))

    corr_1 = jnp.where(both, cn, jnp.where(c1_vary, c1v, 0.0))
    corr_2 = jnp.where(both, cn, jnp.where(c1_vary, 0.0, jnp.where(c2_vary, c2v, 0.0)))
    return corr_1, corr_2


def _corr_cloud_below(rc_1, rc_2, corr_cloud, corr_below):
    """Per-component in-precip correlation by cloud presence: rc_i > rc_tol → cloud value, else below-cloud
    value (the shared body of component_corr_{x_hm,hmx_hmy,w_hm}_n_ip). rc_1/rc_2 are (ngrdcol,nzt)."""
    c1 = jnp.where(rc_1 > rc_tol, corr_cloud, corr_below)
    c2 = jnp.where(rc_2 > rc_tol, corr_cloud, corr_below)
    return c1, c2


def component_corr_w_hm_n_ip(corr_w_hm_1_n_in, rc_1, corr_w_hm_2_n_in, rc_2,
                             corr_w_hm_n_NL_cloud, corr_w_hm_n_NL_below, l_calc_w_corr):
    """Calculates the in-precip correlation of w and the natural logarithm of a
    hydrometeor species within the ith PDF component.

    Fortran declarations (component_corr_w_hm_n_ip in src/CLUBB_core/setup_clubb_pdf_params.F90):

        !--------------------------- Input Variables ---------------------------
        integer, intent(in) :: &
          nzt,     & ! Number of model thermodynamic vertical grid levels
          ngrdcol    ! Number of grid columns

        real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
          corr_w_hm_1_n_in, & ! Correlation of w and ln hm (1st PDF comp.) ip    [-]
          corr_w_hm_2_n_in, & ! Correlation of w and ln hm (2nd PDF comp.) ip    [-]
          rc_1,             & ! Mean cloud water mix. ratio (1st PDF comp.)  [kg/kg]
          rc_2                ! Mean cloud water mix. ratio (2nd PDF comp.)  [kg/kg]

        real( kind = core_rknd ), intent(in) :: &
          corr_w_hm_n_NL_cloud, & ! Corr. of w & ln hm (ith PDF comp.) ip; cloud [-]
          corr_w_hm_n_NL_below    ! Corr. of w & ln hm (ith PDF comp.) ip; clear [-]

        logical, intent(in) :: &
          l_calc_w_corr ! Calculate the correlations between w and the hydrometeors

        !--------------------------- Output Variables ---------------------------
        real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(out) :: &
          corr_w_hm_1_n, & ! Correlation of w and ln hm (1st PDF component) ip   [-]
          corr_w_hm_2_n    ! Correlation of w and ln hm (2nd PDF component) ip   [-]
    """
    if l_calc_w_corr:
        return corr_w_hm_1_n_in, corr_w_hm_2_n_in
    return _corr_cloud_below(rc_1, rc_2, corr_w_hm_n_NL_cloud, corr_w_hm_n_NL_below)


def component_corr_x_hm_n_ip(rc_1, rc_2, corr_x_hm_n_NL_cloud, corr_x_hm_n_NL_below):
    """Calculates the in-precip correlation of x and a hydrometeor species
    within the ith PDF component.  Here, x is a variable with a normally
    distributed individual marginal PDF, such as chi or eta.

    Fortran declarations (component_corr_x_hm_n_ip in src/CLUBB_core/setup_clubb_pdf_params.F90):

        !--------------------------- Input Variables ---------------------------
        integer, intent(in) :: &
          nzt,     & ! Number of model thermodynamic vertical grid levels
          ngrdcol    ! Number of grid columns

        real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
          rc_1,         & ! Mean cloud water mixing ratio (1st PDF comp.)   [kg/kg]
          rc_2            ! Mean cloud water mixing ratio (2nd PDF comp.)   [kg/kg]

        real( kind = core_rknd ), intent(in) :: &
          corr_x_hm_n_NL_cloud, & ! Corr. of x and ln hm (ith PDF comp.) ip     [-]
          corr_x_hm_n_NL_below    ! Corr. of x and ln hm (ith PDF comp.) ip     [-]

        !--------------------------- Output Variables ---------------------------
        real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(out) :: &
          corr_x_hm_1_n, & ! Correlation of x and ln hm (1st PDF component) ip  [-]
          corr_x_hm_2_n    ! Correlation of x and ln hm (2nd PDF component) ip  [-]
    """
    return _corr_cloud_below(rc_1, rc_2, corr_x_hm_n_NL_cloud, corr_x_hm_n_NL_below)


def component_corr_hmx_hmy_n_ip(rc_1, rc_2, corr_hmx_hmy_n_LL_cloud, corr_hmx_hmy_n_LL_below):
    """Calculates the in-precip correlation of the natural logarithms of
    hydrometeor x and hydrometeor y within the ith PDF component.

    Fortran declarations (component_corr_hmx_hmy_n_ip in src/CLUBB_core/setup_clubb_pdf_params.F90):

        !--------------------------- Input Variables ---------------------------
        integer, intent(in) :: &
          nzt,     & ! Number of model thermodynamic vertical grid levels
          ngrdcol    ! Number of grid columns

        real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
          rc_1,         & ! Mean cloud water mixing ratio (ith PDF comp.) [kg/kg]
          rc_2            ! Mean cloud water mixing ratio (ith PDF comp.) [kg/kg]

        real( kind = core_rknd ), intent(in) :: &
          corr_hmx_hmy_n_LL_cloud, & ! Corr.: ln hmx & ln hmy (ith PDF comp.) ip [-]
          corr_hmx_hmy_n_LL_below    ! Corr.: ln hmx & ln hmy (ith PDF comp.) ip [-]

        !--------------------------- Output Variables ---------------------------
        real( kind = core_rknd ), dimension(ngrdcol,nzt) :: &
          corr_hmx_hmy_1_n, & ! Corr. of ln hmx & ln hmy (ith PDF comp.) ip      [-]
          corr_hmx_hmy_2_n    ! Corr. of ln hmx & ln hmy (ith PDF comp.) ip      [-]
    """
    return _corr_cloud_below(rc_1, rc_2, corr_hmx_hmy_n_LL_cloud, corr_hmx_hmy_n_LL_below)


def component_corr_eta_hm_n_ip(corr_chi_eta_1, corr_chi_hm_n_1, corr_chi_eta_2, corr_chi_hm_n_2):
    """Estimates the correlation of eta and the natural logarithm of a
    hydrometeor species using the correlation of chi and eta and the
    correlation of chi and the natural logarithm of the hydrometeor.  This
    facilities the Cholesky decomposability of the correlation array that will
    inevitably be decomposed for SILHS purposes. Without this estimation, we
    have found that the resulting correlation matrix cannot be decomposed.

    Returns:
        (corr_eta_hm_n_1, corr_eta_hm_n_2)

    Fortran declarations (component_corr_eta_hm_n_ip in src/CLUBB_core/setup_clubb_pdf_params.F90):

        !--------------------------- Input Variables ---------------------------
        integer, intent(in) :: &
          nzt,     & ! Number of model thermodynamic vertical grid levels
          ngrdcol    ! Number of grid columns

        real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
          corr_chi_eta_1,  & ! Component correlation of chi and eta              [-]
          corr_chi_eta_2,  & ! Component correlation of chi and eta              [-]
          corr_chi_hm_n_1, & ! Component correlation of chi and ln hm            [-]
          corr_chi_hm_n_2    ! Component correlation of chi and ln hm            [-]

        !--------------------------- Output Variables ---------------------------
        real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(out) :: &
          corr_eta_hm_n_1, & ! Component correlation of eta and ln hm            [-]
          corr_eta_hm_n_2    ! Component correlation of eta and ln hm            [-]
    """
    corr_eta_hm_n_1 = corr_chi_eta_1 * corr_chi_hm_n_1
    corr_eta_hm_n_2 = corr_chi_eta_2 * corr_chi_hm_n_2
    return corr_eta_hm_n_1, corr_eta_hm_n_2


# iiPDF_type enumeration (model_flags.F90:31-37) — imported from its Fortran-home model_flags.py (the PDF types
# whose ADG standards fix corr(w,x)=0).
from clubb_jax.src.CLUBB_core.model_flags import (
    iiPDF_ADG1 as IIPDF_TYPE_ADG1, iiPDF_ADG2 as IIPDF_TYPE_ADG2, iiPDF_new_hybrid as IIPDF_TYPE_NEW_HYBRID,
)


def component_corr_w_x(rc_1, rc_2, corr_w_x_NN_cloud, corr_w_x_NN_below,
                       iiPDF_type, l_follow_ADG1_PDF_standards):
    """Calculates the correlation of w and x within the ith PDF component.
    Here, x is a variable with a normally distributed individual marginal PDF,
    such as chi or eta.

    Fortran declarations (component_corr_w_x in src/CLUBB_core/setup_clubb_pdf_params.F90):

        !--------------------------- Input Variables ---------------------------
        integer, intent(in) :: &
          nzt,     & ! Number of model thermodynamic vertical grid levels
          ngrdcol    ! Number of grid columns

        real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
          rc_1,           & ! Mean cloud water mixing ratio (1st PDF comp.)  [kg/kg]
          rc_2              ! Mean cloud water mixing ratio (2nd PDF comp.)  [kg/kg]

        real( kind = core_rknd ), intent(in) :: &
          corr_w_x_NN_cloud, & ! Corr. of w and x (ith PDF comp.); cloudy levs [-]
          corr_w_x_NN_below    ! Corr. of w and x (ith PDF comp.); clear levs  [-]

        integer, intent(in) :: &
          iiPDF_type    ! Selected option for the two-component normal (double
                        ! Gaussian) PDF type to use for the w, rt, and theta-l (or
                        ! w, chi, and eta) portion of CLUBB's multivariate,
                        ! two-component PDF.

        !--------------------------- Output Variables ---------------------------
        real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(out) :: &
          corr_w_x_1, & ! Correlation of w and x (1st PDF component)  [-]
          corr_w_x_2    ! Correlation of w and x (2nd PDF component)  [-]
    """
    if l_follow_ADG1_PDF_standards and iiPDF_type in (IIPDF_TYPE_ADG1, IIPDF_TYPE_ADG2, IIPDF_TYPE_NEW_HYBRID):
        return jnp.zeros_like(rc_1), jnp.zeros_like(rc_1)
    return _corr_cloud_below(rc_1, rc_2, corr_w_x_NN_cloud, corr_w_x_NN_below)


def component_corr_chi_eta(rc_1, rc_2, corr_chi_eta_NN_cloud, corr_chi_eta_NN_below,
                           l_limit_corr_chi_eta):
    """Calculates the correlation of chi (old s) and eta (old t) within the
    ith PDF component.

    Returns:
        (c1, c2)

    Fortran declarations (component_corr_chi_eta in src/CLUBB_core/setup_clubb_pdf_params.F90):

        !--------------------------- Input Variables ---------------------------
        integer, intent(in) :: &
          nzt,     & ! Number of model thermodynamic vertical grid levels
          ngrdcol    ! Number of grid columns

        real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
          rc_1,               & ! Mean cloud water mix. rat. (1st PDF comp.) [kg/kg]&
          rc_2                  ! Mean cloud water mix. rat. (2nd PDF comp.) [kg/kg]

        real( kind = core_rknd ), intent(in) :: &
          corr_chi_eta_NN_cloud, & ! Corr. of chi & eta (ith PDF comp.); cloudy  [-]
          corr_chi_eta_NN_below    ! Corr. of chi & eta (ith PDF comp.); clear   [-]

        logical, intent(in) :: &
          l_limit_corr_chi_eta    ! We must limit the correlation of chi and eta if
                                  ! we are to take the Cholesky decomposition of the
                                  ! resulting correlation matrix. This is because a
                                  ! perfect correlation of chi and eta was found to
                                  ! be unrealizable.

        !--------------------------- Output Variables ---------------------------
        real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(out) :: &
          corr_chi_eta_1, & ! Correlation of chi and eta (1st PDF component)     [-]
          corr_chi_eta_2    ! Correlation of chi and eta (2nd PDF component)     [-]
    """
    c1, c2 = _corr_cloud_below(rc_1, rc_2, corr_chi_eta_NN_cloud, corr_chi_eta_NN_below)
    if l_limit_corr_chi_eta:
        c1 = jnp.clip(c1, -max_mag_correlation, max_mag_correlation)
        c2 = jnp.clip(c2, -max_mag_correlation, max_mag_correlation)
    return c1, c2


def comp_corr_norm(mu_x_1, mu_x_2, sigma_x_1, sigma_x_2, sigma_x_1_n, sigma_x_2_n,
                   wm_zt, rc_1, rc_2, mixt_frac, precip_frac_1, precip_frac_2,
                   wpNcnp_zt, wphydrometp_zt, corr_array_n_cloud, corr_array_n_below,
                   iiPDF_chi, iiPDF_eta, iiPDF_w, iiPDF_Ncn, pdf2hydromet, hydromet_tol, Ncn_tol_val,
                   iiPDF_type, l_calc_w_corr, l_fix_w_chi_eta_correlations,
                   pdf_params_corr=None):
    """Assemble the normal-space PDF correlation arrays.

    Builds the lower-triangular (then symmetrized) (ngrdcol, nzt, pdf_dim, pdf_dim) correlation arrays for the
    two PDF components from the component_corr_* routines plus calc_corr_w_hm_n (for the w-correlations when
    l_calc_w_corr). The prescribed-array index layout is chi, eta, w, Ncn, <hydrometeors> (iiPDF indices).

    Arrays: mu_x_*/sigma_x_* are (ngrdcol, nzt, pdf_dim); wm_zt/rc_*/mixt_frac/precip_frac_*/wpNcnp_zt are
    (ngrdcol, nzt); wphydrometp_zt is (ngrdcol, nzt, hydromet_dim); corr_array_n_cloud/below are
    (pdf_dim, pdf_dim). pdf2hydromet maps a pdf index to its hydromet index. pdf_params_corr (with keys
    corr_chi_eta_1/2, corr_w_chi_1/2) is required only for the l_fix_w_chi_eta_correlations=False path.

    Returns (corr_array_1_n, corr_array_2_n), each (ngrdcol, nzt, pdf_dim, pdf_dim), symmetric, unit diagonal.

    Faithfulness note: the l_fix_w_chi_eta_correlations=False ("preferred") branch reproduces the Fortran's
    eta–w block exactly, including its quirk (F90:1560) of writing corr_w_chi into the (w, chi) slot a second
    time rather than (w, eta) — so (w, eta) is left at 0 on that path, matching the oracle.

    Fortran declarations (comp_corr_norm in src/CLUBB_core/setup_clubb_pdf_params.F90):

        !------------------------ Input Variables ------------------------
        integer, intent(in) :: &
          nzt,          & ! Number of thermodynamic vertical levels
          pdf_dim,      & ! Number of variables in the corr/mean/stdev arrays
          ngrdcol,      & ! Number of grid columns
          hydromet_dim    ! Number of hydrometeor species

        type (hm_metadata_type), intent(in) :: &
          hm_metadata

        real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
          wm_zt,         & ! Mean vertical velocity, <w>, on thermo. levels    [m/s]
          rc_1,          & ! Mean of r_c (1st PDF component)                 [kg/kg]
          rc_2,          & ! Mean of r_c (2nd PDF component)                 [kg/kg]
          mixt_frac,     & ! Mixture fraction                                    [-]
          precip_frac_1, & ! Precipitation fraction (1st PDF component)          [-]
          precip_frac_2, & ! Precipitation fraction (2nd PDF component)          [-]
          wpNcnp_zt        ! Covariance of w and N_cn on t-levs.      [(m/s) num/kg]

        real( kind = core_rknd ), dimension(ngrdcol,nzt,hydromet_dim), intent(in) :: &
          wphydrometp_zt    ! Covariance of w and hm interp. to t-levs.  [(m/s)u.v.]

        real( kind = core_rknd ), dimension(ngrdcol,nzt,pdf_dim), intent(in) :: &
          mu_x_1,      & ! Mean of x array (1st PDF component)          [units vary]
          mu_x_2,      & ! Mean of x array (2nd PDF component)          [units vary]
          sigma_x_1,   & ! Standard deviation of x array (1st PDF comp.)  [un. vary]
          sigma_x_2,   & ! Standard deviation of x array (2nd PDF comp.)  [un. vary]
          sigma_x_1_n, & ! Std. dev. array (normal space): PDF vars (comp. 1) [u.v.]
          sigma_x_2_n    ! Std. dev. array (normal space): PDF vars (comp. 2) [u.v.]

        real( kind = core_rknd ), dimension(pdf_dim,pdf_dim), &
        intent(in) :: &
          corr_array_n_cloud, & ! Prescribed correlation array in cloud        [-]
          corr_array_n_below    ! Prescribed correlation array below cloud     [-]

        type(pdf_parameter), intent(in) :: &
          pdf_params    ! PDF parameters                               [units vary]

        integer, intent(in) :: &
          iiPDF_type    ! Selected option for the two-component normal (double
                        ! Gaussian) PDF type to use for the w, rt, and theta-l (or
                        ! w, chi, and eta) portion of CLUBB's multivariate,
                        ! two-component PDF.

        logical, intent(in) :: &
          l_calc_w_corr, &                ! Calculate the correlations between w and the hydrometeors
          l_fix_w_chi_eta_correlations    ! Use a fixed correlation for s and t Mellor(chi/eta)

        !------------------------ Output Variables ------------------------
        real( kind = core_rknd ), dimension(ngrdcol,nzt,pdf_dim,pdf_dim), intent(out) :: &
          corr_array_1_n, & ! Corr. array (normal space) of PDF vars. (comp. 1)  [-]
          corr_array_2_n    ! Corr. array (normal space) of PDF vars. (comp. 2)  [-]
    """
    cc = corr_array_n_cloud; cb = corr_array_n_below

    ng, nzt, pdf_dim = mu_x_1.shape
    ones = jnp.ones((ng, nzt))
    hm_indices = list(range(iiPDF_Ncn + 1, pdf_dim))   # hydrometeor pdf indices

    A1 = jnp.zeros((ng, nzt, pdf_dim, pdf_dim))
    A2 = jnp.zeros((ng, nzt, pdf_dim, pdf_dim))
    idx = jnp.arange(pdf_dim)
    A1 = A1.at[:, :, idx, idx].set(1.0)
    A2 = A2.at[:, :, idx, idx].set(1.0)
    # Normal space correlations
    # Initialize corr_w_hm_1_n and corr_w_hm_2_n arrays to 0.
    # Set ones_vector to a vector of 1s.
    # Initialize the normal space correlation arrays
    # The corr_arrays are assumed to be lower triangular matrices
    # Set diagonal elements to 1
    # This code assumes the following order in the prescribed correlation
    # arrays (iiPDF indices):
    # chi, eta, w, Ncn, <hydrometeors> (indices increasing from left to right)

    def _set(A, r, c, val):
        return A.at[:, :, r, c].set(val)

    # Calculate normal space correlations involving w by first calculating total
    # covariances involving w (<w'Ncn'>, etc.) using the down-gradient
    # approximation.
    corr_w_Ncn_1 = jnp.zeros((ng, nzt)); corr_w_Ncn_2 = jnp.zeros((ng, nzt))
    corr_w_hm_1 = {j: jnp.zeros((ng, nzt)) for j in hm_indices}
    corr_w_hm_2 = {j: jnp.zeros((ng, nzt)) for j in hm_indices}
    if l_calc_w_corr:
        # Calculate the correlation of w and ln Ncn in each PDF component.
        # The subroutine calc_corr_w_hm_n can be used to do this as long as a
        # value of 1 is sent in for precip_frac_1 and precip_frac_2.
        corr_w_Ncn_1, corr_w_Ncn_2 = calc_corr_w_hm_n(
            wm_zt, wpNcnp_zt, mu_x_1[:, :, iiPDF_w], mu_x_2[:, :, iiPDF_w],
            mu_x_1[:, :, iiPDF_Ncn], mu_x_2[:, :, iiPDF_Ncn],
            sigma_x_1[:, :, iiPDF_w], sigma_x_2[:, :, iiPDF_w],
            sigma_x_1[:, :, iiPDF_Ncn], sigma_x_2[:, :, iiPDF_Ncn],
            sigma_x_1_n[:, :, iiPDF_Ncn], sigma_x_2_n[:, :, iiPDF_Ncn],
            mixt_frac, ones, ones, Ncn_tol_val)
        for j in hm_indices:
            # Calculate the correlation of w and the natural logarithm of the
            # hydrometeor for each PDF component and each hydrometeor type.
            hm_idx = int(pdf2hydromet[j])
            corr_w_hm_1[j], corr_w_hm_2[j] = calc_corr_w_hm_n(
                wm_zt, wphydrometp_zt[:, :, hm_idx], mu_x_1[:, :, iiPDF_w], mu_x_2[:, :, iiPDF_w],
                mu_x_1[:, :, j], mu_x_2[:, :, j], sigma_x_1[:, :, iiPDF_w], sigma_x_2[:, :, iiPDF_w],
                sigma_x_1[:, :, j], sigma_x_2[:, :, j], sigma_x_1_n[:, :, j], sigma_x_2_n[:, :, j],
                mixt_frac, precip_frac_1, precip_frac_2, hydromet_tol[hm_idx])

    # In order to decompose the normal space correlation matrix,
    # we must not have a perfect correlation of chi and
    # eta. Thus, we impose a limitation.
    if l_fix_w_chi_eta_correlations:
        # Correlation of chi (old s) and eta (old t)
        l_limit_corr_chi_eta = True  # Flag to limit the correlation of chi and eta [-]
        c1, c2 = component_corr_chi_eta(rc_1, rc_2, cc[iiPDF_eta, iiPDF_chi],
                                        cb[iiPDF_eta, iiPDF_chi], l_limit_corr_chi_eta)
    else:
        # Preferred, more accurate version.
        c1 = pdf_params_corr['corr_chi_eta_1']
        c2 = pdf_params_corr['corr_chi_eta_2']
    A1 = _set(A1, iiPDF_eta, iiPDF_chi, c1); A2 = _set(A2, iiPDF_eta, iiPDF_chi, c2)

    if l_fix_w_chi_eta_correlations:
        # Correlation of chi (old s) and w
        c1, c2 = component_corr_w_x(rc_1, rc_2, cc[iiPDF_w, iiPDF_chi], cb[iiPDF_w, iiPDF_chi],
                                    iiPDF_type, l_follow_ADG1_PDF_standards)
    else:
        # Preferred, more accurate version.
        c1 = pdf_params_corr['corr_w_chi_1']; c2 = pdf_params_corr['corr_w_chi_2']
    A1 = _set(A1, iiPDF_w, iiPDF_chi, c1); A2 = _set(A2, iiPDF_w, iiPDF_chi, c2)

    # Correlation of chi (old s) and ln Ncn, corr_array_n_cloud used twice because
    # Ncn is an inherently in-cloud property.
    c1, c2 = component_corr_x_hm_n_ip(rc_1, rc_2, cc[iiPDF_Ncn, iiPDF_chi], cc[iiPDF_Ncn, iiPDF_chi])
    A1 = _set(A1, iiPDF_Ncn, iiPDF_chi, c1); A2 = _set(A2, iiPDF_Ncn, iiPDF_chi, c2)

    # Correlation of chi (old s) and the natural logarithm of the hydrometeors
    for j in hm_indices:
        c1, c2 = component_corr_x_hm_n_ip(rc_1, rc_2, cc[j, iiPDF_chi], cb[j, iiPDF_chi])
        A1 = _set(A1, j, iiPDF_chi, c1); A2 = _set(A2, j, iiPDF_chi, c2)

    # Correlation of eta (old t) and w
    if l_fix_w_chi_eta_correlations:
        # Correlation of chi (old s) and w
        c1, c2 = component_corr_w_x(rc_1, rc_2, cc[iiPDF_w, iiPDF_eta], cb[iiPDF_w, iiPDF_eta],
                                    iiPDF_type, l_follow_ADG1_PDF_standards)
        A1 = _set(A1, iiPDF_w, iiPDF_eta, c1); A2 = _set(A2, iiPDF_w, iiPDF_eta, c2)
    else:
        # Faithful to F90:1560 quirk: re-writes (w, chi), leaving (w, eta) at 0.
        # Preferred, more accurate version.
        c1 = pdf_params_corr['corr_w_chi_1']; c2 = pdf_params_corr['corr_w_chi_2']
        A1 = _set(A1, iiPDF_w, iiPDF_chi, c1); A2 = _set(A2, iiPDF_w, iiPDF_chi, c2)

    # Correlation of eta (old t) and ln Ncn, corr_array_n_cloud used twice because
    # Ncn is an inherently in-cloud property.
    c1, c2 = component_corr_x_hm_n_ip(rc_1, rc_2, cc[iiPDF_Ncn, iiPDF_eta], cc[iiPDF_Ncn, iiPDF_eta])
    A1 = _set(A1, iiPDF_Ncn, iiPDF_eta, c1); A2 = _set(A2, iiPDF_Ncn, iiPDF_eta, c2)

    # Correlation of eta (old t) and the natural logarithm of the hydrometeors
    for j in hm_indices:
        c1, c2 = component_corr_eta_hm_n_ip(A1[:, :, iiPDF_eta, iiPDF_chi], A1[:, :, j, iiPDF_chi],
                                            A2[:, :, iiPDF_eta, iiPDF_chi], A2[:, :, j, iiPDF_chi])
        A1 = _set(A1, j, iiPDF_eta, c1); A2 = _set(A2, j, iiPDF_eta, c2)

    # Correlation of w and ln Ncn
    c1, c2 = component_corr_w_hm_n_ip(corr_w_Ncn_1, rc_1, corr_w_Ncn_2, rc_2,
                                      cc[iiPDF_Ncn, iiPDF_w], cb[iiPDF_Ncn, iiPDF_w], l_calc_w_corr)
    A1 = _set(A1, iiPDF_Ncn, iiPDF_w, c1); A2 = _set(A2, iiPDF_Ncn, iiPDF_w, c2)

    # Correlation of w and the natural logarithm of the hydrometeors
    for j in hm_indices:
        c1, c2 = component_corr_w_hm_n_ip(corr_w_hm_1[j], rc_1, corr_w_hm_2[j], rc_2,
                                          cc[j, iiPDF_w], cb[j, iiPDF_w], l_calc_w_corr)
        A1 = _set(A1, j, iiPDF_w, c1); A2 = _set(A2, j, iiPDF_w, c2)

    # Correlation of ln Ncn and the natural logarithm of the hydrometeors
    for j in hm_indices:
        c1, c2 = component_corr_hmx_hmy_n_ip(rc_1, rc_2, cc[j, iiPDF_Ncn], cb[j, iiPDF_Ncn])
        A1 = _set(A1, j, iiPDF_Ncn, c1); A2 = _set(A2, j, iiPDF_Ncn, c2)

    # Correlation of the natural logarithm of two hydrometeors
    for ii in range(iiPDF_Ncn + 1, pdf_dim - 1):
        for jj in range(ii + 1, pdf_dim):
            c1, c2 = component_corr_hmx_hmy_n_ip(rc_1, rc_2, cc[jj, ii], cb[jj, ii])
            A1 = _set(A1, jj, ii, c1); A2 = _set(A2, jj, ii, c2)

    # For ease of use later in the code, we make the correlation arrays
    # symmetrical
    eye = jnp.zeros((ng, nzt, pdf_dim, pdf_dim)).at[:, :, idx, idx].set(1.0)
    A1 = A1 + jnp.swapaxes(A1, -1, -2) - eye
    A2 = A2 + jnp.swapaxes(A2, -1, -2) - eye
    return A1, A2


def denorm_transform_corr(sigma_x_1_n, sigma_x_2_n, sigma2_on_mu2_ip_1, sigma2_on_mu2_ip_2,
                          corr_array_1_n, corr_array_2_n,
                          iiPDF_chi, iiPDF_eta, iiPDF_w, iiPDF_Ncn):
    """Calculates the true or "real-space" correlations between PDF variables,
    where at least one of the variables that is part of a correlation has an
    assumed lognormal distribution -- which are the precipitating hydrometeors
    (in precipitation) and N_cn.

    As in Fortran, the component-2 transforms involving Ncn reuse the
    component-1 Ncn variance ratio sigma2_on_mu2_ip_1[Ncn]; Ncn is inherently
    in-cloud, so its ratio is shared across components.

    Returns:
        (A1, A2)

    Fortran declarations (denorm_transform_corr in src/CLUBB_core/setup_clubb_pdf_params.F90):

        !--------------------------- Input Variables ---------------------------
        integer, intent(in) :: &
          nzt,     & ! Number of thermodynamic vertical levels
          pdf_dim, & ! Number of variables
          ngrdcol    ! Number of grid columns

        type (hm_metadata_type), intent(in) :: &
          hm_metadata

        real( kind = core_rknd ), dimension(ngrdcol,nzt,pdf_dim), intent(in) :: &
          sigma_x_1_n, & ! Std. dev. array (normal space): PDF vars (comp. 1) [u.v.]
          sigma_x_2_n    ! Std. dev. array (normal space): PDF vars (comp. 2) [u.v.]

        real( kind = core_rknd ), dimension(ngrdcol,nzt,pdf_dim), intent(in) :: &
          sigma2_on_mu2_ip_1, & ! Ratio array sigma_hm_1^2/mu_hm_1^2             [-]
          sigma2_on_mu2_ip_2    ! Ratio array sigma_hm_2^2/mu_hm_2^2             [-]

        real( kind = core_rknd ), dimension(ngrdcol,nzt,pdf_dim,pdf_dim), &
        intent(in) :: &
          corr_array_1_n, & ! Corr. array (normal space) of PDF vars. (comp. 1)  [-]
          corr_array_2_n    ! Corr. array (normal space) of PDF vars. (comp. 2)  [-]

        !--------------------------- Output Variables ---------------------------
        real( kind = core_rknd ), dimension(ngrdcol,nzt,pdf_dim,pdf_dim), &
        intent(out) ::  &
          corr_array_1, & ! Correlation array of PDF vars. (comp. 1)             [-]
          corr_array_2    ! Correlation array of PDF vars. (comp. 2)             [-]
    """
    s1n = sigma_x_1_n; s2n = sigma_x_2_n
    r1 = sigma2_on_mu2_ip_1; r2 = sigma2_on_mu2_ip_2
    Cn1 = corr_array_1_n; Cn2 = corr_array_2_n

    ng, nzt, pdf_dim, _ = Cn1.shape
    idx = jnp.arange(pdf_dim)
    A1 = jnp.zeros((ng, nzt, pdf_dim, pdf_dim)).at[:, :, idx, idx].set(1.0)
    A2 = jnp.zeros((ng, nzt, pdf_dim, pdf_dim)).at[:, :, idx, idx].set(1.0)
    hm_indices = list(range(iiPDF_Ncn + 1, pdf_dim))

    def _s(A, r, c, val):
        return A.at[:, :, r, c].set(val)

    # Initialize diagonal elements to one
    # The correlations in each PDF component between two of w, chi (old s), and
    # eta (old t) do not need to be transformed to standard space, since w, chi,
    # and eta follow assumed normal distributions in each PDF component.  The
    # normal space correlations between any two of these variables are the same
    # as the actual correlations.
    for (rr, cc_) in ((iiPDF_eta, iiPDF_chi), (iiPDF_w, iiPDF_chi), (iiPDF_w, iiPDF_eta)):
        A1 = _s(A1, rr, cc_, Cn1[:, :, rr, cc_]); A2 = _s(A2, rr, cc_, Cn2[:, :, rr, cc_])

    # Calculate the true correlation of variables that have an assumed normal
    # distribution and variables that have an assumed lognormal distribution
    # for the ith PDF component, given their normal space correlation and the
    # normal space standard deviation of the variable with the assumed
    # lognormal distribution.
    # Transform the correlations between chi/eta/w and N_cn to standard space.
    for x in (iiPDF_chi, iiPDF_eta, iiPDF_w):
        A1 = _s(A1, iiPDF_Ncn, x, corr_NN2NL(Cn1[:, :, iiPDF_Ncn, x], s1n[:, :, iiPDF_Ncn], r1[:, :, iiPDF_Ncn]))
        A2 = _s(A2, iiPDF_Ncn, x, corr_NN2NL(Cn2[:, :, iiPDF_Ncn, x], s2n[:, :, iiPDF_Ncn], r1[:, :, iiPDF_Ncn]))

    # Transform the correlations (in-precip) between chi/eta/w and the
    # precipitating hydrometeors to standard space.
    for x in (iiPDF_chi, iiPDF_eta, iiPDF_w):
        for j in hm_indices:
            # Transform the correlation (in-precip) between w, chi, or eta and a
            # precipitating hydrometeor, hm, to standard space in PDF component 1.
            A1 = _s(A1, j, x, corr_NN2NL(Cn1[:, :, j, x], s1n[:, :, j], r1[:, :, j]))
            # Transform the correlation (in-precip) between w, chi, or eta and a
            # precipitating hydrometeor, hm, to standard space in PDF component 2.
            A2 = _s(A2, j, x, corr_NN2NL(Cn2[:, :, j, x], s2n[:, :, j], r2[:, :, j]))

    # Calculate the true correlation of two variables that both have an
    # assumed lognormal distribution for the ith PDF component, given their
    # normal space correlation and both of their normal space standard
    # deviations.
    # Transform the correlations (in-precip) between N_cn and the precipitating
    # hydrometeors to standard space.
    for j in hm_indices:
        # Transform the correlation (in-precip) between N_cn and a precipitating
        # hydrometeor, hm, to standard space in PDF component 1.
        A1 = _s(A1, j, iiPDF_Ncn, corr_NN2LL(Cn1[:, :, j, iiPDF_Ncn], s1n[:, :, iiPDF_Ncn], s1n[:, :, j],
                                             r1[:, :, iiPDF_Ncn], r1[:, :, j]))
        # Transform the correlation (in-precip) between N_cn and a precipitating
        # hydrometeor, hm, to standard space in PDF component 2.
        A2 = _s(A2, j, iiPDF_Ncn, corr_NN2LL(Cn2[:, :, j, iiPDF_Ncn], s2n[:, :, iiPDF_Ncn], s2n[:, :, j],
                                             r1[:, :, iiPDF_Ncn], r2[:, :, j]))

    # Transform the correlations (in-precip) between two precipitating
    # hydrometeors to standard space.
    for ii in range(iiPDF_Ncn + 1, pdf_dim - 1):
        for jj in range(ii + 1, pdf_dim):
            # Transform the correlation (in-precip) between two precipitating
            # hydrometeors (for example, r_r and N_r) to standard space in PDF
            # component 1.
            A1 = _s(A1, jj, ii, corr_NN2LL(Cn1[:, :, jj, ii], s1n[:, :, ii], s1n[:, :, jj],
                                           r1[:, :, ii], r1[:, :, jj]))
            # Transform the correlation (in-precip) between two precipitating
            # hydrometeors (for example, r_r and N_r) to standard space in PDF
            # component 2.
            A2 = _s(A2, jj, ii, corr_NN2LL(Cn2[:, :, jj, ii], s2n[:, :, ii], s2n[:, :, jj],
                                           r2[:, :, ii], r2[:, :, jj]))

    eye = jnp.zeros((ng, nzt, pdf_dim, pdf_dim)).at[:, :, idx, idx].set(1.0)
    A1 = A1 + jnp.swapaxes(A1, -1, -2) - eye
    A2 = A2 + jnp.swapaxes(A2, -1, -2) - eye
    return A1, A2


def calc_corr_norm_and_cholesky_factor(corr_array_n_cloud, corr_array_n_below, rc_1, rc_2,
                                       iiPDF_type, iiPDF_chi, iiPDF_eta, iiPDF_w, iiPDF_Ncn,
                                       l_follow_ADG1_PDF_standards):
    """This subroutine computes the correlation arrays and correlation
       Cholesky matrices of PDF vars for both components. Here, we assume that
       there are only two unique correlation arrays, which allows us to compute
       these two unique arrays and their corresponding Cholesky decompositions,
       then use rc to determine which one to assign to each grid column and
       vertical level. If the correlation arrays vary based on vertically varying
       values, then this subroutine is not appropriate.

    References:
    https://github.com/larson-group/cam/issues/129#issuecomment-816205563

    Returns:
        (corr_1, corr_2, chol_1, chol_2)

    Fortran declarations (calc_corr_norm_and_cholesky_factor in src/CLUBB_core/setup_clubb_pdf_params.F90):

        !--------------------------- Input Variables ---------------------------
        integer, intent(in) :: &
          nzt,     & ! Number of thermodynamic vertical levels
          pdf_dim, & ! Number of variables in the corr/mean/stdev arrays
          ngrdcol    ! Number of grid columns

        real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
          rc_1,   & ! Mean of r_c (1st PDF component)                 [kg/kg]
          rc_2      ! Mean of r_c (2nd PDF component)                 [kg/kg]

        real( kind = core_rknd ), dimension(pdf_dim,pdf_dim), intent(in) :: &
          corr_array_n_cloud, & ! Prescribed correlation array in cloud        [-]
          corr_array_n_below    ! Prescribed correlation array below cloud     [-]

        type (hm_metadata_type), intent(in) :: &
          hm_metadata

        integer, intent(in) :: &
          iiPDF_type    ! Selected option for the two-component normal (double
                        ! Gaussian) PDF type to use for the w, rt, and theta-l (or
                        ! w, chi, and eta) portion of CLUBB's multivariate,
                        ! two-component PDF.

        !--------------------------- Output Variables ---------------------------
        real( kind = core_rknd ), dimension(ngrdcol,nzt,pdf_dim,pdf_dim), intent(out) :: &
          corr_array_1_n,       & ! Corr. array (normal space) of PDF vars. (comp. 1)  [-]
          corr_array_2_n,       & ! Corr. array (normal space) of PDF vars. (comp. 2)  [-]
          corr_cholesky_mtx_1,  & ! Transposed corr. cholesky matrix, 1st comp.        [-]
          corr_cholesky_mtx_2     ! Transposed corr. cholesky matrix, 2nd comp.        [-]
    """
    cc = corr_array_n_cloud
    cb = corr_array_n_below
    pdf_dim = cc.shape[0]
    hm_indices = list(range(iiPDF_Ncn + 1, pdf_dim))

    # ADG standards fix corr(w, chi) = corr(w, eta) = 0.
    if l_follow_ADG1_PDF_standards and iiPDF_type in (IIPDF_TYPE_ADG1, IIPDF_TYPE_ADG2, IIPDF_TYPE_NEW_HYBRID):
        for (r, c) in ((iiPDF_w, iiPDF_chi), (iiPDF_w, iiPDF_eta)):
            cc = cc.at[r, c].set(0.0); cb = cb.at[r, c].set(0.0)

    # Ncn is inherently in-cloud: replace its below-cloud correlations with the in-cloud ones.
    cb = cb.at[iiPDF_Ncn, iiPDF_chi].set(cc[iiPDF_Ncn, iiPDF_chi])
    cb = cb.at[iiPDF_Ncn, iiPDF_eta].set(cc[iiPDF_Ncn, iiPDF_eta])

    # eta–hm correlation estimated as chi–eta * chi–hm (keeps the matrix Cholesky-decomposable).
    for j in hm_indices:
        cc = cc.at[j, iiPDF_eta].set(cc[iiPDF_eta, iiPDF_chi] * cc[j, iiPDF_chi])
        cb = cb.at[j, iiPDF_eta].set(cb[iiPDF_eta, iiPDF_chi] * cb[j, iiPDF_chi])

    # Symmetrize from the (modified) lower triangle, then factor.
    def _symm(M):
        L = jnp.tril(M)
        return L + jnp.swapaxes(L, -1, -2) - jnp.diag(jnp.diagonal(M))

    cc_s = _symm(cc); cb_s = _symm(cb)
    _, chol_cloud, _ = Cholesky_factor(cc_s)
    _, chol_below, _ = Cholesky_factor(cb_s)
    chol_cloud = jnp.tril(chol_cloud)        # zero the upper triangle
    chol_below = jnp.tril(chol_below)

    # Assign per column/level by rc (the "two unique arrays" selection).
    def _select(rc):
        sel = (rc > rc_tol)[:, :, None, None]
        corr = jnp.where(sel, cc_s, cb_s)
        chol = jnp.where(sel, chol_cloud, chol_below)
        return corr, chol

    corr_1, chol_1 = _select(rc_1)
    corr_2, chol_2 = _select(rc_2)
    return corr_1, corr_2, chol_1, chol_2


def pdf_param_hm_stats(nzt, ngrdcol, pdf_dim, hydromet_dim, hm_metadata,
                       hm_1, hm_2,
                       mu_x_1, mu_x_2,
                       sigma_x_1, sigma_x_2,
                       corr_array_1, corr_array_2,
                       stats):
    """Record statistics for standard PDF parameters involving hydrometeors.

    Returns:
        stats

    Fortran declarations (pdf_param_hm_stats in src/CLUBB_core/setup_clubb_pdf_params.F90):

        !--------------------------- Input Variables ---------------------------
        integer, intent(in) :: &
          nzt,          & ! Number of thermodynamic vertical levels
          ngrdcol,      & ! Number of grid columns
          pdf_dim,      & ! Number of variables in the correlation array
          hydromet_dim    ! Number of hydrometeor species

        type (hm_metadata_type), intent(in) :: &
          hm_metadata

        real( kind = core_rknd ), dimension(ngrdcol,nzt,hydromet_dim), intent(in) :: &
          hm_1, & ! Mean of a precip. hydrometeor (1st PDF component)   [units vary]
          hm_2    ! Mean of a precip. hydrometeor (2nd PDF component)   [units vary]

        real( kind = core_rknd ), dimension(ngrdcol,nzt,pdf_dim), intent(in) :: &
          mu_x_1,    & ! Mean array of PDF vars. (1st PDF component)    [units vary]
          mu_x_2,    & ! Mean array of PDF vars. (2nd PDF component)    [units vary]
          sigma_x_1, & ! Standard deviation array of PDF vars (comp. 1) [units vary]
          sigma_x_2    ! Standard deviation array of PDF vars (comp. 2) [units vary]

        real( kind = core_rknd ), dimension(ngrdcol,nzt,pdf_dim,pdf_dim), intent(in) :: &
          corr_array_1, & ! Correlation array of PDF vars. (comp. 1)             [-]
          corr_array_2    ! Correlation array of PDF vars. (comp. 2)             [-]

        type(stats_type), intent(inout) :: &
          stats
    """
    del nzt, ngrdcol

    if not stats.l_sample:
        return stats

    iiPDF_chi = hm_metadata.iiPDF_chi
    iiPDF_eta = hm_metadata.iiPDF_eta
    iiPDF_w = hm_metadata.iiPDF_w
    iiPDF_Ncn = hm_metadata.iiPDF_Ncn

    for ivar in range(hydromet_dim):
        hm_type = hm_metadata.hydromet_list[ivar]

        # Mean of the precipitating hydrometeor in PDF component 1.
        var_name = f"{hm_type[:2].strip()}_1"
        stats = stats.update(var_name, hm_1[:, :, ivar])

        # Mean of the precipitating hydrometeor in PDF component 2.
        var_name = f"{hm_type[:2].strip()}_2"
        stats = stats.update(var_name, hm_2[:, :, ivar])

    for ivar in range(iiPDF_Ncn + 1, pdf_dim):
        hm_idx = pdf2hydromet_idx(ivar, hm_metadata)
        hm_type = hm_metadata.hydromet_list[hm_idx]

        # Mean of the precipitating hydrometeor (in-precip) in PDF component 1.
        var_name = f"mu_{hm_type[:2].strip()}_1"
        stats = stats.update(var_name, mu_x_1[:, :, ivar])

        # Mean of the precipitating hydrometeor (in-precip) in PDF component 2.
        var_name = f"mu_{hm_type[:2].strip()}_2"
        stats = stats.update(var_name, mu_x_2[:, :, ivar])

    # Mean of cloud nuclei concentration in PDF component 1.
    stats = stats.update("mu_Ncn_1", mu_x_1[:, :, iiPDF_Ncn])
    # Mean of cloud nuclei concentration in PDF component 2.
    stats = stats.update("mu_Ncn_2", mu_x_2[:, :, iiPDF_Ncn])

    for ivar in range(iiPDF_Ncn + 1, pdf_dim):
        hm_idx = pdf2hydromet_idx(ivar, hm_metadata)
        hm_type = hm_metadata.hydromet_list[hm_idx]

        # Standard deviation of the precipitating hydrometeor (in-precip) in PDF component 1.
        var_name = f"sigma_{hm_type[:2].strip()}_1"
        stats = stats.update(var_name, sigma_x_1[:, :, ivar])

        # Standard deviation of the precipitating hydrometeor (in-precip) in PDF component 2.
        var_name = f"sigma_{hm_type[:2].strip()}_2"
        stats = stats.update(var_name, sigma_x_2[:, :, ivar])

    # Standard deviation of cloud nuclei concentration in PDF component 1.
    stats = stats.update("sigma_Ncn_1", sigma_x_1[:, :, iiPDF_Ncn])
    # Standard deviation of cloud nuclei concentration in PDF component 2.
    stats = stats.update("sigma_Ncn_2", sigma_x_2[:, :, iiPDF_Ncn])

    # Correlations of w and chi/eta found in the correlation arrays.
    stats = stats.update("corr_w_chi_1_ca", corr_array_1[:, :, iiPDF_w, iiPDF_chi])
    stats = stats.update("corr_w_chi_2_ca", corr_array_2[:, :, iiPDF_w, iiPDF_chi])
    stats = stats.update("corr_w_eta_1_ca", corr_array_1[:, :, iiPDF_w, iiPDF_eta])
    stats = stats.update("corr_w_eta_2_ca", corr_array_2[:, :, iiPDF_w, iiPDF_eta])

    for ivar in range(iiPDF_Ncn + 1, pdf_dim):
        hm_idx = pdf2hydromet_idx(ivar, hm_metadata)
        hm_type = hm_metadata.hydromet_list[hm_idx]

        # Correlation (in-precip) of w and the precipitating hydrometeor in PDF component 1.
        var_name = f"corr_w_{hm_type[:2].strip()}_1"
        stats = stats.update(var_name, corr_array_1[:, :, ivar, iiPDF_w])

        # Correlation (in-precip) of w and the precipitating hydrometeor in PDF component 2.
        var_name = f"corr_w_{hm_type[:2].strip()}_2"
        stats = stats.update(var_name, corr_array_2[:, :, ivar, iiPDF_w])

    # Correlation of w and N_cn in PDF component 1.
    stats = stats.update("corr_w_Ncn_1", corr_array_1[:, :, iiPDF_Ncn, iiPDF_w])
    # Correlation of w and N_cn in PDF component 2.
    stats = stats.update("corr_w_Ncn_2", corr_array_2[:, :, iiPDF_Ncn, iiPDF_w])

    # Correlation of chi and eta in each PDF component found in the correlation arrays.
    stats = stats.update("corr_chi_eta_1_ca", corr_array_1[:, :, iiPDF_eta, iiPDF_chi])
    stats = stats.update("corr_chi_eta_2_ca", corr_array_2[:, :, iiPDF_eta, iiPDF_chi])

    for ivar in range(iiPDF_Ncn + 1, pdf_dim):
        hm_idx = pdf2hydromet_idx(ivar, hm_metadata)
        hm_type = hm_metadata.hydromet_list[hm_idx]

        # Correlation (in-precip) of chi (old s) and the precipitating hydrometeor in PDF component 1.
        var_name = f"corr_chi_{hm_type[:2].strip()}_1"
        stats = stats.update(var_name, corr_array_1[:, :, ivar, iiPDF_chi])

        # Correlation (in-precip) of chi (old s) and the precipitating hydrometeor in PDF component 2.
        var_name = f"corr_chi_{hm_type[:2].strip()}_2"
        stats = stats.update(var_name, corr_array_2[:, :, ivar, iiPDF_chi])

    # Correlation of chi (old s) and N_cn in PDF component 1.
    stats = stats.update("corr_chi_Ncn_1", corr_array_1[:, :, iiPDF_Ncn, iiPDF_chi])
    # Correlation of chi (old s) and N_cn in PDF component 2.
    stats = stats.update("corr_chi_Ncn_2", corr_array_2[:, :, iiPDF_Ncn, iiPDF_chi])

    for ivar in range(iiPDF_Ncn + 1, pdf_dim):
        hm_idx = pdf2hydromet_idx(ivar, hm_metadata)
        hm_type = hm_metadata.hydromet_list[hm_idx]

        # Correlation (in-precip) of eta (old t) and the precipitating hydrometeor in PDF component 1.
        var_name = f"corr_eta_{hm_type[:2].strip()}_1"
        stats = stats.update(var_name, corr_array_1[:, :, ivar, iiPDF_eta])

        # Correlation (in-precip) of eta (old t) and the precipitating hydrometeor in PDF component 2.
        var_name = f"corr_eta_{hm_type[:2].strip()}_2"
        stats = stats.update(var_name, corr_array_2[:, :, ivar, iiPDF_eta])

    # Correlation of eta (old t) and N_cn in PDF component 1.
    stats = stats.update("corr_eta_Ncn_1", corr_array_1[:, :, iiPDF_Ncn, iiPDF_eta])
    # Correlation of eta (old t) and N_cn in PDF component 2.
    stats = stats.update("corr_eta_Ncn_2", corr_array_2[:, :, iiPDF_Ncn, iiPDF_eta])

    for ivar in range(iiPDF_Ncn + 1, pdf_dim):
        hm_idx = pdf2hydromet_idx(ivar, hm_metadata)
        hm_type = hm_metadata.hydromet_list[hm_idx]

        # Correlation (in-precip) of N_cn and the precipitating hydrometeor in PDF component 1.
        var_name = f"corr_Ncn_{hm_type[:2].strip()}_1"
        stats = stats.update(var_name, corr_array_1[:, :, ivar, iiPDF_Ncn])

        # Correlation (in-precip) of N_cn and the precipitating hydrometeor in PDF component 2.
        var_name = f"corr_Ncn_{hm_type[:2].strip()}_2"
        stats = stats.update(var_name, corr_array_2[:, :, ivar, iiPDF_Ncn])

    for ivar in range(iiPDF_Ncn + 1, pdf_dim):
        hm_idx_ivar = pdf2hydromet_idx(ivar, hm_metadata)
        hmx_type = hm_metadata.hydromet_list[hm_idx_ivar]

        for jvar in range(ivar + 1, pdf_dim):
            hm_idx_jvar = pdf2hydromet_idx(jvar, hm_metadata)
            hmy_type = hm_metadata.hydromet_list[hm_idx_jvar]

            # Correlation (in-precip) of two different hydrometeors in PDF component 1.
            var_name = f"corr_{hmx_type[:2].strip()}_{hmy_type[:2].strip()}_1"
            stats = stats.update(var_name, corr_array_1[:, :, jvar, ivar])

            # Correlation (in-precip) of two different hydrometeors in PDF component 2.
            var_name = f"corr_{hmx_type[:2].strip()}_{hmy_type[:2].strip()}_2"
            stats = stats.update(var_name, corr_array_2[:, :, jvar, ivar])

    return stats


def pdf_param_ln_hm_stats(nzt, ngrdcol, pdf_dim, hm_metadata,
                          mu_x_1_n, mu_x_2_n,
                          sigma_x_1_n, sigma_x_2_n,
                          corr_array_1_n,
                          corr_array_2_n,
                          stats):
    """Record statistics for normal space PDF parameters involving hydrometeors.

    Returns:
        stats

    Fortran declarations (pdf_param_ln_hm_stats in src/CLUBB_core/setup_clubb_pdf_params.F90):

        !------------------------ Input Variables ------------------------
        integer, intent(in) :: &
          nzt,      & ! Number of thermodynamic vertical levels
          ngrdcol,  & ! Number of grid columns
          pdf_dim     ! Number of variables in the correlation array

        type (hm_metadata_type), intent(in) :: &
          hm_metadata

        real( kind = core_rknd ), dimension(ngrdcol,nzt,pdf_dim), intent(in) :: &
          mu_x_1_n,    & ! Mean array (normal space): PDF vars. (comp. 1) [un. vary]
          mu_x_2_n,    & ! Mean array (normal space): PDF vars. (comp. 2) [un. vary]
          sigma_x_1_n, & ! Std. dev. array (normal space): PDF vars (comp. 1) [u.v.]
          sigma_x_2_n    ! Std. dev. array (normal space): PDF vars (comp. 2) [u.v.]

        real( kind = core_rknd ), dimension(ngrdcol,nzt,pdf_dim,pdf_dim), &
        intent(in) :: &
          corr_array_1_n, & ! Corr. array (normal space) of PDF vars. (comp. 1)  [-]
          corr_array_2_n    ! Corr. array (normal space) of PDF vars. (comp. 2)  [-]

        type(stats_type), intent(inout) :: &
          stats
    """
    del nzt, ngrdcol

    if not stats.l_sample:
        return stats

    iiPDF_chi = hm_metadata.iiPDF_chi
    iiPDF_eta = hm_metadata.iiPDF_eta
    iiPDF_w = hm_metadata.iiPDF_w
    iiPDF_Ncn = hm_metadata.iiPDF_Ncn

    for ivar in range(iiPDF_Ncn + 1, pdf_dim):
        hm_idx = pdf2hydromet_idx(ivar, hm_metadata)
        hm_type = hm_metadata.hydromet_list[hm_idx]

        var_name = f"mu_{hm_type[:2].strip()}_1_n"
        if stats.var_on_stats_list(var_name):
            stats = stats.update(var_name, mu_x_1_n[:, :, ivar])

        var_name = f"mu_{hm_type[:2].strip()}_2_n"
        if stats.var_on_stats_list(var_name):
            stats = stats.update(var_name, mu_x_2_n[:, :, ivar])

    if stats.var_on_stats_list("mu_Ncn_1_n"):
        stats = stats.update("mu_Ncn_1_n", mu_x_1_n[:, :, iiPDF_Ncn])

    if stats.var_on_stats_list("mu_Ncn_2_n"):
        stats = stats.update("mu_Ncn_2_n", mu_x_2_n[:, :, iiPDF_Ncn])

    for ivar in range(iiPDF_Ncn + 1, pdf_dim):
        hm_idx = pdf2hydromet_idx(ivar, hm_metadata)
        hm_type = hm_metadata.hydromet_list[hm_idx]

        var_name = f"sigma_{hm_type[:2].strip()}_1_n"
        stats = stats.update(var_name, sigma_x_1_n[:, :, ivar])

        var_name = f"sigma_{hm_type[:2].strip()}_2_n"
        stats = stats.update(var_name, sigma_x_2_n[:, :, ivar])

    stats = stats.update("sigma_Ncn_1_n", sigma_x_1_n[:, :, iiPDF_Ncn])
    stats = stats.update("sigma_Ncn_2_n", sigma_x_2_n[:, :, iiPDF_Ncn])

    for ivar in range(iiPDF_Ncn + 1, pdf_dim):
        hm_idx = pdf2hydromet_idx(ivar, hm_metadata)
        hm_type = hm_metadata.hydromet_list[hm_idx]

        # Correlation (in-precip) of w and ln hm in PDF component 1.
        var_name = f"corr_w_{hm_type[:2].strip()}_1_n"
        stats = stats.update(var_name, corr_array_1_n[:, :, ivar, iiPDF_w])

        # Correlation (in-precip) of w and ln hm in PDF component 2.
        var_name = f"corr_w_{hm_type[:2].strip()}_2_n"
        stats = stats.update(var_name, corr_array_2_n[:, :, ivar, iiPDF_w])

    # Correlation of w and ln N_cn in PDF component 1.
    stats = stats.update("corr_w_Ncn_1_n", corr_array_1_n[:, :, iiPDF_Ncn, iiPDF_w])
    # Correlation of w and ln N_cn in PDF component 2.
    stats = stats.update("corr_w_Ncn_2_n", corr_array_2_n[:, :, iiPDF_Ncn, iiPDF_w])

    for ivar in range(iiPDF_Ncn + 1, pdf_dim):
        hm_idx = pdf2hydromet_idx(ivar, hm_metadata)
        hm_type = hm_metadata.hydromet_list[hm_idx]

        # Correlation (in-precip) of chi (old s) and ln hm in PDF component 1.
        var_name = f"corr_chi_{hm_type[:2].strip()}_1_n"
        stats = stats.update(var_name, corr_array_1_n[:, :, ivar, iiPDF_chi])

        # Correlation (in-precip) of chi (old s) and ln hm in PDF component 2.
        var_name = f"corr_chi_{hm_type[:2].strip()}_2_n"
        stats = stats.update(var_name, corr_array_2_n[:, :, ivar, iiPDF_chi])

    # Correlation of chi (old s) and ln N_cn in PDF component 1.
    stats = stats.update("corr_chi_Ncn_1_n", corr_array_1_n[:, :, iiPDF_Ncn, iiPDF_chi])
    # Correlation of chi (old s) and ln N_cn in PDF component 2.
    stats = stats.update("corr_chi_Ncn_2_n", corr_array_2_n[:, :, iiPDF_Ncn, iiPDF_chi])

    for ivar in range(iiPDF_Ncn + 1, pdf_dim):
        hm_idx = pdf2hydromet_idx(ivar, hm_metadata)
        hm_type = hm_metadata.hydromet_list[hm_idx]

        # Correlation (in-precip) of eta (old t) and ln hm in PDF component 1.
        var_name = f"corr_eta_{hm_type[:2].strip()}_1_n"
        stats = stats.update(var_name, corr_array_1_n[:, :, ivar, iiPDF_eta])

        # Correlation (in-precip) of eta (old t) and ln hm in PDF component 2.
        var_name = f"corr_eta_{hm_type[:2].strip()}_2_n"
        stats = stats.update(var_name, corr_array_2_n[:, :, ivar, iiPDF_eta])

    # Correlation of eta (old t) and ln N_cn in PDF component 1.
    stats = stats.update("corr_eta_Ncn_1_n", corr_array_1_n[:, :, iiPDF_Ncn, iiPDF_eta])
    # Correlation of eta (old t) and ln N_cn in PDF component 2.
    stats = stats.update("corr_eta_Ncn_2_n", corr_array_2_n[:, :, iiPDF_Ncn, iiPDF_eta])

    for ivar in range(iiPDF_Ncn + 1, pdf_dim):
        hm_idx = pdf2hydromet_idx(ivar, hm_metadata)
        hm_type = hm_metadata.hydromet_list[hm_idx]

        # Correlation (in-precip) of ln N_cn and ln hm in PDF component 1.
        var_name = f"corr_Ncn_{hm_type[:2].strip()}_1_n"
        stats = stats.update(var_name, corr_array_1_n[:, :, ivar, iiPDF_Ncn])

        # Correlation (in-precip) of ln N_cn and ln hm in PDF component 2.
        var_name = f"corr_Ncn_{hm_type[:2].strip()}_2_n"
        stats = stats.update(var_name, corr_array_2_n[:, :, ivar, iiPDF_Ncn])

    for ivar in range(iiPDF_Ncn + 1, pdf_dim):
        hm_idx_ivar = pdf2hydromet_idx(ivar, hm_metadata)
        hmx_type = hm_metadata.hydromet_list[hm_idx_ivar]

        for jvar in range(ivar + 1, pdf_dim):
            hm_idx_jvar = pdf2hydromet_idx(jvar, hm_metadata)
            hmy_type = hm_metadata.hydromet_list[hm_idx_jvar]

            # Correlation (in-precip) of ln hmx and ln hmy in PDF component 1.
            var_name = f"corr_{hmx_type[:2].strip()}_{hmy_type[:2].strip()}_1_n"
            stats = stats.update(var_name, corr_array_1_n[:, :, jvar, ivar])

            # Correlation (in-precip) of ln hmx and ln hmy in PDF component 2.
            var_name = f"corr_{hmx_type[:2].strip()}_{hmy_type[:2].strip()}_2_n"
            stats = stats.update(var_name, corr_array_2_n[:, :, jvar, ivar])

    return stats


def pack_hydromet_pdf_params(nzt, ngrdcol,
                             hydromet_dim, hm_metadata,
                             hm_1, hm_2, pdf_dim, mu_x_1,
                             mu_x_2, sigma_x_1, sigma_x_2,
                             corr_array_1, corr_array_2,
                             hydromet_pdf_params):
    """Pack the standard means and variances involving hydrometeors, as well as a
    few other variables, into the structure hydromet_pdf_params.

    Fortran packs an (ngrdcol, nzt) array of derived types; the JAX object
    stores those per-cell fields as leading array dimensions on each field of
    one HydrometPdfParameter pytree.

    Fortran declarations (pack_hydromet_pdf_params in src/CLUBB_core/setup_clubb_pdf_params.F90):

        !-------------------------- Input Variables --------------------------
        integer, intent(in) :: &
          nzt,          & ! Number of thermodynamic vertical grid levels
          ngrdcol,      & ! Number of grid columns
          hydromet_dim    ! Number of hydrometeor species

        real( kind = core_rknd ), dimension(ngrdcol,nzt,hydromet_dim), intent(in) :: &
          hm_1, & ! Mean of a precip. hydrometeor (1st PDF component)  [units vary]
          hm_2    ! Mean of a precip. hydrometeor (2nd PDF component)  [units vary]

        integer, intent(in) :: &
          pdf_dim   ! Number of variables in the mean/stdev arrays

        type (hm_metadata_type), intent(in) :: &
          hm_metadata

        real( kind = core_rknd ), dimension(ngrdcol,nzt,pdf_dim), intent(in) :: &
          mu_x_1,    & ! Mean array of PDF vars. (1st PDF component)    [units vary]
          mu_x_2,    & ! Mean array of PDF vars. (2nd PDF component)    [units vary]
          sigma_x_1, & ! Standard deviation array of PDF vars (comp. 1) [units vary]
          sigma_x_2    ! Standard deviation array of PDF vars (comp. 2) [units vary]

        real( kind = core_rknd ), dimension(ngrdcol,nzt,pdf_dim,pdf_dim), intent(in) :: &
          corr_array_1, & ! Correlation array of PDF vars. (comp. 1)    [-]
          corr_array_2    ! Correlation array of PDF vars. (comp. 2)    [-]

        !-------------------------- Output Variable --------------------------
        type(hydromet_pdf_parameter), dimension(ngrdcol,nzt), intent(out) :: &
          hydromet_pdf_params    ! Hydrometeor PDF parameters        [units vary]
    """
    del nzt, ngrdcol, pdf_dim, hydromet_pdf_params

    iiPDF_chi = hm_metadata.iiPDF_chi
    iiPDF_eta = hm_metadata.iiPDF_eta
    iiPDF_w = hm_metadata.iiPDF_w
    iiPDF_Ncn = hm_metadata.iiPDF_Ncn

    base_shape = hm_1.shape[:2]
    hm_shape = base_shape + (MAX_HYDROMET_DIM,)
    corr_shape = base_shape + (MAX_HYDROMET_DIM, MAX_HYDROMET_DIM)

    hm_1_out = jnp.zeros(hm_shape)
    hm_2_out = jnp.zeros(hm_shape)
    mu_hm_1 = jnp.zeros(hm_shape)
    mu_hm_2 = jnp.zeros(hm_shape)
    sigma_hm_1 = jnp.zeros(hm_shape)
    sigma_hm_2 = jnp.zeros(hm_shape)
    corr_w_hm_1 = jnp.zeros(hm_shape)
    corr_w_hm_2 = jnp.zeros(hm_shape)
    corr_chi_hm_1 = jnp.zeros(hm_shape)
    corr_chi_hm_2 = jnp.zeros(hm_shape)
    corr_eta_hm_1 = jnp.zeros(hm_shape)
    corr_eta_hm_2 = jnp.zeros(hm_shape)
    corr_hmx_hmy_1 = jnp.zeros(corr_shape)
    corr_hmx_hmy_2 = jnp.zeros(corr_shape)

    # Pack remaining means and standard deviations into hydromet_pdf_params.
    for ivar in range(hydromet_dim):
        pdf_idx = hydromet2pdf_idx(ivar, hm_metadata)

        # Mean of a hydrometeor (overall) in the 1st PDF component.
        hm_1_out = hm_1_out.at[:, :, ivar].set(hm_1[:, :, ivar])
        # Mean of a hydrometeor (overall) in the 2nd PDF component.
        hm_2_out = hm_2_out.at[:, :, ivar].set(hm_2[:, :, ivar])

        # Mean of a hydrometeor (in-precip) in the 1st PDF component.
        mu_hm_1 = mu_hm_1.at[:, :, ivar].set(mu_x_1[:, :, pdf_idx])
        # Mean of a hydrometeor (in-precip) in the 2nd PDF component.
        mu_hm_2 = mu_hm_2.at[:, :, ivar].set(mu_x_2[:, :, pdf_idx])

        # Standard deviation of a hydrometeor (in-precip) in the 1st PDF component.
        sigma_hm_1 = sigma_hm_1.at[:, :, ivar].set(sigma_x_1[:, :, pdf_idx])
        # Standard deviation of a hydrometeor (in-precip) in the 2nd PDF component.
        sigma_hm_2 = sigma_hm_2.at[:, :, ivar].set(sigma_x_2[:, :, pdf_idx])

        # Correlation (in-precip) of w and a hydrometeor in the 1st PDF component.
        corr_w_hm_1 = corr_w_hm_1.at[:, :, ivar].set(corr_array_1[:, :, pdf_idx, iiPDF_w])

        # Correlation (in-precip) of w and a hydrometeor in the 2nd PDF component.
        corr_w_hm_2 = corr_w_hm_2.at[:, :, ivar].set(corr_array_2[:, :, pdf_idx, iiPDF_w])

        # Correlation (in-precip) of chi and a hydrometeor in the 1st PDF component.
        corr_chi_hm_1 = corr_chi_hm_1.at[:, :, ivar].set(corr_array_1[:, :, pdf_idx, iiPDF_chi])

        # Correlation (in-precip) of chi and a hydrometeor in the 2nd PDF component.
        corr_chi_hm_2 = corr_chi_hm_2.at[:, :, ivar].set(corr_array_2[:, :, pdf_idx, iiPDF_chi])

        # Correlation (in-precip) of eta and a hydrometeor in the 1st PDF component.
        corr_eta_hm_1 = corr_eta_hm_1.at[:, :, ivar].set(corr_array_1[:, :, pdf_idx, iiPDF_eta])

        # Correlation (in-precip) of eta and a hydrometeor in the 2nd PDF component.
        corr_eta_hm_2 = corr_eta_hm_2.at[:, :, ivar].set(corr_array_2[:, :, pdf_idx, iiPDF_eta])

        # Correlation (in-precip) of two hydrometeors, hmx and hmy, in the 1st PDF component.
        corr_hmx_hmy_1 = corr_hmx_hmy_1.at[:, :, ivar, ivar].set(1.0)

        for jvar in range(ivar + 1, hydromet_dim):
            jpdf_idx = hydromet2pdf_idx(jvar, hm_metadata)
            corr_hmx_hmy_1 = corr_hmx_hmy_1.at[:, :, jvar, ivar].set(
                corr_array_1[:, :, jpdf_idx, pdf_idx]
            )
            corr_hmx_hmy_1 = corr_hmx_hmy_1.at[:, :, ivar, jvar].set(
                corr_hmx_hmy_1[:, :, jvar, ivar]
            )

        # Correlation (in-precip) of two hydrometeors, hmx and hmy, in the 2nd PDF component.
        corr_hmx_hmy_2 = corr_hmx_hmy_2.at[:, :, ivar, ivar].set(1.0)

        for jvar in range(ivar + 1, hydromet_dim):
            jpdf_idx = hydromet2pdf_idx(jvar, hm_metadata)
            corr_hmx_hmy_2 = corr_hmx_hmy_2.at[:, :, jvar, ivar].set(
                corr_array_2[:, :, jpdf_idx, pdf_idx]
            )
            corr_hmx_hmy_2 = corr_hmx_hmy_2.at[:, :, ivar, jvar].set(
                corr_hmx_hmy_2[:, :, jvar, ivar]
            )

    # Mean of Ncn (overall) in the 1st PDF component.
    mu_Ncn_1 = mu_x_1[:, :, iiPDF_Ncn]
    # Mean of Ncn (overall) in the 2nd PDF component.
    mu_Ncn_2 = mu_x_2[:, :, iiPDF_Ncn]

    # Standard deviation of Ncn (overall) in the 1st PDF component.
    sigma_Ncn_1 = sigma_x_1[:, :, iiPDF_Ncn]
    # Standard deviation of Ncn (overall) in the 2nd PDF component.
    sigma_Ncn_2 = sigma_x_2[:, :, iiPDF_Ncn]

    return HydrometPdfParameter(
        hm_1=hm_1_out,
        hm_2=hm_2_out,
        mu_hm_1=mu_hm_1,
        mu_hm_2=mu_hm_2,
        sigma_hm_1=sigma_hm_1,
        sigma_hm_2=sigma_hm_2,
        corr_w_hm_1=corr_w_hm_1,
        corr_w_hm_2=corr_w_hm_2,
        corr_chi_hm_1=corr_chi_hm_1,
        corr_chi_hm_2=corr_chi_hm_2,
        corr_eta_hm_1=corr_eta_hm_1,
        corr_eta_hm_2=corr_eta_hm_2,
        corr_hmx_hmy_1=corr_hmx_hmy_1,
        corr_hmx_hmy_2=corr_hmx_hmy_2,
        mu_Ncn_1=mu_Ncn_1,
        mu_Ncn_2=mu_Ncn_2,
        sigma_Ncn_1=sigma_Ncn_1,
        sigma_Ncn_2=sigma_Ncn_2,
    )


def compute_rtp2_from_chi(sigma_chi_1, sigma_chi_2, sigma_eta_1, sigma_eta_2,
                          rt_1, rt_2, crt_1, crt_2, mixt_frac,
                          corr_chi_eta_1, corr_chi_eta_2):
    """Compute the variance of rt from the distribution of chi and eta. The
    resulting variance will be consistent with CLUBB's extended PDF
    involving chi and eta, including if l_fix_w_chi_eta_correlations = .true..

    Returns:
        rtp2_zt_from_chi

    Fortran declarations (compute_rtp2_from_chi in src/CLUBB_core/setup_clubb_pdf_params.F90):

        !--------------------------- Input Variables ---------------------------
        real( kind = core_rknd ), intent(in) :: &
          sigma_chi_1,         & ! Standard deviation of chi (1st PDF comp.) [kg/kg]
          sigma_chi_2,         & ! Standard deviation of chi (2nd PDF comp.) [kg/kg]
          sigma_eta_1,         & ! Standard deviation of eta (1st PDF comp.) [kg/kg]
          sigma_eta_2,         & ! Standard deviation of eta (2nd PDF comp.) [kg/kg]
          crt_1,               & ! Coef. of r_t in chi/eta eqns. (1st comp.) [-]
          crt_2,               & ! Coef. of r_t in chi/eta eqns. (2nd comp.) [-]
          rt_1,                & ! Mean of rt (1st PDF component)            [kg/kg]
          rt_2,                & ! Mean of rt (2nd PDF component)            [kg/kg]
          mixt_frac              ! Weight of 1st gaussian PDF component      [-]

        real( kind = core_rknd ), intent(in) :: &
          corr_chi_eta_1, & ! Correlation of chi and eta in 1st PDF component [-]
          corr_chi_eta_2    ! Correlation of chi and eta in 2nd PDF component [-]

        !--------------------------- Output Variables ---------------------------
        real( kind = core_rknd ) :: &
          rtp2_zt_from_chi    ! Grid-box variance of rtp2 on thermo. levels  [kg/kg]
    """
    #----- Begin Code -----
    varnce_rt_1 = (corr_chi_eta_1 * sigma_chi_1 * sigma_eta_1
                   + 0.5 * sigma_chi_1 ** 2 + 0.5 * sigma_eta_1 ** 2) / (2.0 * crt_1 ** 2)
    varnce_rt_2 = (corr_chi_eta_2 * sigma_chi_2 * sigma_eta_2
                   + 0.5 * sigma_chi_2 ** 2 + 0.5 * sigma_eta_2 ** 2) / (2.0 * crt_2 ** 2)
    rtm = mixt_frac * rt_1 + (1.0 - mixt_frac) * rt_2
    sigma_rt_1 = sqrt_clipped(varnce_rt_1)
    sigma_rt_2 = sqrt_clipped(varnce_rt_2)
    return compute_variance_binormal(rtm, rt_1, rt_2, sigma_rt_1, sigma_rt_2, mixt_frac)
