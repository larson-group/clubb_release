"""JAX port of `src/CLUBB_core/model_flags.F90`.

Description:
  Various model options that can be toggled off and on as desired.

References:
  None

Porting deviations:
- The Fortran `clubb_config_flags_type` declaration is represented by
  `clubb_jax.src.CLUBB_core.config_flags.ConfigFlags`; this module constructs
  that type directly instead of using Fortran out-arguments.
- `set_default_clubb_config_flags_api` is exposed as `get_default_config_flags`
  because Python callers expect a returned value.
- `initialize_clubb_config_flags_type_api` is omitted because direct
  `ConfigFlags(...)` construction replaces it.
- `print_clubb_config_flags_api` is omitted because the JAX path does not use
  Fortran unit-based formatted output.
"""
from __future__ import annotations

from clubb_jax.src.CLUBB_core.config_flags import ConfigFlags


# Advance subroutine ordering variables
order_xm_wpxp = 1
order_xp2_xpyp = 2
order_wp2_wp3 = 3
order_windm = 4

# Options for the two component normal (double Gaussian) PDF type to use for
# the w, rt, and theta-l (or w, chi, and eta) portion of CLUBB's multivariate,
# two-component PDF.
iiPDF_ADG1 = 1        # ADG1 PDF
iiPDF_ADG2 = 2        # ADG2 PDF
iiPDF_3D_Luhar = 3    # 3D Luhar PDF
iiPDF_new = 4         # new PDF
iiPDF_TSDADG = 5      # new TSDADG PDF
iiPDF_LY93 = 6        # Lewellen and Yoh (1993)
iiPDF_new_hybrid = 7  # new hybrid PDF

# Options for the placement of the call to CLUBB's PDF.
ipdf_pre_advance_fields = 1       # Call before advancing predictive fields
ipdf_post_advance_fields = 2      # Call after advancing predictive fields
ipdf_pre_post_advance_fields = 3  # Call both before and after advancing
                                  # predictive fields

# Options to set which algorithm the fill_holes routine uses to correct below threshold values
# in field solutions. The fill_holes method attempts to fill in a mass preserving way, in hopes
# of avoiding the need to perform blunt clipping, which can cause surious sources/sinks.
# An important consideration with these method is the locality - moving mass from one grid level
# to a far away one can create unintended non-local effects, so most methods attempt
# to fill with some degree of locality before relying on a global fill.
global_fill = 1          # Fast but minimally local, most methods use this as a last resort.

sliding_window = 2       # Attempt a highly local fill with a sliding window technique,
                         # and falls back to global if local pass failed.

widening_windows = 3     # Attempt to fill within fixed windows of a certain size, then
                         # repeat with increasaingly larger window sizes until all holes
                         # are filled. Window size can increase to entire domain, which
                         # is equivalent to a global fill.

smart_window = 4         # Uses lightweight hueristics to determine ranges to fill in one
                         # pass. This is highly local when possible, maintains some
                         # locality when wide hole ranges are encountered (if possible),
                         # and range can be the whole domain, which is equivalent to
                         # a global fallback. The gauranteed "one pass" feature seems
                         # to cause this to be the fastest method overall, at least with
                         # the current common hole patterns observed in CLUBB.

smart_window_smooth = 5  # Same as smart window, but contains fancy smoothing features.
                         # This could fail if the field is average above (but close to)
                         # threshold. The efficacy of the smoothing features is (currently)
                         # untested, and this is about 25% slower than smart_window without
                         # the smoothing, and has no global fallack (we could add one),
                         # but the smoothing could matter in theory, and looks nice.

parallel_fill = 6        # A parallelizable method that limits the mass each hole can
                         # steal, then considers each hole independently.
                         # Despite the parallizability being an attractive GPU
                         # feature, current timing results suggest this is the slowest
                         # method on a GPU (and CPU).

lapack = 1          # Use lapack library for matrix solves
penta_lu = 2        # Use penta_lu solver for 5 banded matrices
tridiag_lu = 2      # Use tridiag_lu solver for 3 banded matrices
penta_bicgstab = 3  # Use bicgstab to solve 5 banded matrices

# Options for grid_remap_method, define
# the remapping technique to remap the values from one grid to another
# starts at 1, so 0 is an invalid option for this flag
cons_ullrich_remap = 1  # uses the remapping method proposed by Ullrich et al. in
                        # 'Arbitrary-Order Conservative and Consistent Remapping and a
                        #  Theory of Linear Maps: Part II' (Formula (30))
ppm_remap = 2           # uses piecewise parabolic method from E3SM implementation

# Options for grid_adapt_in_time_method, either don't use this setup at all (0) or define
# the variables and the way the grid density function is formed
no_grid_adaptation = 0  # the grid gets initialized once at the start and
                        # stays constant over every timestep (default)
Lscale_and_wp2 = 1      # uses Lscale and wp2 to form a grid density function

l_pos_def = False        # Flux limiting positive definite scheme on rtm
l_clip_turb_adv = False  # Corrects thlm/rtm when w'th_l'/w'r_t' is clipped

l_gamma_Skw = True       # Use a Skw dependent gamma parameter

l_use_boussinesq = False  # Flag to use the Boussinesq form of the
                          # predictive equations.  The predictive
                          # equations are anelastic by default.

# Flag to use explicit turbulent advection in the wp3 predictive equation.
l_explicit_turbulent_adv_wp3 = False
# Flag to use explicit turbulent advection in the wpxp predictive equation.
l_explicit_turbulent_adv_wpxp = False
# Flag to use explicit turbulent advection in the xp2 and xpyp predictive
# equations.
l_explicit_turbulent_adv_xpyp = False

# Flag to advance xp3 using a simplified version of the d(xp3)/dt predictive
# equation or calculate it using a steady-state approximation.  When the flag
# is turned off, the Larson and Golaz (2005) ansatz to calculate xp3 after
# calculating Skx using the ansatz.
l_advance_xp3 = False

l_morr_xp2_mc = False  # Flag to include the effects of rain evaporation
                       # on rtp2 and thlp2.  The moister (rt_1 or rt_2)
                       # and colder (thl_1 or thl_2) will be fed into
                       # the morrison microphys, and rain evaporation will
                       # be allowed to increase variances

l_evaporate_cold_rcm = False  # Flag to evaporate cloud water at temperatures
                              # colder than -37C.  This is to be used for
                              # Morrison microphysics, to prevent excess ice

l_cubic_interp = False  # Flag to convert grid points with cubic monotonic
                        # spline interpolation as opposed to linear interpolation.

l_upwind_Kh_dp_term = False

# These are the integer constants that represent the various saturation
# formulas. To add a new formula, add an additional constant here,
# add the logic to check the strings for the new formula in clubb_core and
# this module, and add logic in saturation to call the proper function--
# the control logic will be based on these named constants.
saturation_bolton = 1  # Constant for Bolton approximations of saturation
saturation_gfdl = 2    # Constant for the GFDL approximation of saturation
saturation_flatau = 3  # Constant for Flatau approximations of saturation
saturation_lookup = 4  # Use a lookup table for mixing length
                       # saturation vapor pressure calculations

# -----------------------------------------------------------------------------
# Options that can be changed at runtime
# The default values are chosen below and overwritten if desired by the user
# -----------------------------------------------------------------------------

# Use a quintic polynomial in mono_cubic_interp
l_quintic_poly_interp = False

l_silhs_rad = False  # Resolve radiation over subcolumns using SILHS

# Previously used within 'ifdef GFDL'
I_sat_sphum = False  # h1g, 2010-06-15

# This flag is only enabled when performing a generalized grid
# (ascending vs. descending grid) test.
l_test_grid_generalization = False

# Forces our matrices to be solved in descending mode, useful for the grid_generalization test
l_force_descending_solves = False


def get_default_config_flags() -> ConfigFlags:
    """Sets all CLUBB flags to a default setting.

    References:
      None
    """
    # Begin code
    # WARNING: THE DEFAULT VALUES OF THE FLAGS BELOW MAY BE OVERWRITTEN
    #    BY NAMELIST VALUES FROM, E.G., configurable_clubb_flags_nl!!!
    return ConfigFlags(
        iiPDF_type=iiPDF_ADG1,
        ipdf_call_placement=ipdf_post_advance_fields,
        penta_solve_method=penta_lu,
        tridiag_solve_method=tridiag_lu,
        saturation_formula=saturation_flatau,
        grid_remap_method=ppm_remap,
        grid_adapt_in_time_method=no_grid_adaptation,
        fill_holes_type=sliding_window,
        l_use_precip_frac=True,
        l_predict_upwp_vpwp=True,
        l_ho_nontrad_coriolis=False,
        l_ho_trad_coriolis=False,
        l_min_wp2_from_corr_wx=False,
        l_min_xp2_from_corr_wx=True,
        l_C2_cloud_frac=False,
        l_diffuse_rtm_and_thlm=False,
        l_stability_correct_Kh_N2_zm=False,
        l_calc_thlp2_rad=True,
        l_upwind_xpyp_ta=True,
        l_upwind_xm_ma=True,
        l_uv_nudge=False,
        l_rtm_nudge=False,
        l_vert_avg_closure=False,
        l_trapezoidal_rule_zt=False,
        l_trapezoidal_rule_zm=False,
        l_call_pdf_closure_twice=False,
        l_standard_term_ta=False,
        l_partial_upwind_wp3=False,
        l_godunov_upwind_wpxp_ta=False,
        l_godunov_upwind_xpyp_ta=False,
        l_use_cloud_cover=False,
        l_diagnose_correlations=False,
        l_calc_w_corr=False,
        l_const_Nc_in_cloud=False,
        l_fix_w_chi_eta_correlations=True,
        l_stability_correct_tau_zm=False,
        l_damp_wp2_using_em=True,
        l_do_expldiff_rtm_thlm=False,
        l_Lscale_plume_centered=False,
        l_diag_Lscale_from_tau=False,
        l_use_C7_Richardson=True,
        l_use_C11_Richardson=False,
        l_use_shear_Richardson=False,
        l_brunt_vaisala_freq_moist=False,
        l_use_thvm_in_bv_freq=False,
        l_rcm_supersat_adj=True,
        l_damp_wp3_Skw_squared=True,
        l_prescribed_avg_deltaz=False,
        l_lmm_stepping=False,
        l_e3sm_config=False,
        l_vary_convect_depth=False,
        l_use_tke_in_wp3_pr_turb_term=True,
        l_use_tke_in_wp2_wp3_K_dfsn=False,
        l_use_wp3_lim_with_smth_Heaviside=False,
        l_smooth_Heaviside_tau_wpxp=False,
        l_modify_limiters_for_cnvg_test=False,
        l_enable_relaxed_clipping=False,
        l_linearize_pbl_winds=False,
        l_mono_flux_lim_thlm=True,
        l_mono_flux_lim_rtm=True,
        l_mono_flux_lim_um=True,
        l_mono_flux_lim_vm=True,
        l_mono_flux_lim_spikefix=True,
        l_host_applies_sfc_fluxes=False,
        l_wp2_fill_holes_tke=True,
        l_add_dycore_grid=False,
    )
