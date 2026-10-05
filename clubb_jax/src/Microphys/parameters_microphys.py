"""Parameters for microphysical schemes; mirrors parameters_microphys.F90.

Module configuration is initialized on the host before tracing physics kernels.
"""

# Constant parameters: sampled-microphysics feedback modes.
lh_microphys_interactive = 1  # Feed the samples into the microphysics and allow feedback
lh_microphys_non_interactive = 2  # Feed the samples into the microphysics with no feedback
lh_microphys_disabled = 3  # Disable Latin hypercube entirely

# Morrison aerosol options.
morrison_no_aerosol = 0
morrison_power_law = 1
morrison_lognormal = 2

# Module-owned scheme settings, initialized before tracing.
l_cloud_sed = False  # Cloud water sedimentation (K&K/No microphysics)
l_ice_microphys = False  # Compute ice (COAMPS/Morrison)
l_upwind_diff_sed = False  # Use upwind differencing approx for sedimentation (K&K/COAMPS)
l_graupel = False  # Compute graupel (COAMPS/Morrison)
l_hail = False  # Assumption about graupel/hail? (Morrison)
l_seifert_beheng = False  # Use Seifert and Behneng warm drizzle (Morrison)
l_predict_Nc = False  # Predict cloud droplet concentration (Morrison)
l_subgrid_w = True  # Use subgrid w (Morrison)
l_arctic_nucl = False  # Use MPACE observations (Morrison)
l_fix_pgam = False  # Fix pgam (Morrison)
l_in_cloud_Nc_diff = True  # Use in cloud values of Nc for diffusion
l_var_covar_src = False  # Upscaled microphysics sources for predictive variances/covariances

# KK source adjustment normally limits cloud-water depletion at each sample.
# To test convergence to analytic upscaling, disable the per-sample adjustment
# and adjust only the sampled mean instead (source ticket 558).
l_silhs_KK_convergence_adj_mean = False

l_cloud_edge_activation = False  # Activate on cloud edges (Morrison)
l_local_kk = False  # Local drizzle for Khairoutdinov & Kogan microphysics
specify_aerosol = "morrison_lognormal"  # Specify aerosol (Morrison)

# Sampling configuration; the native JAX random key uses the retained seed.
lh_num_samples = 2  # Number of Latin-hypercube samples passed to microphysics
lh_sequence_length = 1  # Number of timesteps before the latin hypercube seq. repeats
lh_seed = 5489  # Seed for native JAX randomness

# Select disabled, interactive, or diagnostic-only sampled microphysics.
lh_microphys_type = 3

# Scheme selection, species sedimentation flags, and activation option.
microphys_scheme = "none"  # khairoutdinov_kogan, simplified_ice, coamps, etc.
l_hydromet_sed = ()  # Flag to sediment each mean hydrometeor field
l_gfdl_activation = False  # GFDL activation (currently gated during initialization)

# Timing and initial cloud-droplet properties.
microphys_start_time = 0.0  # When to start the microphysics      [s]
Nc0_in_cloud = 1.0e8  # Initial cloud droplet concentration    [num/m^3]
sigma_g = 1.5  # Geometric std. dev. of cloud droplets falling in a stokes regime.
