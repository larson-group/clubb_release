"""Parameters for microphysical schemes; mirrors parameters_microphys.F90.

Module configuration is initialized on the host before tracing physics kernels.
"""

lh_microphys_interactive     = 1  # Feed the samples into the microphysics and allow feedback
lh_microphys_non_interactive = 2  # Feed the samples into the microphysics with no feedback
lh_microphys_disabled        = 3  # Disable Latin hypercube entirely
morrison_no_aerosol = 0
morrison_power_law  = 1
morrison_lognormal  = 2
l_cloud_sed = False  # Cloud water sedimentation (K&K/No microphysics)
l_ice_microphys  = False  # Compute ice (COAMPS/Morrison)
l_upwind_diff_sed = False  # Use upwind differencing approx for sedimentation (K&K/COAMPS)
l_graupel  = False  # Compute graupel (COAMPS/Morrison)
l_hail  = False  # Assumption about graupel/hail? (Morrison)
l_seifert_beheng  = False  # Use Seifert and Behneng warm drizzle (Morrison)
l_predict_Nc = False  # Predict cloud droplet conconcentration (Morrison)
l_subgrid_w = True  # Use subgrid w (Morrison)
l_arctic_nucl = False  # Use MPACE observations (Morrison)
l_fix_pgam = False  # Fix pgam (Morrison)
l_in_cloud_Nc_diff = True  # Use in cloud values of Nc for diffusion
l_var_covar_src = False  # Flag for using upscaled microphysics source terms
l_silhs_KK_convergence_adj_mean = False
l_cloud_edge_activation = False  # Activate on cloud edges (Morrison)
l_local_kk              = False  # Local drizzle for Khairoutdinov & Kogan microphysics
specify_aerosol = "morrison_lognormal"  # Specify aerosol (Morrison)
lh_num_samples = 2  # Number of latin hypercube samples to call the microphysics (or
lh_sequence_length = 1  # Number of timesteps before the latin hypercube seq. repeats
lh_seed = 5489  # Seed for the Mersenne
lh_microphys_type = 3
microphys_scheme = "none"  # khairoutdinv_kogan, simplified_ice, coamps, etc.
l_hydromet_sed = ()
l_gfdl_activation = False
microphys_start_time = 0.  # When to start the microphysics      [s]
Nc0_in_cloud = 1.0e8  # Initial cloud droplet concentration    [num/m^3]
sigma_g  = 1.5  # Geometric std. dev. of cloud droplets falling in a stokes regime.
