"""CLUBB constants and parameter indices for the Python port.

Constants mirror src/CLUBB_core/constants_clubb.F90.
Parameter indices mirror src/CLUBB_core/parameter_indices.F90.
Model flags mirror src/CLUBB_core/model_flags.F90.
"""

import numpy as np

# --------------------------------------------------------------------------- #
# Physical / numerical constants  (constants_clubb.F90)
# --------------------------------------------------------------------------- #
grav = 9.81
Cp = 1004.67
Lv = 2.5e6
Rd = 287.04
Rv = 461.5
ep = Rd / Rv
ep1 = (1.0 - ep) / ep
ep2 = 1.0 / ep
kappa = Rd / Cp
p0 = 1.0e5

w_tol = 2.0e-2
w_tol_sqd = w_tol ** 2
thl_tol = 1.0e-2
rt_tol = 1.0e-8
em_min = 1.5 * w_tol_sqd   # minimum TKE

zero = 0.0
one = 1.0
two = 2.0
three_halves = 1.5
one_half = 0.5

eps = max(1.0e-10, np.finfo(np.float64).eps)
zero_threshold = 0.0
unused_var = -999.0
min_max_smth_mag = 1.0e-9
fstderr = 0   # Fortran stderr unit (unused in Python)

cloud_frac_min = 0.005

# --------------------------------------------------------------------------- #
# Parameter indices  (parameter_indices.F90)
# --------------------------------------------------------------------------- #
# NOTE: Fortran uses 1-based indices into the clubb_params(ngrdcol, nparams)
# array.  The Python clubb_params array uses the SAME 1-based convention
# (column 0 is unused / padding) so these constants are used directly.
nparams = 102

ic_K = 37
ic_K1 = 38
ic_K2 = 40
ic_K6 = 42
ic_K8 = 44
ic_K9 = 46
ic_K10 = 74
ic_K10h = 75

iC1 = 1
iC1b = 2
iC2rt = 3
iC2thl = 5
iC2rtthl = 7
iC4 = 8
iC5 = 9
iC6rt = 10
iC6rtb = 11
iC6rtc = 12
iC6thl = 13
iC6thlb = 14
iC6thlc = 15
iC7 = 16
iC7b = 17
iC8 = 18
iC8b = 19
iC10 = 20
iC11 = 21
iC11b = 22
iC11c = 23
iC12 = 24
iC13 = 25
iC14 = 26
iC_wp2_splat = 32

igamma_coef = 57
igamma_coefb = 58
igamma_coefc = 59
imu = 60
ilmin_coef = 62
iomicron = 63
izeta_vrnce_rat = 64
ilambda0_stability_coef = 66

itaumin = 68
itaumax = 69

iup2_sfc_coef = 78

ixp3_coef_base = 90
ixp3_coef_slope = 91

ia3_coef_min = 98
ia_const = 99
ibv_efold = 100

# --------------------------------------------------------------------------- #
# Model flag constants  (model_flags.F90)
# --------------------------------------------------------------------------- #
ipdf_pre_advance_fields = 1
ipdf_post_advance_fields = 2
ipdf_pre_post_advance_fields = 3

iiPDF_ADG1 = 1

l_gamma_Skw = True       # Use Skw-dependent gamma parameter
l_advance_xp3 = False    # Use predictive xp3 equation

# Default advance ordering
order_xm_wpxp = 1
order_xp2_xpyp = 2
order_wp2_wp3 = 3
order_windm = 4

# Local constants used in advance_clubb_core
tau_const = 1000.0
ufmin = 0.01
l_use_invrs_tau_N2_iso = False
l_smooth_min_max = False
smth_type = 2
below_grnd_val = 0.01
