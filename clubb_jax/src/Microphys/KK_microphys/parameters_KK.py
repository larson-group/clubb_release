"""Pure-Python port of Microphys/KK_microphys/parameters_KK.F90.

Parameters for KK (Khairoutdinov & Kogan 2000) microphysics: the fixed exponents in the
autoconversion / accretion / evaporation / mean-volume-radius equations, plus the two
namelist-changeable parameters (the assumed new-drop radius r_0 and the evaporation
coefficient C_evap). Values mirror parameters_KK.F90 exactly.
"""
from clubb_jax.src.CLUBB_core.constants_clubb import one, one_third, two_thirds

# ── Values of exponents in KK microphysics (parameters_KK.F90:24-34, `parameter, public`) ──
KK_evap_Supersat_exp = one          # Exponent on Supersaturation (S) in KK evap. eq.;  1
KK_evap_rr_exp       = one_third    # Exponent on r_r in KK evaporation eq.;            1/3
KK_evap_Nr_exp       = two_thirds   # Exponent on N_r in KK evaporation eq.;            2/3
KK_auto_rc_exp       = 2.47         # Exponent on r_c in KK autoconversion eq.;         2.47
KK_auto_Nc_exp       = -1.79        # Exponent on N_c in KK autoconversion eq.;        -1.79
KK_accr_rc_exp       = 1.15         # Exponent on r_c in KK accretion eq.;              1.15
KK_accr_rr_exp       = 1.15         # Exponent on r_r in KK accretion eq.;              1.15
KK_mvr_rr_exp        = one_third    # Exponent on r_r in KK mean volume radius;         1/3
KK_mvr_Nr_exp        = -one_third   # Exponent on N_r in KK mean volume radius;        -1/3
KK_Nrm_evap_nu       = one          # Exponent (parameter) in <N_r> evaporation eq.;    1

# ── Namelist-changeable (parameters_KK.F90:42-48; not `parameter` because the microphysics_setting
#    namelist can override them — the standalone SCM uses these defaults). ──
r_0    = 25.0e-6   # Assumed radius of all new drops [m]
C_evap = 0.86      # K&K (2000) ratio of drizzle drop mean geometric radius to mean volume radius
