"""Physical/math constants for BUGSrad — port of bugsrad_physconst.F90.

These are BUGSrad's OWN constant values (e.g. R_d=287.0, gravity=9.80665) — distinct from
constants_clubb; the radiation port must use these to stay faithful.
"""

GRAVITY = 9.80665          # m s^-2
CP_DRY_AIR = 1.004e+03     # J kg^-1 K^-1
R_D = 287.0                # J K^-1 kg^-1
R_STAR = 8.3143e+03        # J K^-1 kmol^-1
SOL_CONST = 1.367e+03      # W m^-2
P_STD = 1.01325e+05        # Pa
T_STD = 273.15             # K
PI = 3.1415926535          # (-)
MOLAR_VOLUME = 2.2421e+4   # cm3-atm

# Molecular weights (g/mol)
MW_DRY_AIR = 28.964
MW_H2O = 18.016
MW_CO2 = 44.010
MW_O3 = 48.000
MW_CH4 = 16.042
MW_N2O = 44.016

EPSILON = MW_H2O / MW_DRY_AIR
F_VIRT = (1.0 - EPSILON) / EPSILON
