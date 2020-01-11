"""
-------------------------------------------------------------------------------
    G E N E R A L   I N F O R M A T I O N
-------------------------------------------------------------------------------
        | CLUBB BUDGET PLOTTING SETUP |
        -------------------------------
This file contains general constants and information about the CLUBB variables saved
in the netCDF file needed for plotgen.py.

The "plots" array contains setup information about each individual plots,
while the each variable listed in plots lists the lines in that particular plot,
and the netCDF variables used to compute its data.

To add text to a plot, replace the string listed at index 3 with the text you want to appear,
and put the text coordinates as a tuple at index 4 in plots.
The point (0,0) is the origin of the plot, while the point (1,1) is at the upper right corner of the plot.

To generate new plots, add an entry to the plots variable,
create a new list of lines, and put it into index 5 of the new plots entry.

Note for CLUBB: The variables are divided between multiple netCDF files.
We are currently using 1: zm
So newly added variables will have to go into the correct lists.
The user should be careful not to use the same plot name in sortPlots_zm and _zt,
as these will be used as keys in a dictionary.

TODO:   - Construct plot name from long name in netcdf instead
        - Add style entry to lines, 
        - Choice: Keep plot entries as list and use alias for indexing OR Change from lists to dicts OR use pandas
"""
#-------------------------------------------------------------------------------
#   I M P O R T S
#-------------------------------------------------------------------------------
from numpy import nan

#-------------------------------------------------------------------------------
#   C O N S T A N T S
#-------------------------------------------------------------------------------
# Definition of commonly used conversion factors
DAY = 24                                                    # 1d  = 24h
HOUR = 3600                                                 # 1h  = 3600s
KG = 1000.                                                  # 1kg = 1000g
g_per_second_to_kg_per_day = 1. / (DAY * HOUR * KG)
kg_per_second_to_kg_per_day = 1. / (DAY * HOUR)

filler = nan                                                # Define the fill value which should replace invalid values in the data
header = 'CLUBB budgets'                                    # Plot description used for the header on the html page (NOTE: redundant with name?)
name = 'clubb_budgets'                                      # Plot description used for file names and identifying plot types, MUST CONTAIN MODEL NAMES!
prefix = 'CLUBB, '                                          # Prefix identifier used in plot titles
nc_files = ['clubb_zm']                                     # List of NETCDF files containing the data needed for creating the plots listed belo. Paths are defined in case setup files

#-------------------------------------------------------------------------------
# P L O T S
#-------------------------------------------------------------------------------
# Define lines of each plot:
# Define plot line lists
# Line list entry description:
# index | description
#---------------------
#   0   | line name
#   1   | line visibility setting
#   2   | line variable or expression
#   3   | line scaling factor
#   4   | line xaxis selection

# zm

upwp_detailed = [\
    ['upwp_bt', True, 'upwp_bt', 1., 0],\
    ['upwp_ma', True, 'upwp_ma', 1., 0],\
    ['upwp_ta', True, 'upwp_ta', 1., 0],\
    ['upwp_tp', True, 'upwp_tp', 1. ,0],\
    ['upwp_ac', True, 'upwp_ac', 1., 0],\
    ['upwp_bp', True, 'upwp_bp', 1., 0],\
    ['upwp_pr1', True, 'upwp_pr1', 1., 0],\
    ['upwp_pr2', True, 'upwp_pr2', 1., 0],\
    ['upwp_pr3', True, 'upwp_pr3', 1., 0],\
    ['upwp_pr4', True, 'upwp_pr4', 1., 0],\
    ['upwp_dp1', True, 'upwp_dp1', 1., 0],\
    ['upwp_mfl', True, 'upwp_mfl', 1., 0],\
    ['upwp_cl', True, 'upwp_cl', 1., 0],\
    ['upwp_res', True, 'upwp_bt - (upwp_ma + upwp_ta + upwp_tp + upwp_ac + upwp_bp + upwp_pr1 + upwp_pr2 + upwp_pr3 + upwp_pr4 + upwp_dp1 + upwp_mfl + upwp_cl)', 1., 0],\
    ]

vpwp_detailed = [\
    ['vpwp_bt', True, 'vpwp_bt', 1., 0],\
    ['vpwp_ma', True, 'vpwp_ma', 1., 0],\
    ['vpwp_ta', True, 'vpwp_ta', 1., 0],\
    ['vpwp_tp', True, 'vpwp_tp', 1. ,0],\
    ['vpwp_ac', True, 'vpwp_ac', 1., 0],\
    ['vpwp_bp', True, 'vpwp_bp', 1., 0],\
    ['vpwp_pr1', True, 'vpwp_pr1', 1., 0],\
    ['vpwp_pr2', True, 'vpwp_pr2', 1., 0],\
    ['vpwp_pr3', True, 'vpwp_pr3', 1., 0],\
    ['vpwp_pr4', True, 'vpwp_pr4', 1., 0],\
    ['vpwp_dp1', True, 'vpwp_dp1', 1., 0],\
    ['vpwp_mfl', True, 'vpwp_mfl', 1., 0],\
    ['vpwp_cl', True, 'vpwp_cl', 1., 0],\
    ['vpwp_res', True, 'vpwp_bt - (vpwp_ma + vpwp_ta + vpwp_tp + vpwp_ac + vpwp_bp + vpwp_pr1 + vpwp_pr2 + vpwp_pr3 + vpwp_pr4 + vpwp_dp1 + vpwp_mfl + vpwp_cl)', 1., 0],\
    ]

up2_detailed = [
    ['up2_bt', True, 'up2_bt', 1., 0],\
    ['up2_ma', True, 'up2_ma', 1., 0],\
    ['up2_ta', True, 'up2_ta', 1., 0],\
    ['up2_tp', True, 'up2_tp', 1., 0],\
    ['up2_dp1', True, 'up2_dp1', 1., 0],\
    ['up2_dp2', True, 'up2_dp2', 1., 0],\
    ['up2_pr1', True, 'up2_pr1', 1., 0],\
    ['up2_pr2', True, 'up2_pr2', 1., 0],\
    ['up2_sdmp', True, 'up2_sdmp', 1., 0],\
    ['up2_cl', True, 'up2_cl', 1., 0],\
    ['up2_pd', True, 'up2_pd', 1., 0],\
    ['up2_sf', True, 'up2_sf', 1., 0],\
    ['up2_splat', True, 'up2_splat', 1., 0],\
    ['up2_res', True, 'up2_bt - (up2_ma + up2_ta + up2_tp + up2_dp1 + up2_dp2 + up2_pr1 + up2_pr2 + up2_sdmp + up2_cl+ up2_pd + up2_sf + up2_splat)', 1., 0],\
    ]

vp2_detailed = [
    ['vp2_bt', True, 'vp2_bt', 1., 0],\
    ['vp2_ma', True, 'vp2_ma', 1., 0],\
    ['vp2_ta', True, 'vp2_ta', 1., 0],\
    ['vp2_tp', True, 'vp2_tp', 1., 0],\
    ['vp2_dp1', True, 'vp2_dp1', 1., 0],\
    ['vp2_dp2', True, 'vp2_dp2', 1., 0],\
    ['vp2_pr1', True, 'vp2_pr1', 1., 0],\
    ['vp2_pr2', True, 'vp2_pr2', 1., 0],\
    ['vp2_sdmp', True, 'vp2_sdmp', 1., 0],\
    ['vp2_cl', True, 'vp2_cl', 1., 0],\
    ['vp2_pd', True, 'vp2_pd', 1., 0],\
    ['vp2_sf', True, 'vp2_sf', 1., 0],\
    ['vp2_splat', True, 'vp2_splat', 1., 0],\
    ['vp2_res', True, 'vp2_bt - (vp2_ma + vp2_ta + vp2_tp + vp2_dp1 + vp2_dp2 + vp2_pr1 + vp2_pr2 + vp2_sdmp + vp2_cl + vp2_pd + vp2_sf + vp2_splat)', 1., 0],\
    ]

wp2_detailed = [
    ['wp2_bt', True, 'wp2_bt', 1., 0],\
    ['wp2_ma', True, 'wp2_ma', 1., 0],\
    ['wp2_ta', True, 'wp2_ta', 1., 0],\
    ['wp2_ac', True, 'wp2_ac', 1., 0],\
    ['wp2_bp', True, 'wp2_bp', 1., 0],\
    ['wp2_dp1', True, 'wp2_dp1', 1., 0],\
    ['wp2_dp2', True, 'wp2_dp2', 1., 0],\
    ['wp2_pr1', True, 'wp2_pr1', 1., 0],\
    ['wp2_pr2', True, 'wp2_pr2', 1., 0],\
    ['wp2_pr3', True, 'wp2_pr3', 1., 0],\
    ['wp2_sdmp', True, 'wp2_sdmp', 1., 0],\
    ['wp2_cl', True, 'wp2_cl', 1., 0],\
    ['wp2_pd', True, 'wp2_pd', 1., 0],\
    ['wp2_sf', True, 'wp2_sf', 1., 0],\
    ['wp2_splat', True, 'wp2_splat', 1., 0],\
    ['wp2_res', True, 'wp2_bt - (wp2_ma + wp2_ta + wp2_ac + wp2_bp + wp2_dp1 + wp2_dp2 + wp2_pr1 + wp2_pr2 + wp2_pr3  + wp2_sdmp + wp2_cl+ wp2_pd + wp2_sf + wp2_splat)', 1., 0],\
    ]

# Define plots with reduced number of lines:
# In order to preserve color mapping, lines have to be in the following order:
# advection (ta)
# buoyancy (bp)
# dissipation (sum of dp)
# pressure (sum of pr, splat)
# turbulent production (tp)
# time tendency (bt)
# residual (sum of small terms)
# For missing terms, insert dummy entries: ['<line name>', True, None, 1., 0]


upwp_reduced = [\
    ['upwp_dp1', False, 'upwp_dp1', 1., 0],\
    ['upwp_pr1', False, 'upwp_pr1', 1., 0],\
    ['upwp_pr2', False, 'upwp_pr2', 1., 0],\
    ['upwp_pr3', False, 'upwp_pr3', 1., 0],\
    ['upwp_pr4', False, 'upwp_pr4', 1., 0],\
    ['upwp_ac', False, 'upwp_ac', 1., 0],\
    ['upwp_cl', False, 'upwp_cl', 1., 0],\
    ['upwp_ma', False, 'upwp_ma', 1., 0],\
    ['upwp_mfl', False, 'upwp_mfl', 1., 0],\
    ['advection', True, 'upwp_ta', 1., 0],\
    ['buoyancy', True, 'upwp_bp', 1., 0],\
    ['dissipation', True, 'upwp_dp1', 1., 0],\
    ['pressure', True, 'upwp_pr1 + upwp_pr2 + upwp_pr3 + upwp_pr4', 1., 0],\
    ['turb. prod.', True, 'upwp_tp', 1. ,0],\
    ['time tndcy', True, 'upwp_bt', 1., 0],\
    ['residual', True, 'upwp_ac + upwp_cl + upwp_ma + upwp_mfl', 1., 0],\
    ]

vpwp_reduced = [\
    ['vpwp_dp1', False, 'vpwp_dp1', 1., 0],\
    ['vpwp_pr1', False, 'vpwp_pr1', 1., 0],\
    ['vpwp_pr2', False, 'vpwp_pr2', 1., 0],\
    ['vpwp_pr3', False, 'vpwp_pr3', 1., 0],\
    ['vpwp_pr4', False, 'vpwp_pr4', 1., 0],\
    ['vpwp_ac', False, 'vpwp_ac', 1., 0],\
    ['vpwp_cl', False, 'vpwp_cl', 1., 0],\
    ['vpwp_ma', False, 'vpwp_ma', 1., 0],\
    ['vpwp_mfl', False, 'vpwp_mfl', 1., 0],\
    ['advection', True, 'vpwp_ta', 1., 0],\
    ['buoyancy', True, 'vpwp_bp', 1., 0],\
    ['dissipation', True, 'vpwp_dp1', 1., 0],\
    ['pressure', True, 'vpwp_pr1 + vpwp_pr2 + vpwp_pr3 + vpwp_pr4', 1., 0],\
    ['turb. prod.', True, 'vpwp_tp', 1. ,0],\
    ['time tndcy', True, 'vpwp_bt', 1., 0],\
    ['residual', True, 'vpwp_ac + vpwp_cl + vpwp_ma + vpwp_mfl', 1., 0],\
    ]

up2_reduced = [
    ['up2_dp1', False, 'up2_dp1', 1., 0],\
    ['up2_dp2', False, 'up2_dp2', 1., 0],\
    ['up2_pr1', False, 'up2_pr1', 1., 0],\
    ['up2_pr2', False, 'up2_pr2', 1., 0],\
    ['up2_splat', False, 'up2_splat', 1., 0],\
    ['up2_cl', False, 'up2_cl', 1., 0],\
    ['up2_ma', False, 'up2_ma', 1., 0],\
    ['up2_pd', False, 'up2_pd', 1., 0],\
    ['up2_sdmp', False, 'up2_sdmp', 1., 0],\
    ['up2_sf', False, 'up2_sf', 1., 0],\
    ['advection', True, 'up2_ta', 1., 0],\
    ['buoyancy', True, None, 1., 0],\
    ['dissipation', True, 'up2_pr1 + up2_dp2', 1., 0],\
    ['pressure', True, 'up2_dp1 + up2_pr2 + up2_splat', 1., 0],\
    ['turb. prod.', True, 'up2_tp', 1., 0],\
    ['time tndcy', True, 'up2_bt', 1., 0],\
    ['residual', True, 'up2_cl + up2_ma + up2_pd + up2_sdmp + up2_sf', 1., 0],\
    ]

vp2_reduced = [
    ['vp2_dp1', False, 'vp2_dp1', 1., 0],\
    ['vp2_dp2', False, 'vp2_dp2', 1., 0],\
    ['vp2_pr1', False, 'vp2_pr1', 1., 0],\
    ['vp2_pr2', False, 'vp2_pr2', 1., 0],\
    ['vp2_splat', False, 'vp2_splat', 1., 0],\
    ['vp2_cl', False, 'vp2_cl', 1., 0],\
    ['vp2_ma', False, 'vp2_ma', 1., 0],\
    ['vp2_pd', False, 'vp2_pd', 1., 0],\
    ['vp2_sdmp', False, 'vp2_sdmp', 1., 0],\
    ['vp2_sf', False, 'vp2_sf', 1., 0],\
    ['advection', True, 'vp2_ta', 1., 0],\
    ['buoyancy', True, None, 1., 0],\
    ['dissipation', True, 'vp2_pr1 + vp2_dp2', 1., 0],\
    ['pressure', True, 'vp2_dp1 + vp2_pr2 + vp2_splat', 1., 0],\
    ['turb. prod.', True, 'vp2_tp', 1., 0],\
    ['time tndcy', True, 'vp2_bt', 1., 0],\
    ['residual', True, 'vp2_cl + vp2_ma + vp2_pd + vp2_sdmp + vp2_sf', 1., 0],\
    ]

wp2_reduced = [
    ['wp2_dp1', False, 'wp2_dp1', 1., 0],\
    ['wp2_dp2', False, 'wp2_dp2', 1., 0],\
    ['wp2_pr1', False, 'wp2_pr1', 1., 0],\
    ['wp2_pr2', False, 'wp2_pr2', 1., 0],\
    ['wp2_pr3', False, 'wp2_pr3', 1., 0],\
    ['wp2_splat', False, 'wp2_splat', 1., 0],\
    ['wp2_ac', False, 'wp2_ac', 1., 0],\
    ['wp2_cl', False, 'wp2_cl', 1., 0],\
    ['wp2_ma', False, 'wp2_ma', 1., 0],\
    ['wp2_pd', False, 'wp2_pd', 1., 0],\
    ['wp2_sdmp', False, 'wp2_sdmp', 1., 0],\
    ['wp2_sf', False, 'wp2_sf', 1., 0],\
    ['advection', True, 'wp2_ta', 1., 0],\
    ['buoyancy', True, 'wp2_bp', 1., 0],\
    ['dissipation', True, 'wp2_dp1 + wp2_dp2', 1., 0],\
    ['pressure', True, 'wp2_pr1 + wp2_pr2 + wp2_pr3 + wp2_splat', 1., 0],\
    ['turb. prod.', True, None, 1., 0],\
    ['time tndcy', True, 'wp2_bt', 1., 0],\
    ['residual', True, 'wp2_ac + wp2_cl + wp2_ma + wp2_pd + wp2_sdmp + wp2_sf', 1., 0],\
    ]

# List of combined setup parameters for each single plot
# List entry description:
# index | description
#--------------------
#   0   | plot id
#   1   | plot title
#   2   | plot x-label
#   3   | plot text
#   4   | plot text position
#   5   | plot lines
plots_zm = [
    ['upwp_detailed', r"$\overline{u'w'}$ budget", r"$\overline{u'w'}$ budget terms $\mathrm{\left[m^2\,s^{-3} \right]}$", '', (.1,.9), upwp_detailed],
    ['vpwp_detailed', r"$\overline{v'w'}$ budget", r"$\overline{v'w'}$ budget terms $\mathrm{\left[m^2\,s^{-3} \right]}$", '', (.1,.9), vpwp_detailed],
    ['up2_detailed', r"$\overline{u'^2}$ budget", r"$\overline{u'^2}$ budget terms $\mathrm{\left[m^2\,s^{-3} \right]}$", '', (.1,.9), up2_detailed],
    ['vp2_detailed', r"$\overline{v'^2}$ budget", r"$\overline{v'^2}$ budget terms $\mathrm{\left[m^2\,s^{-3} \right]}$", '', (.1,.9), vp2_detailed],
    ['wp2_detailed', r"$\overline{w'^2}$ budget", r"$\overline{w'^2}$ budget terms $\mathrm{\left[m^2\,s^{-3} \right]}$", '', (.1,.9), wp2_detailed],
    ['upwp_reduced', r"$\overline{u'w'}$ budget", r"$\overline{u'w'}$ budget terms $\mathrm{\left[m^2\,s^{-3} \right]}$", '', (.1,.9), upwp_reduced],
    ['vpwp_reduced', r"$\overline{v'w'}$ budget", r"$\overline{v'w'}$ budget terms $\mathrm{\left[m^2\,s^{-3} \right]}$", '', (.1,.9), vpwp_reduced],
    ['up2_reduced', r"$\overline{u'^2}$ budget", r"$\overline{u'^2}$ budget terms $\mathrm{\left[m^2\,s^{-3} \right]}$", '', (.1,.9), up2_reduced],
    ['vp2_reduced', r"$\overline{v'^2}$ budget", r"$\overline{v'^2}$ budget terms $\mathrm{\left[m^2\,s^{-3} \right]}$", '', (.1,.9), vp2_reduced],
    ['wp2_reduced', r"$\overline{w'^2}$ budget", r"$\overline{w'^2}$ budget terms $\mathrm{\left[m^2\,s^{-3} \right]}$", '', (.1,.9), wp2_reduced]
    ]

plots_zt = []

plots = plots_zm + plots_zt