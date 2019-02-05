#-------------------------------------------------------------------------------
#   C O N S T A N T S
#-------------------------------------------------------------------------------
DAY = 24
HOUR = 3600
KG = 1000
g_per_second_to_kg_per_day = 1. / (DAY * HOUR * KG)
kg_per_second_to_kg_per_day = 1. / (DAY * HOUR)
header = 'CLUBB budgets'
name = 'clubb_budgets'
nc_files = ['clubb_zm']

#-------------------------------------------------------------------------------
# P L O T S
#-------------------------------------------------------------------------------
# Names of variables
sortPlots_zm = ['upwp_detailed', 'vpwp_detailed', 'up2_detailed', 'vp2_detailed', 'wp2_detailed', 'upwp_reduced', 'vpwp_reduced', 'up2_reduced', 'vp2_reduced', 'wp2_reduced']
sortPlots_zt = []
sortPlots = sortPlots_zm + sortPlots_zt

# Construct plot name from long name in netcdf instead
plotNames_zm = [\
    ["Vertical east-west momentum flux", "u'w' [m^2/s^2]"],\
    ["Vertical north-south momentum flux", "v'w' [m^2/s^2]"],\
    ["Variance of east-west air velocity", "u'^2[m^2/s^2]"],\
    ["Variance of north-south air velocity", "v'^2 [m^2/s^2]"],\
    ["Variance of vertical air velocity", "w'^2 [m^2/s^2]"],\
    ["Vertical east-west momentum flux", "u'w' [m^2/s^2]"],\
    ["Vertical north-south momentum flux", "v'w' [m^2/s^2]"],\
    ["Variance of east-west air velocity", "u'^2[m^2/s^2]"],\
    ["Variance of north-south air velocity", "v'^2 [m^2/s^2]"],\
    ["Variance of vertical air velocity", "w'^2 [m^2/s^2]"],\
        ]

plotNames_zt = []

plotNames = plotNames_zm + plotNames_zt

# Define plot lists

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
    ]

# Define plots with reduced complexity
upwp_reduced = [\
    ['upwp_bt', False, 'upwp_bt', 1., 0],\
    ['upwp_ma', False, 'upwp_ma', 1., 0],\
    ['upwp_ta', True, 'upwp_ta', 1., 0],\
    ['upwp_tp', True, 'upwp_tp', 1. ,0],\
    ['upwp_ac', True, 'upwp_ac', 1., 0],\
    ['upwp_bp', True, 'upwp_bp', 1., 0],\
    ['upwp_pr1', False, 'upwp_pr1', 1., 0],\
    ['upwp_pr2', False, 'upwp_pr2', 1., 0],\
    ['upwp_pr3', False, 'upwp_pr3', 1., 0],\
    ['upwp_pr4', False, 'upwp_pr4', 1., 0],\
    ['upwp_dp1', False, 'upwp_dp1', 1., 0],\
    ['upwp_mfl', False, 'upwp_mfl', 1., 0],\
    ['upwp_cl', False, 'upwp_cl', 1., 0],\
    ['upwp_pres', True, 'upwp_pr1 + upwp_pr2 + upwp_pr3 + upwp_pr4', 1., 0],\
    ['upwp_diss', True, 'upwp_dp1', 1., 0],\
    ['upwp_res', True, 'upwp_bt + upwp_ma + upwp_cl + upwp_mfl', 1., 0],\
    ]

vpwp_reduced = [\
    ['vpwp_bt', False, 'vpwp_bt', 1., 0],\
    ['vpwp_ma', False, 'vpwp_ma', 1., 0],\
    ['vpwp_ta', True, 'vpwp_ta', 1., 0],\
    ['vpwp_tp', True, 'vpwp_tp', 1. ,0],\
    ['vpwp_ac', True, 'vpwp_ac', 1., 0],\
    ['vpwp_bp', True, 'vpwp_bp', 1., 0],\
    ['vpwp_pr1', False, 'vpwp_pr1', 1., 0],\
    ['vpwp_pr2', False, 'vpwp_pr2', 1., 0],\
    ['vpwp_pr3', False, 'vpwp_pr3', 1., 0],\
    ['vpwp_pr4', False, 'vpwp_pr4', 1., 0],\
    ['vpwp_dp1', False, 'vpwp_dp1', 1., 0],\
    ['vpwp_mfl', False, 'vpwp_mfl', 1., 0],\
    ['vpwp_cl', False, 'vpwp_cl', 1., 0],\
    ['vpwp_pres', True, 'vpwp_pr1 + vpwp_pr2 + vpwp_pr3 + vpwp_pr4', 1., 0],\
    ['vpwp_diss', True, 'vpwp_dp1', 1., 0],\
    ['vpwp_res', True, 'vpwp_bt + vpwp_ma + vpwp_cl + vpwp_mfl', 1., 0],\
    ]

up2_reduced = [
    ['up2_bt', False, 'up2_bt', 1., 0],\
    ['up2_ma', False, 'up2_ma', 1., 0],\
    ['up2_ta', True, 'up2_ta', 1., 0],\
    ['up2_tp', True, 'up2_tp', 1., 0],\
    ['up2_dp1', False, 'up2_dp1', 1., 0],\
    ['up2_dp2', False, 'up2_dp2', 1., 0],\
    ['up2_pr1', False, 'up2_pr1', 1., 0],\
    ['up2_pr2', False, 'up2_pr2', 1., 0],\
    ['up2_sdmp', False, 'up2_sdmp', 1., 0],\
    ['up2_cl', False, 'up2_cl', 1., 0],\
    ['up2_pd', False, 'up2_pd', 1., 0],\
    ['up2_sf', False, 'up2_sf', 1., 0],\
    ['up2_splat', True, 'up2_splat', 1., 0],\
    ['up2_pres', True, 'up2_dp1 + up2_pr2', 1., 0],\
    ['up2_diss', True, 'up2_dp2 + up2_pr1', 1., 0],\
    ['up2_res', True, 'up2_sf + up2_pd + up2_bt + up2_ma + up2_sdmp + up2_cl', 1., 0],\
    ]

vp2_reduced = [
    ['vp2_bt', False, 'vp2_bt', 1., 0],\
    ['vp2_ma', False, 'vp2_ma', 1., 0],\
    ['vp2_ta', True, 'vp2_ta', 1., 0],\
    ['vp2_tp', True, 'vp2_tp', 1., 0],\
    ['vp2_dp1', False, 'vp2_dp1', 1., 0],\
    ['vp2_dp2', False, 'vp2_dp2', 1., 0],\
    ['vp2_pr1', False, 'vp2_pr1', 1., 0],\
    ['vp2_pr2', False, 'vp2_pr2', 1., 0],\
    ['vp2_sdmp', False, 'vp2_sdmp', 1., 0],\
    ['vp2_cl', False, 'vp2_cl', 1., 0],\
    ['vp2_pd', False, 'vp2_pd', 1., 0],\
    ['vp2_sf', False, 'vp2_sf', 1., 0],\
    ['vp2_splat', True, 'vp2_splat', 1., 0],\
    ['vp2_pres', True, 'vp2_dp1 + vp2_pr2', 1., 0],\
    ['vp2_diss', True, 'vp2_dp2 + vp2_pr1', 1., 0],\
    ['vp2_res', True, 'vp2_sf + vp2_pd + vp2_bt + vp2_ma + vp2_sdmp + vp2_cl', 1., 0],\
    ]

wp2_reduced = [
    ['wp2_bt', False, 'wp2_bt', 1., 0],\
    ['wp2_ma', False, 'wp2_ma', 1., 0],\
    ['wp2_ta', True, 'wp2_ta', 1., 0],\
    ['wp2_ac', True, 'wp2_ac', 1., 0],\
    ['wp2_bp', True, 'wp2_bp', 1., 0],\
    ['wp2_dp1', False, 'wp2_dp1', 1., 0],\
    ['wp2_dp2', False, 'wp2_dp2', 1., 0],\
    ['wp2_pr1', False, 'wp2_pr1', 1., 0],\
    ['wp2_pr2', False, 'wp2_pr2', 1., 0],\
    ['wp2_pr3', False, 'wp2_pr3', 1., 0],\
    ['wp2_sdmp', False, 'wp2_sdmp', 1., 0],\
    ['wp2_cl', False, 'wp2_cl', 1., 0],\
    ['wp2_pd', False, 'wp2_pd', 1., 0],\
    ['wp2_sf', False, 'wp2_sf', 1., 0],\
    ['wp2_splat', True, 'wp2_splat', 1., 0],\
    ['wp2_pres', True, 'wp2_pr1 + wp2_pr2 + wp2_pr3', 1., 0],\
    ['wp2_diss', True, 'wp2_dp1 + wp2_dp2', 1., 0],\
    ['wp2_res', True, 'wp2_sf + wp2_pd + wp2_bt + wp2_ma + wp2_sdmp + wp2_cl', 1., 0],\
    ]

#lines_zm = [upwp, vpwp, up2, vp2, wp2]
lines_zm = [upwp_detailed, vpwp_detailed, up2_detailed, vp2_detailed, wp2_detailed, upwp_reduced, vpwp_reduced, up2_reduced, vp2_reduced, wp2_reduced]
lines_zt = []
lines = lines_zm + lines_zt