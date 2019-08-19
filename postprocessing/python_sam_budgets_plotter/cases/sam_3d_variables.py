"""
-------------------------------------------------------------------------------
G E N E R A L   I N F O R M A T I O N
-------------------------------------------------------------------------------
This file contains general constants and information about the SAM variables saved
in the netCDF file needed for plotgen.py.

The list variables sortPlots, plotNames and lines are sorted identically in
order to relate the individual variables.
"""
#-------------------------------------------------------------------------------
#   I M P O R T S
#-------------------------------------------------------------------------------
from numpy import nan

#-------------------------------------------------------------------------------
#   C O N S T A N T S
#-------------------------------------------------------------------------------
DAY = 24
HOUR = 3600
KG = 1000.
g_per_second_to_kg_per_day = 1. / (DAY * HOUR * KG)
kg_per_second_to_kg_per_day = 1. / (DAY * HOUR)
filler = nan
header = 'SAM horizontal plots'
name = 'sam_3d_{wt}'
nc_files = ['sam', 'sam_3d']
title_template = '{{wt}}, {case}, {x}x{y}x{z}, dx={dx:.0f}m, dz={dz:.1f}m, t={t:.0f}min, h={{h:.1f}}m'
wind_types = ['total_horizontal_wind+w_map', 'total_horizontal_wind+up_map', 'total_horizontal_wind+vp_map', 'total_horizontal_wind+uw_map', 'total_horizontal_wind+vw_map', 'horizontal_wind_perturbation+w_map', 'horizontal_wind_perturbation+up_map', 'horizontal_wind_perturbation+vp_map', 'horizontal_wind_perturbation+uw_map', 'horizontal_wind_perturbation+vw_map']

#-------------------------------------------------------------------------------
# P L O T   S E T U P
#-------------------------------------------------------------------------------
# Names of the variables
cloud_frac_cmap = 'Blues'
profiles = ['uw', 'vw']

#-------------------------------------------------------------------------------
# P L O T S
#-------------------------------------------------------------------------------
# Names of the variables (NOTE: THETAV NOT IN 3D DATA!)
#sortPlots_3d = ['qn', 'thetav', 'u', 'v', 'w']
sortPlots_3d = ['qn_3d', 'qv_3d', 'u_3d', 'v_3d', 'w_3d']
sortPlots_std = ['uw', 'vw', 'u', 'v', 'w', 'ucld', 'vcld', 'wcld']
# TODO: additional std data: U, V (,U2, V2, THV, UTHV, VTHV)
sortPlots = sortPlots_3d + sortPlots_std

qn_3d = [\
    ['QN', True, 'QN', 1., 0],\
    ]

qv_3d = [\
    ['QV', True, 'QV', 1., 0],\
    ]

#thetav = [\
    #['THETAV', True, 'THETAV', 1., 0],\
    #]

u_3d = [\
    ['U', True, 'U', 1., 0],\
    ]

v_3d = [\
    ['V', True, 'V', 1., 0],\
    ]

w_3d = [\
    ['W', True, 'W', 1., 0],\
    ]

uw = [\
    [r"$\overline{u'w'}$", True, 'UW', 1., 0],\
    ]

vw = [\
    [r"$\overline{v'w'}$", True, 'VW', 1., 0],\
    ]

u = [\
    # variables of U
    ['U', True, 'U', 1., 0],\
    ]

v = [\
    # variables of V
    ['V', True, 'V', 1., 0],\
    ]

w = [\
    # variables of W
    ['WM', True, 'WM', 1., 0],\
    ]

ucld = [\
    ['In-cloud mean of U', True, 'UCLD', 1., 0 ],\
    ]

vcld = [\
    ['In-cloud mean of V', True, 'VCLD', 1., 0 ],\
    ]

wcld = [\
    ['In-cloud mean of W', True, 'WCLD', 1., 0 ],\
    ]

#lines_3d = [qn, thetav, u, v, w]
lines_3d = [qn_3d, qv_3d, u_3d, v_3d, w_3d]
lines_std = [uw, vw, u, v, w, ucld, vcld, wcld]
lines = lines_3d + lines_std
