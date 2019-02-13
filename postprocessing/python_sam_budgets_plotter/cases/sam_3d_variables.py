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
nc_files = ['sam_3d']
title_template = '{wt} {case} {x}x{y}x{z} dx={dx}m dz={dz}m t={t} height level {h}m'
wind_types = ['total_horizontal_wind', 'horizontal_wind_perturbation']

#-------------------------------------------------------------------------------
# P L O T   S E T U P
#-------------------------------------------------------------------------------
# Names of the variables
cloud_frac_cmap = 'Blues'
quiver_cmap = 'copper'


#-------------------------------------------------------------------------------
# P L O T S
#-------------------------------------------------------------------------------
# Names of the variables
sortPlots = ['qn', 'u', 'v', 'w']

qn = [\
    ['QN', True, 'QN', 1., 0],\
    ]

u = [\
    ['U', True, 'U', 1., 0],\
    ]

v = [\
    ['V', True, 'V', 1., 0],\
    ]

w = [\
    ['W', True, 'W', 1., 0],\
    ]

lines = [qn, u, v, w]