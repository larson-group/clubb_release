"""
-------------------------------------------------------------------------------
G E N E R A L   I N F O R M A T I O N
-------------------------------------------------------------------------------
This file contains general constants and information about the CLUBB variables saved
in the netCDF file needed for plotgen.py.

The list variables sortPlots, plotNames and lines are sorted identically in
order to relate the individual variables.
"""
# TODO: How to store zm and zt parts?
#-------------------------------------------------------------------------------
#   C O N S T A N T S
#-------------------------------------------------------------------------------
DAY = 24
HOUR = 3600
KG = 1000.
g_per_second_to_kg_per_day = 1. / (DAY * HOUR * KG)
kg_per_second_to_kg_per_day = 1. / (DAY * HOUR)
header = 'CLUBB standalone profiles'
name = 'clubb_standalone'
nc_files = ['clubb_zm', 'clubb_zt']

#-------------------------------------------------------------------------------
# P L O T S
#-------------------------------------------------------------------------------
# Names of the variables
# zm
sortPlots_zm = ['uprcp', 'vprcp', 'upthvp', 'vpthvp', 'uprtp', 'vprtp', 'upthlp', 'vpthlp']
# zt
sortPlots_zt = ['um', 'vm']

sortPlots = sortPlots_zm + sortPlots_zt

# Construct plot name from long name in netcdf instead
plotNames_zm = [\
    ["Eastward liquid water flux", "u'rc' [(m/s)(kg/kg)]"],\
    ["Northward liquid water flux", "v'rc' [(m/s)(kg/kg)]"],\
    ["Eastward theta_v flux", "u'thv' [(m/s)K]"],\
    ["Northward theta_v flux","v'thv' [(m/s)K]"],\
    ["Eastward total water flux", "u'rt' [(m/s)(kg/kg)]"],\
    ["Northward total water flux", "v'rt' [(m/s)(kg/kg)]"],\
    ["Eastward theta_l flux", "u'thl' [(m/s)K]"],\
    ["Northward theta_l flux", "v'thl' [(m/s)K]"],\
    ]

plotNames_zt = [\
    ['East-west (u) wind', 'um [m/s]'],\
    ['North-south (v) wind', 'vm [m/s]'],\
    ]
plotNames = plotNames_zm + plotNames_zt

# Define plots
# zm

uprcp = [\
    ['uprcp', True, 'uprcp', 1., 0],\
    ]

vprcp = [\
    ['vprcp', True, 'vprcp', 1., 0],\
    ]

upthvp = [\
    ['upthvp', True, 'upthvp', 1., 0],\
    ]

vpthvp = [\
    ['vpthvp', True, 'vpthvp',1., 0],\
    ]

uprtp = [\
    ['uprtp', True, 'uprtp', 1., 0],\
    ]

vprtp = [\
    ['vprtp', True, 'vprtp', 1., 0],\
    ]

upthlp = [\
    ['upthlp', True, 'upthlp', 1., 0],\
    ]

vpthlp = [\
    ['vpthlp', True, 'vpthlp',1., 0],\
    ]

# zt

um = [\
    ['um', True, 'um', 1., 0],\
    ]

vm = [\
    ['vm', True, 'vm', 1., 0],\
    ]

# Gather plots in list

lines_zm = [uprcp, vprcp, upthvp, vpthvp, uprtp, vprtp, upthlp, vpthlp]
lines_zt = [um, vm]
lines = lines_zm + lines_zt