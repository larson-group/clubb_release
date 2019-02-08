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
    ["Eastward liquid water flux", r"$\overline{u'r_c'}\ \left[\frac{m\ kg}{s\ kg}\right]$"],\
    ["Northward liquid water flux", r"$\overline{v'r_c'}\ \left[\frac{m\ kg}{s\ kg}\right]$"],\
    ["Eastward theta_v flux", r"$\overline{u'\theta_v'}\ \left[\frac{m\ K}{s}\right]$"],\
    ["Northward theta_v flux",r"$\overline{v'\theta_v'}\ \left[\frac{m\ K}{s}\right]$"],\
    ["Eastward total water flux", r"$\overline{u'r_t'}\ \left[\frac{m\ kg}{s\ kg}\right]$"],\
    ["Northward total water flux", r"$\overline{v'r_t'}\ \left[\frac{m\ kg}{s\ kg}\right]$"],\
    ["Eastward theta_l flux", r"$\overline{u'\theta_l'}\ \left[\frac{m\ K}{s}\right]$"],\
    ["Northward theta_l flux", r"$\overline{v'\theta_l'}\ \left[\frac{m\ K}{s}\right]$"],\
    ]

plotNames_zt = [\
    ['East-west (u) wind', r"$\bar{u} \left[\frac{m}{s}\right]$"],\
    ['North-south (v) wind', r"$\bar{v} \left[\frac{m}{s}\right]$"],\
    ]
plotNames = plotNames_zm + plotNames_zt

# Define plots
# zm

uprcp = [\
    [r"$\overline{u'r_c'}$", True, 'uprcp', 1., 0],\
    ]

vprcp = [\
    [r"$\overline{v'r_c'}$", True, 'vprcp', 1., 0],\
    ]

upthvp = [\
    [r"$\overline{u'\theta_v'}$", True, 'upthvp', 1., 0],\
    ]

vpthvp = [\
    [r"$\overline{v'\theta_v'}$", True, 'vpthvp',1., 0],\
    ]

uprtp = [\
    [r"$\overline{u'r_t'}$", True, 'uprtp', 1., 0],\
    ]

vprtp = [\
    [r"$\overline{v'r_t'}$", True, 'vprtp', 1., 0],\
    ]

upthlp = [\
    [r"$\overline{u'\theta_l'}$", True, 'upthlp', 1., 0],\
    ]

vpthlp = [\
    [r"$\overline{v'\theta_l'}$", True, 'vpthlp',1., 0],\
    ]

# zt

um = [\
    [r"$\bar{u}$", True, 'um', 1., 0],\
    ]

vm = [\
    [r"$\bar{v}$", True, 'vm', 1., 0],\
    ]

# Gather plots in list

lines_zm = [uprcp, vprcp, upthvp, vpthvp, uprtp, vprtp, upthlp, vpthlp]
lines_zt = [um, vm]
lines = lines_zm + lines_zt