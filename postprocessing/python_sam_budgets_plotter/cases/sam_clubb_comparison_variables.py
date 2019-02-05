"""
-------------------------------------------------------------------------------
G E N E R A L   I N F O R M A T I O N
-------------------------------------------------------------------------------
This file contains general constants and information about the CLUBB variables saved
in the netCDF file needed for plotgen.py.

The list variables sortPlots, plotNames and lines are sorted identically in
order to relate the individual variables.
"""

#-------------------------------------------------------------------------------
#   C O N S T A N T S
#-------------------------------------------------------------------------------
DAY = 24
HOUR = 3600
KG = 1000.
g_per_second_to_kg_per_day = 1. / (DAY * HOUR * KG)
kg_per_second_to_kg_per_day = 1. / (DAY * HOUR)
header = 'SAM CLUBB comparison'
name = 'sam_clubb_comparison'
nc_files = ['clubb_zm', 'clubb_zt', 'sam']

#-------------------------------------------------------------------------------
# P L O T S
#-------------------------------------------------------------------------------
# Names of the variables
# zm
sortPlots_zm = ['upwp', 'vpwp', 'up2', 'vp2', 'wp2', 'uprcp', 'vprcp', 'upthvp', 'vpthvp', 'uprtp', 'vprtp', 'upthlp', 'vpthlp']
# zt
sortPlots_zt = ['um', 'vm']

sortPlots = sortPlots_zm + sortPlots_zt

# Construct plot name from long name in netcdf instead
plotNames_zm = [\
    ["Vertical east-west momentum flux", "u'w' [m^2/s^2]"],\
    ["Vertical north-south momentum flux", "v'w' [m^2/s^2]"],\
    ["Variance of east-west air velocity", "u'^2[m^2/s^2]"],\
    ["Variance of north-south air velocity", "v'^2 [m^2/s^2]"],\
    ["Variance of vertical air velocity", "w'^2 [m^2/s^2]"],\
    ["Eastward liquid water flux", "u'rc' [(m/s)(kg/kg)]"],\
    ["Northward liquid water flux", "v'rc' [(m/s)(kg/kg)]"],\
    ["Eastward theta_v flux", "u'thv' [(m/s)K]"],\
    ["Northward theta_v flux", "v'thv' [(m/s)K]"],\
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
##CLUBB
# zt

um_clubb = [\
    ['um', True, 'um', 1., 0],\
    ]

vm_clubb = [\
    ['vm', True, 'vm', 1., 0],\
    ]

# zm

upwp_clubb = [\
    ['upwp', True, 'upwp', 1., 0],\
    ]

vpwp_clubb = [\
    ['vpwp', True, 'vpwp', 1., 0],\
    ]

up2_clubb = [\
    ['up2', True, 'up2', 1., 0],\
    ]

vp2_clubb = [\
    ['vp2', True, 'vp2', 1., 0],\
    ]

wp2_clubb = [\
    ['wp2', True, 'wp2', 1., 0],\
    ]

uprcp_clubb = [\
    ['uprcp', True, 'uprcp', 1., 0],\
    ]

vprcp_clubb = [\
    ['vprcp', True, 'vprcp', 1., 0],\
    ]

upthvp_clubb = [\
    ['upthvp', True, 'upthvp', 1., 0],\
    ]

vpthvp_clubb = [\
    ['vpthvp', True, 'vpthvp',1., 0],\
    ]

uprtp_clubb = [\
    ['uprtp', True, 'uprtp', 1., 0],\
    ]

vprtp_clubb = [\
    ['vprtp', True, 'vprtp', 1., 0],\
    ]

upthlp_clubb = [\
    ['upthlp', True, 'upthlp', 1., 0],\
    ]

vpthlp_clubb = [\
    ['vpthlp', True, 'vpthlp',1., 0],\
    ]

# Gather plots in list
lines_zm = [upwp_clubb, vpwp_clubb, up2_clubb, vp2_clubb, wp2_clubb, uprcp_clubb, vprcp_clubb, upthvp_clubb, vpthvp_clubb, uprtp_clubb, vprtp_clubb, upthlp_clubb, vpthlp_clubb]
lines_zt = [um_clubb, vm_clubb]
lines_clubb = lines_zm + lines_zt

##SAM
um_sam = [\
    ['U', True, 'U', 1., 0],\
    ]

vm_sam = [\
    ['V', True, 'V', 1., 0],\
    ]

upwp_sam = [\
    ['UW', True, 'UW', 1., 0],\
    ]

vpwp_sam = [\
    ['VW', True, 'VW', 1., 0],\
    ]

up2_sam = [\
    ['U2', True, 'U2', 1., 0],\
    ]

vp2_sam = [\
    ['V2', True, 'V2', 1., 0],\
    ]

wp2_sam = [\
    ['W2', True, 'W2', 1., 0],\
    ]

uprcp_sam = [\
    ['UPRCP', True, 'UPRCP', 1., 0],\
    ]

vprcp_sam = [\
    ['VPRCP', True, 'VPRCP', 1., 0],\
    ]

upthvp_sam = [\
    ['UPTHVP', True, 'UPTHVP', 1., 0],\
    ]

vpthvp_sam = [\
    ['VPTHVP', True, 'VPTHVP',1., 0],\
    ]

uprtp_sam = [\
    ['UPRTP', True, 'UPRTP', 1., 0],\
    ]

vprtp_sam = [\
    ['VPRTP', True, 'VPRTP', 1., 0],\
    ]

upthlp_sam = [\
    ['UPTHLP', True, 'UPTHLP', 1., 0],\
    ]

vpthlp_sam = [\
    ['VPTHLP', True, 'VPTHLP',1., 0],\
    ]

lines_sam = [upwp_sam, vpwp_sam, up2_sam, vp2_sam, wp2_sam, uprcp_sam, vprcp_sam, upthvp_sam, vpthvp_sam, uprtp_sam, vprtp_sam, upthlp_sam, vpthlp_sam, um_sam, vm_sam]