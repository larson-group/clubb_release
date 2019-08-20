"""
-------------------------------------------------------------------------------
G E N E R A L   I N F O R M A T I O N
-------------------------------------------------------------------------------
This file contains general constants and information about the CLUBB variables saved
in the netCDF file needed for plotgen.py.

The list variables sortPlots, plotNames and lines are split up,
depending on which file they are contained in,
and sorted identically in order to relate the individual variables.

The user should be careful not to use the same plot name in sortPlots_zm and _zt,
as these will be used as keys in a dictionary.
"""
#-------------------------------------------------------------------------------
#   I M P O R T S
#-------------------------------------------------------------------------------
from numpy import nan

#-------------------------------------------------------------------------------
#   C O N S T A N T S
#-------------------------------------------------------------------------------
DAY = 24                                                    # 1d  = 24h
HOUR = 3600                                                 # 1h  = 3600s
KG = 1000.                                                  # 1kg = 1000g
g_per_second_to_kg_per_day = 1. / (DAY * HOUR * KG)
kg_per_second_to_kg_per_day = 1. / (DAY * HOUR)
filler = nan                                                # Define the fill value which should replace invalid values in the data
startLevel = 0                                              # Set the lower height level at which the plots should begin. For example, startLevel=2 would cut off the lowest 2 data points for each line. (NOTE: Redundant with startHeight entry in case setup files)
header = 'CLUBB standalone profiles'
name = 'clubb_standalone'                                   # String used as part of the output file name
nc_files = ['clubb_zm', 'clubb_zt']                         # NetCDF files needed for plots, paths are defined
# Put additional text entry into plot (TODO: Create lists for texts and positions for each plot)
plotText = 'a)'                                             # Additional text entry to be put into plot
textPos = (.1,.9)                                           # Position of text within plot in data coordinates (x,y)

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
    ["Eastward liquid water flux", r"$\overline{u'r_c'}\ \mathrm{\left[kg\,kg^{-1}\,m\,s^{-1}\right]}$"],\
    ["Northward liquid water flux", r"$\overline{v'r_c'}\ \mathrm{\left[kg\,kg^{-1}\,m\,s^{-1}\right]}$"],\
    ["Eastward theta_v flux", r"$\overline{u'\theta_v'}\ \mathrm{\left[K\,m\,s^{-1}\right]}$"],\
    ["Northward theta_v flux",r"$\overline{v'\theta_v'}\ \mathrm{\left[K\,m\,s^{-1}\right]}$"],\
    ["Eastward total water flux", r"$\overline{u'r_t'}\ \mathrm{\left[kg\,kg^{-1}\,m\,s^{-1}\right]}$"],\
    ["Northward total water flux", r"$\overline{v'r_t'}\ \mathrm{\left[kg\,kg^{-1}\,m\,s^{-1}\right]}$"],\
    ["Eastward theta_l flux", r"$\overline{u'\theta_l'}\ \mathrm{\left[K\,m\,s^{-1}\right]}$"],\
    ["Northward theta_l flux", r"$\overline{v'\theta_l'}\ \mathrm{\left[K\,m\,s^{-1}\right]}$"],\
    ]

plotNames_zt = [\
    ['East-west (u) wind', r"$\overline{u} \mathrm{\left[m\,s^{-1}\right]}$"],\
    ['North-south (v) wind', r"$\overline{v} \mathrm{\left[m\,s^{-1}\right]}$"],\
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
    [r"$\overline{u}$", True, 'um', 1., 0],\
    ]

vm = [\
    [r"$\overline{v}$", True, 'vm', 1., 0],\
    ]

# Gather plots in list

lines_zm = [uprcp, vprcp, upthvp, vpthvp, uprtp, vprtp, upthlp, vpthlp]
lines_zt = [um, vm]
lines = lines_zm + lines_zt