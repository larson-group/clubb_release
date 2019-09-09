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

TODO: Redo setup lists so noone has to count lines to get the number of entries right
Structure of list entries:
- plot name
- title
- units
- text + position
- list of variables and lines
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
#startLevel = 0                                              # Set the lower height level at which the plots should begin. For example, startLevel=2 would cut off the lowest 2 data points for each line. (DEPRECATED: Redundant with startHeight entry in case setup files)
header = 'CLUBB standalone profiles'                        # Plot description used for the header on the html page (NOTE: redundant with name?)
name = 'clubb_standalone'                                   # Plot description used for file names and identifying plot types, MUST CONTAIN MODEL NAMES!
prefix = ''                                                 # Prefix identifier used in plot titles
nc_files = ['clubb_zm', 'clubb_zt']                         # List of NETCDF files containing the data needed for creating the plots listed belo. Paths are defined in case setup files
# Put additional text entry into plot (TODO: Create lists for texts and positions for each plot)
plotText = ''                                               # Additional text entry to be put into plot
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

# settings of each plot:
# (plot number, )plot title, x-axis label
# TODO: Construct plot name from long name in netcdf instead
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

# Text and position (specified in axis coords: (0,0)=lower left corner, (1,1)=upper right) inserted into each plot)
text_zm = [\
    ['',(.1,.9)],\
    ['',(.1,.9)],\
    ['',(.1,.9)],\
    ['',(.1,.9)],\
    ['',(.1,.9)],\
    ['',(.1,.9)],\
    ['',(.1,.9)],\
    ['',(.1,.9)],\
    ]

text_zt = [\
    ['',(.1,.9)],\
    ['',(.1,.9)],\
    ]

text = text_zm + text_zt

# Define plot lists
# Line list entry description:
#0: line name
#1: line visibility setting
#2: line variable or expression
#3: line scaling factor
#4: line xaxis selection

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

# Combine all setup parameters into one list of lists/dicts(?)
# TODO: Add style entry to lines, 
# Choice: Keep plot entries as list and use alias for indexing OR Change from lists to dicts OR use pandas
# List entry description:
#0: plot id
#1: plot title
#2: plot x-label
#3: plot text
#4: plot text position
#5: plot lines
plots_zm = [
    ['uprcp', "Eastward liquid water flux", r"$\overline{u'r_c'}\ \mathrm{\left[kg\,kg^{-1}\,m\,s^{-1}\right]}$", '',(.1,.9), uprcp],
    ['vprcp', "Northward liquid water flux", r"$\overline{v'r_c'}\ \mathrm{\left[kg\,kg^{-1}\,m\,s^{-1}\right]}$", '',(.1,.9), vprcp],
    ['upthvp', "Eastward theta_v flux", r"$\overline{u'\theta_v'}\ \mathrm{\left[K\,m\,s^{-1}\right]}$", '',(.1,.9), upthvp],
    ['vpthvp', "Northward theta_v flux",r"$\overline{v'\theta_v'}\ \mathrm{\left[K\,m\,s^{-1}\right]}$", '',(.1,.9), vpthvp],
    ['uprtp', "Eastward total water flux", r"$\overline{u'r_t'}\ \mathrm{\left[kg\,kg^{-1}\,m\,s^{-1}\right]}$", '',(.1,.9), uprtp],
    ['vprtp', "Northward total water flux", r"$\overline{v'r_t'}\ \mathrm{\left[kg\,kg^{-1}\,m\,s^{-1}\right]}$", '',(.1,.9), vprtp],
    ['upthlp', "Eastward theta_l flux", r"$\overline{u'\theta_l'}\ \mathrm{\left[K\,m\,s^{-1}\right]}$", '',(.1,.9), upthlp],
    ['vpthlp', "Northward theta_l flux", r"$\overline{v'\theta_l'}\ \mathrm{\left[K\,m\,s^{-1}\right]}$", '',(.1,.9), vpthlp],
    ]

plots_zt = [
    ['um', 'East-west (u) wind', r"$\overline{u} \mathrm{\left[m\,s^{-1}\right]}$", '',(.1,.9), um],
    ['vm', 'North-south (v) wind', r"$\overline{v} \mathrm{\left[m\,s^{-1}\right]}$", '',(.1,.9), vm],
    ]

plots = plots_zm + plots_zt