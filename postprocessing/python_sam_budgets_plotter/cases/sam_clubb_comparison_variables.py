"""
-------------------------------------------------------------------------------
G E N E R A L   I N F O R M A T I O N
-------------------------------------------------------------------------------
This file contains general constants and information about the CLUBB and SAM variables saved
in the respective netCDF files needed for plotgen.py.

The list variables sortPlots, plotNames and lines are split up for CLUBB,
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
TODO: Improve data structures, as there is lots of duplicate data in here:
Needed: 1x line id, 2x var/expression, 2x conversion, 1x visibility flag, 0-1x x-axis selection
There should never be more than one varibale plotted per plot.
Does that make creating appropriate data structures easier?
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
header = 'SAM CLUBB comparison'                             # Plot description used for the header on the html page (NOTE: redundant with name?)
name = 'sam_clubb_comparison'                               # Plot description used for file names and identifying plot types, MUST CONTAIN MODEL NAMES!
prefix = ''                                                 # Prefix identifier used in plot titles
nc_files = ['clubb_zm', 'clubb_zt', 'sam']                  # List of NETCDF files containing the data needed for creating the plots listed belo. Paths are defined in case setup files
# Put additional text entry into plot (TODO: Create lists for texts and positions for each plot)
plotText = ''                                               # Additional text entry to be put into plot
textPos = (.93,.9)                                          # Position of text within plot in data coordinates (x,y)

#-------------------------------------------------------------------------------
# P L O T S
#-------------------------------------------------------------------------------
# Names of the variables
# zm
sortPlots_zm = ['upwp', 'vpwp', 'up2', 'vp2', 'wp2', 'uprcp', 'vprcp', 'upthvp', 'vpthvp', 'uprtp', 'vprtp', 'upthlp', 'vpthlp']
# zt
sortPlots_zt = ['um', 'vm', 'cld', 'theta_l', 'r_t']

sortPlots = sortPlots_zm + sortPlots_zt

# settings of each plot:
# (plot number, )plot title, x-axis label
# TODO: Construct plot name from long name in netcdf instead
#plotNames_zm = [\
    #["Vertical eastward momentum flux", r"$\mathrm{\overline{u'w'}\ \left[\frac{m^2}{s^2}\right]}$"],\
    #["Vertical northward momentum flux", r"$\mathrm{\overline{v'w'}\ \left[\frac{m^2}{s^2}\right]}$"],\
    #["Variance of eastward air velocity", r"$\mathrm{\overline{u'^2}\ \left[\frac{m^2}{s^2}\right]}$"],\
    #["Variance of northward air velocity", r"$\mathrm{\overline{v'^2}\ \left[\frac{m^2}{s^2}\right]}$"],\
    #["Variance of vertical air velocity", r"$\mathrm{\overline{w'^2}\ \left[\frac{m^2}{s^2}\right]}$"],\
    #["Eastward liquid water flux", r"$\mathrm{\overline{u'r_c'}\ \left[\frac{m\ kg}{s\ kg}\right]}$"],\
    #["Northward liquid water flux", r"$\mathrm{\overline{v'r_c'}\ \left[\frac{m\ kg}{s\ kg}\right]}$"],\
    #["Eastward theta_v flux", r"$\mathrm{\overline{u'\theta_v'}\ \left[\frac{m\ K}{s}\right]}$"],\
    #["Northward theta_v flux", r"$\mathrm{\overline{v'\theta_v'}\ \left[\frac{m\ K}{s}\right]}$"],\
    #["Eastward total water flux", r"$\mathrm{\overline{u'r_t'}\ \left[\frac{m\ kg}{s\ kg}\right]}$"],\
    #["Northward total water flux", r"$\mathrm{\overline{v'r_t'}\ \left[\frac{m\ kg}{s\ kg}\right]}$"],\
    #["Eastward theta_l flux", r"$\mathrm{\overline{u'\theta_l'}\ \left[\frac{m\ K}{s}\right]}$"],\
    #["Northward theta_l flux", r"$\mathrm{\overline{v'\theta_l'}\ \left[\frac{m\ K}{s}\right]}$"],\
    #]


plotNames_zm = [\
    [r"$\overline{u'w'}$", r"Momentum flux, $\overline{u'w'}\ \mathrm{\left[m^2\,s^{-2}\right]}$"],\
    [r"$\overline{v'w'}$", r"Momentum flux, $\overline{v'w'}\ \mathrm{\left[m^2\,s^{-2}\right]}$"],\
    [r"$\overline{u'^2}$", r"Momentum variance, $\overline{u'^2}\ \mathrm{\left[m^2\,s^{-2}\right]}$"],\
    [r"$\overline{v'^2}$", r"Momentum variance, $\overline{v'^2}\ \mathrm{\left[m^2\,s^{-2}\right]}$"],\
    [r"$\overline{w'^2}$", r"Momentum variance, $\overline{w'^2}\ \mathrm{\left[m^2\,s^{-2}\right]}$"],\
    [r"$\overline{u'r_c'}$", r"Liquid water flux, $\overline{u'r_c'}\ \mathrm{\left[kg\,kg^{-1}\,m\,s^{-1}\right]}$"],\
    [r"$\overline{v'r_c'}$", r"Liquid water flux, $\overline{v'r_c'}\ \mathrm{\left[kg\,kg^{-1}\,m\,s^{-1}\right]}$"],\
    [r"$\overline{u'\theta_v'}$", r"Virt. pot. temp. flux, $\overline{u'\theta_v'}\ \mathrm{\left[K\,m\,s^{-1}\right]}$"],\
    [r"$\overline{v'\theta_v'}$", r"Virt. pot. temp. flux, $\overline{v'\theta_v'}\ \mathrm{\left[K\,m\,s^{-1}\right]}$"],\
    [r"$\overline{u'r_t'}$", r"Total water flux, $\overline{u'r_t'}\ \mathrm{\left[kg\,kg^{-1}\,m\,s^{-1}\right]}$"],\
    [r"$\overline{v'r_t'}$", r"Total water flux, $\overline{v'r_t'}\ \mathrm{\left[kg\,kg^{-1}\,m\,s^{-1}\right]}$"],\
    [r"$\overline{u'\theta_l'}$", r"Liq. water pot. temp. flux, $\overline{u'\theta_l'}\ \mathrm{\left[K\,m\,s^{-1}\right]}$"],\
    [r"$\overline{v'\theta_l'}$", r"Liq. water pot. temp. flux, $\overline{v'\theta_l'}\ \mathrm{\left[K\,m\,s^{-1}\right]}$"],\
    ]

plotNames_zt = [\
    [r"$\overline{u}}$", r"Eastward mean wind, $\overline{u}\ \mathrm{\left[m\,s^{-1}\right]}$"],\
    [r"$\overline{v}}$", r"Northward mean wind, $\overline{v}\ \mathrm{\left[m\,s^{-1}\right]}$"],\
    [r"Cloud fraction", r"Cloud fraction [-]"],\
    [r"$\theta_l$", r"Liquid water potential temperature, $\theta_l\ \mathrm{\left[K\right]}$"],\
    [r"$r_t$", r"Total water mixing ratio, $r_t\ \mathrm{\left[kg\,kg^{-1}\right]}$"],\
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
    ['',(.1,.9)],\
    ['',(.1,.9)],\
    ['',(.1,.9)],\
    ['',(.1,.9)],\
    ['',(.1,.9)],\
    ]

text_zt = [\
    ['',(.1,.9)],\
    ['',(.1,.9)],\
    ['',(.1,.9)],\
    ['',(.1,.9)],\
    ['',(.1,.9)],\
    ]

# Define plot lists
# Line list entry description:
#0: line name
#1: line visibility setting
#2: line variable or expression
#3: line scaling factor
#4: line xaxis selection

##CLUBB
# zt

um_clubb = [\
    ['um', True, 'um', 1., 0],\
    ]

vm_clubb = [\
    ['vm', True, 'vm', 1., 0],\
    ]

cld_clubb = [\
    ['cloud_frac', True, 'cloud_frac', 1., 0],\
    ]

thetal_clubb = [\
    ['thlm', True, 'thlm', 1., 0],\
    ]

rt_clubb = [\
    ['rtm', True, 'rtm', 1., 0],\
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
lines_clubb_zm = [upwp_clubb, vpwp_clubb, up2_clubb, vp2_clubb, wp2_clubb, uprcp_clubb, vprcp_clubb, upthvp_clubb, vpthvp_clubb, uprtp_clubb, vprtp_clubb, upthlp_clubb, vpthlp_clubb]
lines_clubb_zt = [um_clubb, vm_clubb, cld_clubb, thetal_clubb, rt_clubb]
lines_clubb = lines_clubb_zm + lines_clubb_zt

##SAM
# Define plot lists
# Line list entry description:
#0: line name
#1: line visibility setting
#2: line variable or expression
#3: line scaling factor
#4: line xaxis selection
um_sam = [\
    ['U', True, 'U', 1., 0],\
    ]

vm_sam = [\
    ['V', True, 'V', 1., 0],\
    ]

cld_sam = [\
    ['CLD', True, 'CLD', 1., 0],\
    ]

thetal_sam = [\
    ['THETAL', True, 'THETAL', 1., 0],\
    ]

rt_sam = [\
    # variables of rt
    ['RT', True, '(QT-QI)/1000.', 1., 0],\
    ['QI', False, 'QI', 1., 0],\
    ['QT', False, 'QT', 1., 0],\
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

# Gather plots in list
lines_sam_zm = [upwp_sam, vpwp_sam, up2_sam, vp2_sam, wp2_sam, uprcp_sam, vprcp_sam, upthvp_sam, vpthvp_sam, uprtp_sam, vprtp_sam, upthlp_sam, vpthlp_sam]
lines_sam_zt = [um_sam, vm_sam, cld_sam, thetal_sam, rt_sam]
lines_sam = lines_sam_zm + lines_sam_zt

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
    ['upwp', r"$\overline{u'w'}$", r"Momentum flux, $\overline{u'w'}\ \mathrm{\left[m^2\,s^{-2}\right]}$", '', (0.1, 0.9), upwp_sam, upwp_clubb],
    ['vpwp', r"$\overline{v'w'}$", r"Momentum flux, $\overline{v'w'}\ \mathrm{\left[m^2\,s^{-2}\right]}$", '', (0.1, 0.9), vpwp_sam, vpwp_clubb],
    ['up2', r"$\overline{u'^2}$", r"Momentum variance, $\overline{u'^2}\ \mathrm{\left[m^2\,s^{-2}\right]}$", '', (0.1, 0.9), up2_sam, up2_clubb],
    ['vp2', r"$\overline{v'^2}$", r"Momentum variance, $\overline{v'^2}\ \mathrm{\left[m^2\,s^{-2}\right]}$", '', (0.1, 0.9), vp2_sam, vp2_clubb],
    ['wp2', r"$\overline{w'^2}$", r"Momentum variance, $\overline{w'^2}\ \mathrm{\left[m^2\,s^{-2}\right]}$", '', (0.1, 0.9), wp2_sam, wp2_clubb],
    ['uprcp', r"$\overline{u'r_c'}$", r"Liquid water flux, $\overline{u'r_c'}\ \mathrm{\left[kg\,kg^{-1}\,m\,s^{-1}\right]}$", '', (0.1, 0.9), uprcp_sam, uprcp_clubb],
    ['vprcp', r"$\overline{v'r_c'}$", r"Liquid water flux, $\overline{v'r_c'}\ \mathrm{\left[kg\,kg^{-1}\,m\,s^{-1}\right]}$", '', (0.1, 0.9), vprcp_sam, vprcp_clubb],
    ['upthvp', r"$\overline{u'\theta_v'}$", r"Virt. pot. temp. flux, $\overline{u'\theta_v'}\ \mathrm{\left[K\,m\,s^{-1}\right]}$", '', (0.1, 0.9), upthvp_sam, upthvp_clubb],
    ['vpthvp', r"$\overline{v'\theta_v'}$", r"Virt. pot. temp. flux, $\overline{v'\theta_v'}\ \mathrm{\left[K\,m\,s^{-1}\right]}$", '', (0.1, 0.9), vpthvp_sam, vpthvp_clubb],
    ['uprtp', r"$\overline{u'r_t'}$", r"Total water flux, $\overline{u'r_t'}\ \mathrm{\left[kg\,kg^{-1}\,m\,s^{-1}\right]}$", '', (0.1, 0.9), uprtp_sam, uprtp_clubb],
    ['vprtp', r"$\overline{v'r_t'}$", r"Total water flux, $\overline{v'r_t'}\ \mathrm{\left[kg\,kg^{-1}\,m\,s^{-1}\right]}$", '', (0.1, 0.9), vprtp_sam, vprtp_clubb],
    ['upthlp', r"$\overline{u'\theta_l'}$", r"Liq. water pot. temp. flux, $\overline{u'\theta_l'}\ \mathrm{\left[K\,m\,s^{-1}\right]}$", '', (0.1, 0.9), upthlp_sam, upthlp_clubb],
    ['vpthlp', r"$\overline{v'\theta_l'}$", r"Liq. water pot. temp. flux, $\overline{v'\theta_l'}\ \mathrm{\left[K\,m\,s^{-1}\right]}$", '', (0.1, 0.9), vpthlp_sam, vpthlp_clubb],
    ]

plots_zt = [
    ['um', r"$\overline{u}}$", r"Eastward mean wind, $\overline{u}\ \mathrm{\left[m\,s^{-1}\right]}$", '', (0.1, 0.9), um_sam, um_clubb],
    ['vm', r"$\overline{v}}$", r"Northward mean wind, $\overline{v}\ \mathrm{\left[m\,s^{-1}\right]}$", '', (0.1, 0.9), vm_sam, vm_clubb],
    ['cld', r"Cloud fraction", r"Cloud fraction [-]", '', (0.1, 0.9), cld_sam, cld_clubb],
    ['theta_l', r"$\theta_l$", r"Liquid water potential temperature, $\theta_l\ \mathrm{\left[K\right]}$", '', (0.1, 0.9), thetal_sam, thetal_clubb],
    ['r_t', r"$r_t$", r"Total water mixing ratio, $r_t\ \mathrm{\left[kg\,kg^{-1}\right]}$", '', (0.1, 0.9), rt_sam, rt_clubb],
    ]

# TODO: Split up SAM/CLUBB, zm/zt or combine into one structure??
#tmp = plots_zm+plots_zt

#plots_sam = [tmp[i]+[lines[i]] for i in range(len(lines_sam))]
#tmp = [tmp[i]+[lines[i]] for i in range(len(lines_sam))]

#plots_clubb = [ plots_zm[i]+[lines_zm[i]] for i in range(len(plots_zm))] + [ plots_zt[i]+[lines_zt[i]] for i in range(len(plots_zt))]

# Append line definitions for SAM and CLUBB onto each entry, then combine zm and zt lists
plots = plots_zm + plots_zt