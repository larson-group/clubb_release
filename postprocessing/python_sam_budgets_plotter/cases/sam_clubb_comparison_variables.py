"""
-------------------------------------------------------------------------------
G E N E R A L   I N F O R M A T I O N
-------------------------------------------------------------------------------
This file contains general constants and information about the CLUBB and SAM variables saved
in the respective netCDF files needed for plotgen.py.

The list variables sortPlots, plotNames and lines are split up for CLUBB,
depending on which file they are contained in,
and sorted identically in order to relate the individual variables.
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
startLevel = 0                                              # Set the lower height level at which the plots should begin. For example, startLevel=2 would cut off the lowest 2 data points for each line.
header = 'SAM CLUBB comparison'
name = 'sam_clubb_comparison'                               # String used as part of the output file name
nc_files = ['clubb_zm', 'clubb_zt', 'sam']                  # NetCDF files needed for plots, paths are defined in case setup files

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