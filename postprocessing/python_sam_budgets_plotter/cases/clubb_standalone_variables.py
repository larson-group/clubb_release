"""
-------------------------------------------------------------------------------
    G E N E R A L   I N F O R M A T I O N
-------------------------------------------------------------------------------
        | CLUBB STANDALONE PLOTTING SETUP |
        -----------------------------------
This file contains general constants and information about the CLUBB variables saved
in the netCDF file needed for plotgen.py.

The "plots" array contains setup information about each individual plots,
while the each variable listed in plots lists the lines in that particular plot,
and the netCDF variables used to compute its data.

To add text to a plot, replace the string listed at index 3 with the text you want to appear,
and put the text coordinates as a tuple at index 4 in plots.
The point (0,0) is the origin of the plot, while the point (1,1) is at the upper right corner of the plot.

To generate new plots, add an entry to the plots variable,
create a new list of lines, and put it into index 5 of the new plots entry.

Note for CLUBB: The variables are divided between multiple netCDF files.
We are currently using 2: zm and zt
So newly added variables will have to go into the correct lists.
The user should be careful not to use the same plot name in sortPlots_zm and _zt,
as these will be used as keys in a dictionary.

TODO:   - Construct plot name from long name in netcdf instead
        - Add style entry to lines, 
        - Choice: Keep plot entries as list and use alias for indexing OR Change from lists to dicts OR use pandas
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
header = 'CLUBB standalone profiles'                        # Plot description used for the header on the html page (NOTE: redundant with name?)
name = 'clubb_standalone'                                   # Plot description used for file names and identifying plot types, MUST CONTAIN MODEL NAMES!
prefix = ''                                                 # Prefix identifier used in plot titles
nc_files = ['clubb_zm', 'clubb_zt']                         # List of NETCDF files containing the data needed for creating the plots listed belo. Paths are defined in case setup files

#-------------------------------------------------------------------------------
# P L O T S
#-------------------------------------------------------------------------------
# Define lines of each plot:
# Define plot line lists
# Line list entry description:
# index | description
#---------------------
#   0   | line name
#   1   | line visibility setting
#   2   | line variable or expression
#   3   | line scaling factor
#   4   | line xaxis selection

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

# List of combined setup parameters for each single plot
# List entry description:
# index | description
#--------------------
#   0   | plot id
#   1   | plot title
#   2   | plot x-label
#   3   | plot text
#   4   | plot text position
#   5   | plot lines
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