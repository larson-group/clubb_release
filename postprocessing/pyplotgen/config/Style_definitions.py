"""
:author: Nicolas Strike
:date: Late 2019

Style_definitions.py contains the style formatting parameters used by Pyplotgen throughout the source code.
If you wish to change a line color, font size, etc., then this is the place to do it.
"""

# Line Styles
# --------------

# Color rotation
# Multiple color options are given, uncomment the one that's most useful and comment out the rest. Feel free to make your own as well.
# Pyplot requires that style rotation's and color rotation's have the same length, so be sure to properly (un)comment the style rotations according
# to the desired color scheme.
# Color/style rotations are used on budget plots and for input-folder lines

# Colors are taken from the colorbrewer qualitative 10-color set | http://colorbrewer2.org/#type=qualitative&scheme=Paired&n=10
# Default option. Not colorblind or print safe.
COLOR_ROTATION = ['#1f78b4', '#e31a1c', '#ff7f00', '#6a3d9a', '#33a02c', '#a6cee3', '#b2df8a', '#fb9a99', '#fdbf6f', '#cab2d6']
STYLE_ROTATION = ['-', '--', ':', '-.','-', '--', ':', '-.','-', '--']

# Color-blind and print safe color rotation. These requirements only allow up to 4 colors
# COLOR_ROTATION = ['#1f78b4', '#33a02c', '#a6cee3', '#b2df8a']
# STYLE_ROTATION = ['-', '--', ':', '-.']

# Print safe colors (not color-blind friendly). From colorbrewer (removed yellow)
# COLOR_ROTATION = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#a65628', '#f781bf', '#999999']
# STYLE_ROTATION = ['-', '--', ':', '-.', '-', '--', ':', '-.']

# Details on line styles here: https://matplotlib.org/2.1.2/api/_as_gen/matplotlib.pyplot.plot.html
LES_LINE_STYLE = "k-"
GOLAZ_BEST_R408_LINE_STYLE = "g-"
HOC_LINE_STYLE = 'm-.'

# Line Thicknesses
CLUBB_LINE_THICKNESS = 2
LES_LINE_THICKNESS = 6
ARCHIVED_CLUBB_LINE_THICKNESS = 4
THIN_LINE_THICKNESS = 1  # Used if --thin-lines is passed, overrides every other size

# Line Labels
# This overrides using foldernames for line input folder line labels. Must be a list of strings.
# Labels must be in the same order as the input folders. To use foldernames, remove all elements from the array.
CLUBB_LABEL_OVERRIDE = [] # To override names, uncomment the line below
# CLUBB_LABEL_OVERRIDE = ["current clubb", "previous clubb"] # To use folder names, comment this line out


SAM_LABEL = "SAM-LES"
COAMPS_LABEL = "COAMPS-LES"
GOLAZ_LABEL = "CLUBB r408 'best ever'"
HOC_LABEL = 'HOC 12/17/2005'

# Text styles
# -------------
# Text/Font
SMALL_FONT_SIZE = 12
MEDIUM_FONT_SIZE = 14
LARGE_FONT_SIZE = 18

DEFAULT_TEXT_SIZE = SMALL_FONT_SIZE
TITLE_TEXT_SIZE = LARGE_FONT_SIZE
AXES_LABEL_FONT_SIZE = LARGE_FONT_SIZE
AXES_TITLE_FONT_SIZE = LARGE_FONT_SIZE
X_TICKMARK_FONT_SIZE = MEDIUM_FONT_SIZE
Y_TICKMARK_FONT_SIZE = MEDIUM_FONT_SIZE
LEGEND_FONT_SIZE = SMALL_FONT_SIZE

# Misc
SHOW_GRID = True
