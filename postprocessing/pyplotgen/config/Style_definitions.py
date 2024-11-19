"""
:author: Nicolas Strike
:date: Late 2019

Style_definitions.py contains the style formatting parameters used by Pyplotgen throughout the source code.
If you wish to change a line color, font size, etc., then this is the place to do it.
"""

# For information on making plots "publish ready" see the `Making Publish Ready Plots` portion of the README.md
# https://github.com/larson-group/clubb/blob/master/postprocessing/pyplotgen/README.md#making-publish-ready-plots

# Line Styles
# --------------

# Color rotation
# Multiple color options are given, uncomment the one that's most useful and comment out the rest.
# Feel free to make your own as well.
# Pyplot requires that style rotation's and color rotation's have the same length, so be sure to properly
# (un)comment the style rotations according to the desired color scheme.
# Color/style rotations are used on budget plots and for clubb-folder lines
# If paired_plots is True when plotting a panel, the STYLE_ROTATION will be ignored,
# as it is written to oscillate between dashed and solid lines.
# This currently happens by default for non-budget plots.

# Colors are taken from the colorbrewer qualitative 10-color set:
# http://colorbrewer2.org/#type=qualitative&scheme=Paired&n=10
# Max amount of colors (includes pastel)
# Not colorblind or print safe
# COLOR_ROTATION = ['#1f78b4', '#e31a1c', '#ff7f00', '#6a3d9a', '#33a02c',
#                   '#a6cee3', '#b2df8a', '#fb9a99', '#fdbf6f', '#cab2d6']
# STYLE_ROTATION = ['-', '--', ':', '-.','-', '--', ':', '-.','-', '--']

# Color-blind and print safe color rotation. These requirements only allow up to 4 colors
# COLOR_ROTATION = ['#1f78b4', '#33a02c', '#a6cee3', '#b2df8a']
# STYLE_ROTATION = ['-', '--', ':', '-.']

# Print safe colors (not color-blind friendly). From colorbrewer (removed yellow)
COLOR_ROTATION = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#a65628', '#f781bf', '#999999']
STYLE_ROTATION = ['-', '--', ':', '-.', '-', '--', ':', '-.']

# Details on line styles here: https://matplotlib.org/2.1.2/api/_as_gen/matplotlib.pyplot.plot.html
BENCHMARK_LINE_STYLES = {
    'sam': "k-",
    'coamps': "grey", # dark grey = #404040
    'r408': "g-",
    'hoc': "m-.",
    'e3sm': "c-",
    'wrf': "k-"
}

# Line Thicknesses
CLUBB_LINE_THICKNESS = 2
LES_LINE_THICKNESS = 8
ARCHIVED_CLUBB_LINE_THICKNESS = 4
E3SM_LINE_THICKNESS = 4
THIN_LINE_THICKNESS = 1  # Used if --thin-lines is passed, overrides every other size
# These are used when pyplotgen is oscillating between flat/dashed lines
FLAT_LINE_THICKNESS = 3
DASHED_LINE_THICKNESS = 5

# Line Labels
BENCHMARK_LABELS = {
    'sam': "SAM-LES",
    'coamps': "COAMPS-LES",
    'r408': "CLUBB r408 'best ever'",
    'hoc': 'HOC 12/17/2005',
    'e3sm': "E3SM",
    'wrf': "WRF-LASSO"
}

# Text styles
# -------------
# Text/Font
SMALL_FONT_SIZE = 14
MEDIUM_FONT_SIZE = 16
LARGE_FONT_SIZE = 20

DEFAULT_TEXT_SIZE = SMALL_FONT_SIZE
TITLE_TEXT_SIZE = LARGE_FONT_SIZE
AXES_LABEL_FONT_SIZE = LARGE_FONT_SIZE
AXES_TITLE_FONT_SIZE = LARGE_FONT_SIZE
X_TICKMARK_FONT_SIZE = MEDIUM_FONT_SIZE
Y_TICKMARK_FONT_SIZE = MEDIUM_FONT_SIZE
LEGEND_FONT_SIZE = SMALL_FONT_SIZE

# Misc
SHOW_GRID = True
POW_LIMS = (-2,3)
BUDGET_XAXIS_SCALE_FACTOR = 1.1 # THIS VALUE MUST BE >= 1
MOVIE_XAXIS_SCALE_FACTOR = 0.05 # THIS VALUE SHOULD BE << 1
FRAMES_PER_SECOND = 10 # FOR ANIMATIONS

# Image output configuration
HQ_DPI = 100 # This value overwrites the default IMG_OUTPUT_DPI when --high-quality is passed in via command line
IMG_OUTPUT_DPI = 45 # Recommended default value: 45
FIGSIZE = (10,6) # The x,y proportions (in inches) to pass into matplotlib. It's not recommended to change this.

# Time-height plot settings
# An overview of all colormaps available in matplotlib can be found here:
# https://matplotlib.org/1.4.2/examples/color/colormaps_reference.html
CONTOUR_CMAP_CORR = 'coolwarm' # Default cmap used in matplotlib
CONTOUR_CMAP_NORMED = 'Oranges' # Default cmap used in matplotlib
CONTOUR_CMAP_GENERAL = 'Oranges' # Default cmap used in matplotlib

# Animation settings
ANIM_OUTPUT_FORMAT = 'gif' # or 'mp4'