"""
:author: Nicolas Strike
:date: Late 2019

Style_definitions.py contains the style formatting parameters used by Pyplotgen throughout the source code.
If you wish to change a line color, font size, etc., then this is the place to do it.
"""

# Line Styles
DEFAULT_LINE_STYLE = "r-"
LES_LINE_STYLE = "k-"
GOLAZ_BEST_R408_LINE_STYLE = "g-"

DEFAULT_LINE_THICKNESS = 2
LES_LINE_THICKNESS = 4
THIN_LINE_THICKNESS = 1 # Used if --thin-lines is passed, overrides every other size


# Line Labels
DEFAULT_LABEL = "current clubb"
SAM_LABEL = "SAM-LES"
COAMPS_LABEL = "COAMPS-LES"
GOLAZ_LABEL = "CLUBB r408 'best ever'"

# Text/Font
SMALL_FONT_SIZE = 10
MEDIUM_SIZE = 12
LARGE_FONT_SIZE = 16

DEFAULT_TEXT_SIZE = SMALL_FONT_SIZE
TITLE_TEXT_SIZE = LARGE_FONT_SIZE
AXES_LABEL_FONT_SIZE = LARGE_FONT_SIZE
AXES_TITLE_FONT_SIZE = LARGE_FONT_SIZE
X_TICKMARK_FONT_SIZE = MEDIUM_SIZE
Y_TICKMARK_FONT_SIZE = MEDIUM_SIZE
LEGEND_FONT_SIZE = SMALL_FONT_SIZE

# Misc
SHOW_GRID = True