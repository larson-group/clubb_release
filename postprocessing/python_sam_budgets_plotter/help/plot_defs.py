"""
-------------------------------------------------------------------------------
    DESCRIPTION
-------------------------------------------------------------------------------
Define globally required variables, lists, and functions
"""
#-------------------------------------------------------------------------------
#    I M P O R T S
#-------------------------------------------------------------------------------
import os
import logging

#-------------------------------------------------------------------------------
#   L O G G E R
#-------------------------------------------------------------------------------
logger = logging.getLogger('plotgen.help.defs')
#logger.setLevel(logging.INFO)
logger.setLevel(logging.DEBUG)
#logger.setLevel(logging.CRITICAL)

#-------------------------------------------------------------------------------
#    F U N C T I O N S
#-------------------------------------------------------------------------------
def isFunction(value):
    logger.info('__isFunction__')
    isFunc = False
    if '+' in value:
        isFunc = True
    elif '-' in value:
        isFunc = True
    elif '*' in value:
        isFunc = True
    elif '/' in value:
        isFunc = True
    else:
        isFunc = False
    return isFunc

def makeDirectory(pathToFile):
    if not os.path.exists(pathToFile):
        os.makedirs(pathToFile)

#-------------------------------------------------------------------------------
#    GLOBAL DEFINITIONS
#-------------------------------------------------------------------------------
## CONSTANTS
cld_lim = 0 # lower bound of water content at which grid point is in cloud
## PROGRAM DEFINITIONS
ntrials = 3 # maximum number of input trials before closing the program
affirmatives = ['y', 'yes', 'aye', 'yay', 'pos', '1']
negatives = ['n', 'no', 'nay', 'nope', 'neg', '0']

## CASE FILE DEFINITIONS
# define constants to use as indexer instead of integers
PLOT_ID = 0
PLOT_TITLE = 1
PLOT_XLABEL = 2
PLOT_TEXT = 3
PLOT_TEXTPOS = 4
PLOT_LINES = 5
PLOT_LINES_SAM = 5
PLOT_LINES_CLUBB = 6

# define line list indexer constants
LINE_NAME = 0
LINE_VISIB = 1
LINE_EXPRESSION = 2
LINE_FACTOR = 3
LINE_XAX = 4

## PLOTTING CHOICES
# List of output types
# TODO: How to define??
type_dict = {
    "SAM budgets" : "sam_budget",
    "SAM standalone" : "sam_standalone",
    "SAM correlations and covariances (not implemented)" : "sam_corr_covar",
    "SAM 3D plots" : "sam_3d",
    "CLUBB budgets" : "clubb_budget",
    "CLUBB standalone" : "clubb_standalone",
    "SAM CLUBB comparison" : "sam_clubb_comparison",
    } # list of types, append as needed
type_names = type_dict.keys()
type_modules = type_dict.values()

# Template naming convention for help files
type_name_template = '{}_variables'

# List of cases
case_dict = {
    "ARM9707"       : "arm9707_case",
    "BOMEX 64x64"   : "bomex_case",
    "BOMEX 128x128" : "bomex_large_case",
    "DYCOMS_RF01"   : "dycoms2_rf01_case",
    "DYCOMS_RF02"   : "dycoms2_rf02_case",
    "RICO"          : "rico_case",
    "LBA"           : "lba_case",
    "GCSSARM"       : "gcssarm_case",
    "IOP"           : "iop_case",
    "GABLS3_night"  : "gabls3_night_case",
    } # list of cases, append as needed
case_names = case_dict.keys()
case_modules = case_dict.values()

## FORMATTING
# Define date format to be used in file names by datetime.datetime.strftime
date_file_format = '%Y%m%d'

## FILE IMPORT PATTERNS
# Define search patterns
dt_pattern = "dt = ([0-9]*[.[0-9]*]?)"
prm_patterns = {
    'dt' : "dt\s*=\s*([0-9]*[.[0-9]*]?)",
    'dx' : "dx\s*=\s*([0-9]*[.[0-9]*]?)",
    'dy' : "dy\s*=\s*([0-9]*[.[0-9]*]?)",
    'nstop' : "nstop\s*=\s*([0-9]*[.[0-9]*]?)",
    'nsave3Dstart' : "nsave3Dstart\s*=\s*([0-9]*[.[0-9]*]?)",
    'nsave3Dend' : "nsave3Dend\s*=\s*([0-9]*[.[0-9]*]?)",
    }

## PLOT FORMATTING
quiver_scale_factor = {
    "BOMEX"         : 500,
    "DYCOMS_RF01"   : 750,
    "DYCOMS_RF02"   : 750,
    "RICO"          : 1500,      # grid is huge, individual arrows can barely be distinguished
    "LBA"           : 500,
    "ARM9707"       : 500,
    } # list of cases, append as needed

figure_scale = {
    "BOMEX"         : 1,
    "DYCOMS_RF01"   : 1,
    "DYCOMS_RF02"   : 1,
    "RICO"          : 2,      # grid is huge, individual arrows can barely be distinguished
    "LBA"           : 1,
    "ARM9707"       : 1,
    } # list of cases, append as needed

# Power limit for tick floating point formatting
pow_lims = (-2,3)

# List of fontsizes
fontsizes = {
    'labels' : 25,
    'ticks' : 22,
    'title' : 30,
    'legend' : 24,
    }

# Line formatting for comparison plots
comp_style = {
    'clubb' : {
        'color' : 'red',
        'lw'    : 3,
        'ls'    : '--',
        'label' : 'prog. mom. flux'
        },
    'sam'   : {
        'color' : 'black',
        'lw'    : 5,
        'ls'    : '-',
        'label' : 'SAM-LES'
        },
    'old'   : {
        'color' : 'green',
        'lw'    : 3,
        'ls'    : ':',
        'label' : 'downgradient dfsn.'
        },
    }

# List of line styles (cycling) (BUG: ls '-.' is not displayed correctly in legend!)
styles = ['-','--',':']
#styles = [':', '-','--']

# Legend positions for both x-axes
#legend_pos = [2,1]
#legend_pos = [1,2]
legend_pos = [0,0]

# Distinguishing titles for legends
legend_title = ['bottom', 'top']

# List of distinguishable colors for multiple plot lines
#             red       blue      green     purple    brown     black     grey      orange    magenta
color_arr = ['#e41a1c','#377eb8','#4daf4a','#984ea3','#a65628','#000000','#999999','#ffa600','#f750bf']

# Alternative styles for conditional averages
#               orange      blue      black   green   purple      brown      grey    magenta
#color_arr = ['darkorange', '#377eb8', 'k', '#4daf4a','#984ea3','#a65628','#999999','#f750bf']
#styles = [':', '--','-']
# Conditional style 2
##             red       green     blue      bla             black     brown     grey      orange    magenta
#color_arr = ['#e41a1c','#4daf4a','#377eb8','lightseagreen','#000000','#a65628','#999999','#ffa600','#f750bf']
#styles = ['--','--','--',':','-']

## Heinze line colors (solid, dashed): 
#   - lime: pressure redistribution (Pi (resolved), P (subgrid scale))
#   - dodgerblue/deepskyblue: turbulent transport (T^t), pressure transport (T^p)
#   - red: buoyancy (B)
#   - black: mean-velocity shear (G), viscous dissipation (D)
# Budget color array
# advection, buoyancy, dissipation, isotropy, pressure, turb. prod. time tndncy, residual
#color_arr = ['dodgerblue', 'red', 'black', 'lime', 'dodgerblue', 'black', 'grey', 'orange']
#styles = ['-', '-', '--', '-', '--', '-', '-', ':']

## Horizontal cloud slice plots
# Define base figsize
#figsize = np.array(mpp.rcParams['figure.figsize'])*4
profile_figsize=(11,10)
figsize = (32.5,20)

# Line colors for secondary plots in cloud plot
col_3d = {
    'uw' : 'g',
    'vw' : 'r',
    }

# Width ratios of primary to secondary plot
width_ratio = [5,1]

# Interpolation algorithm for imshow of vertical wind component
cloud_interpolation = 'quadric'

# Legend location in seconary subplot
legend_pos_3d = 1

# Set dilation length for halo definition (1 is not enough!)
dil_len = 4

# Set quantile for outlier cutoff: 0<=quant<=100, enter 0 for min/max
# The amount of outliers seems to be approximately constant with array size, but in order to get a visually acceptable distribution, we need to take quantiles
quant = 1e-2
#trash = 20

# Set point skip for horizontal plots
skip = 4

# Hatching linewidth for cloud contours
hw = 2

# Cloud contour and hatching color
cld_col = 'blue'
hatch_col = (0,0,1,.3)
#cld_col = 'skyblue' # This one was used for coloring patches, changed to two different lines

# Cloud halo contour color
#halo_col = 'mediumblue' # This is too similar to the cloud color, try grey/silver
halo_col = 'grey'

# Segment colors for background image
red = (.9,0,0)
green = (0,.75,0)
zero_col = 1

# TODO: Define dictionary for fontsizes in 3d plot
