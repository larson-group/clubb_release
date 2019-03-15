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
## PROGRAM DEFINITIONS
ntrials = 3 # maximum number of input trials before closing the program
affirmatives = ['y', 'yes', 'aye', 'yay', 'pos', '1']
negatives = ['n', 'no', 'nay', 'nope', 'neg', '0']

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
    "BOMEX 64x64" : "bomex_case",
    "BOMEX 128x128" : "bomex_large_case",
    "DYCOMS_RF01" : "dycoms2_rf01_case",
    "DYCOMS_RF02" : "dycoms2_rf02_case",
    "RICO" : "rico_case",
    "LBA" : "lba_case",
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
    "BOMEX" : 500,
    "DYCOMS_RF01" : 750,
    "DYCOMS_RF02" : 750,
    "RICO" : 1500,      # grid is huge, individual arrows can barely be distinguished
    "LBA" : 500,
    } # list of cases, append as needed

figure_scale = {
    "BOMEX" : 1,
    "DYCOMS_RF01" : 1,
    "DYCOMS_RF02" : 1,
    "RICO" : 2,      # grid is huge, individual arrows can barely be distinguished
    "LBA" : 1,
    } # list of cases, append as needed

# Power limit for tick floating point formatting
pow_lim = 1

# List of fontsizes
fontsizes = {
    'labels' : 25,
    'ticks' : 20,
    'title' : 30,
    'legend' : 18,
    }

# Line formatting for comparison plots
comp_style = {
    'clubb' : {
        'color' : 'red',
        'lw'    : 3,
        'ls'    : '--',
        'label' : 'new CLUBB'
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
        'label' : 'old CLUBB'
        },
    }

# List of line styles (cycling)
styles = ['-','--','-.',':']

# Legend positions for both x-axes
legend_pos = [2,1]

# Distinguishing titles for legends
legend_title = ['bottom', 'top']

# List of distinguishable colors for multiple plots lines
#             red       blue      green     purple    brown     black     grey      orange    magenta
color_arr = ['#e41a1c','#377eb8','#4daf4a','#984ea3','#a65628','#000000','#999999','#ffa600','#f750bf']

## Horizontal cloud slice plots
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

# TODO: Define dictionary for fontsizes in 3d plot
