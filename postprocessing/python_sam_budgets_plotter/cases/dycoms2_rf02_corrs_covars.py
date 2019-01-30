#-------------------------------------------------------------------------------
# G E N E R E L   D E F A U L T   S E T T I N G S
#-------------------------------------------------------------------------------
name = 'DYCOMS2_RF02_Corrs_Covars'
enabled = True
type = 'standard'
startTime = 61.0      # [minutes]
endTime = 180.0      # [minutes]
startHeight = 0.0    # [m]
endHeight = 16000.0   # [m]
headerText = 'DYCOMS II RF02 Correlations and Covariances Minutes ' + str(startTime) + '-' + str(endTime) + ', ' + str(startHeight) + 'm-' + str(endHeight) + 'm'

case = 'DYCOMS2_RF02_corrs_covars'
out_dir = './clubb/plotsSam/%s/'%(case)
sam_file = '/home/ckrome/clubb/outputSam/DYCOMS_RF02_128x128x96_dr_nosed.nc'
plot_case_name = '%s_'%(case)

#-------------------------------------------------------------------------------
# G E N E R E L   P L O T   S E T T I N G S
#-------------------------------------------------------------------------------
lineWidth = 2
color = 'nipy_spectral'
yLabel = 'Height [m]'