#-------------------------------------------------------------------------------
# G E N E R E L   D E F A U L T   S E T T I N G S
#-------------------------------------------------------------------------------
name = 'LBA_Corrs_Covars'
enabled = True
type = 'standard'
startTime = 120.0      # [minutes]
endTime = 360.0      # [minutes]
startHeight = 0.0    # [m]
endHeight = 27000.0   # [m]
headerText = 'LBA Correlations and Covariances Minutes ' + str(startTime) + '-' + str(endTime) + ', ' + str(startHeight) + 'm-' + str(endHeight) + 'm'

case = 'LBA_corrs_covars'
out_dir = './clubb/plotsSam/%s/'%(case)
sam_file = '/home/ckrome/clubb/outputSam/LBA_128kmx128kmx128_1km_Morrison.nc'
plot_case_name = '%s_'%(case)

#-------------------------------------------------------------------------------
# G E N E R E L   P L O T   S E T T I N G S
#-------------------------------------------------------------------------------
lineWidth = 2
color = 'nipy_spectral'
yLabel = 'Height [m]'