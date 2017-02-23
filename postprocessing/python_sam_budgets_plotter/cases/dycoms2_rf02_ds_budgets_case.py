#-------------------------------------------------------------------------------
# G E N E R E L   D E F A U L T   S E T T I N G S
#-------------------------------------------------------------------------------
name = 'DYCOMS2_RF02_BUDGETS'
enabled = True
type = 'budget'
startTime = 0.0      # [minutes]
endTime = 30.0      # [minutes]
startHeight = 0.0    # [m]
endHeight = 1200.0   # [m]
headerText = 'DYCOMS II RF02 BUDGETS Budgets Minutes ' + str(startTime) + '-' + str(endTime)

case = 'DYCOMS2_RF02_BT_30'
out_dir = './clubb/plotsSam/%s/'%(case)
sam_file = '/home/ckrome/clubb/outputSam/DYCOMS_RF02_128x128x96_dr_nosed.nc'
plot_case_name = '%s_'%(case)

#-------------------------------------------------------------------------------
# G E N E R E L   P L O T   S E T T I N G S
#-------------------------------------------------------------------------------
lineWidth = 2
color = 'nipy_spectral'
yLabel = 'Height [m]'