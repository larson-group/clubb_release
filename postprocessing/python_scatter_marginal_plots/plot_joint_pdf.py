# $Id$
#
# plot_joint_pdf.py
# 
# Description: Includes a definition that will read in 1D data from two
#               different models and plot the scatter-plot and joint pdfs
#               of the two variates. This is callable from outside the file.
#               Example code for extracting data is provided, however. 
# An example call to the definition:

# plot_joint_pdf(sam_chi,clb_chi,'chi','[kg kg^-1]',sam_rr,clb_rr,'rrm','[kg s^-1]','Default','Double gamma_coef',Title)

import netCDF4
import numpy as np
import matplotlib.pyplot as plt

# User input
#-------------------------------------------------------------------------------
#CLUBB's samples
clb_silhs_filename = '../../output/lba_nl_lh_sample_points_2D.nc'
#SAM's standard output
sam_filename = '../../output/3D_Output/LBA_128kmx1kmx128_1km_Morrison_64_000000????.nc'
#Output from Morrison Micro.
sam_micro_filename = '../../output/3D_Output/LBA_128kmx1kmx128_1km_Morrison_64_000000????_micro.nc'

analysis_time = 200 #    [min]
analysis_height = 1000 # [m]

Title = "First iteration of a script, under construction"
x_unit = 'kg kg^-1'
y_unit = 'm s^-1'

model1 = 'SAM'
SAM_VARx = 'CHI'
SAM_VARy = 'W'

model2 = 'CLUBB'
CLUBB_VARx = 'chi'
CLUBB_VARy = 'w'

def plot_joint_pdf(x1,x2,x_name,x_unit,
                   y1,y2,y_name,y_unit,
                    model1,model2,Title):
    # Description:
    #   This definition takes in 1D data and plots the joint pdf, 
    #   including the marginals
    
    import numpy as np
    import pylab as pl
    import matplotlib.pyplot as plt

    fig = plt.figure()
    fig.text(.5,.95,Title,horizontalalignment='center',)
    print "Created figure"
    ax1 = plt.subplot2grid((3, 3), (1, 0),rowspan=2,colspan=2)
    ax2 = plt.subplot2grid((3, 3), (0, 0),colspan=2,sharex=ax1)
    ax3 = plt.subplot2grid((3, 3), (1, 2),rowspan=2,sharey=ax1)
    print "Created subplots"
    ax1.scatter(x1,y1,alpha=.3,color='r',label=model1)
    ax1.scatter(x2,y2,alpha=.3,color='b',label=model2)
    ax1.set_xlabel(x_name + x_unit)
    ax1.set_ylabel(y_name + y_unit)
    ax1.legend(loc=2)
    print "Plot1"
    ax2.hist(x1,bins=100,normed=True,alpha=.5,color='r')
    ax2.hist(x2,bins=100,normed=True,alpha=.5,color='b')
    pl.setp( ax2.get_xticklabels(), visible=False)
    print "Plot2"
    ax3.hist(y1,bins=100,normed=True,orientation='horizontal',alpha=.5,color='r')
    ax3.hist(y2,bins=100,normed=True,orientation='horizontal',alpha=.5,color='b')
    pl.setp( ax3.get_yticklabels(), visible=False)
    pl.setp( ax3.get_xticklabels(), rotation=-90)
    print "Plot3"
    fig.show()
    
# Fetch data. Should not have to edit below here
#-------------------------------------------------------------------------------

# Time and height and relevant indicies for CLUBB and SAM
nl_nc = netCDF4.Dataset(clb_silhs_filename)
clb_z = nl_nc.variables['altitude']
clb_time = nl_nc.variables['time']
indx_cz = (np.abs(clb_z[:] - analysis_height)).argmin()
indx_ct = (np.abs(clb_time[:]/60. - analysis_time)).argmin()

nl_nc = netCDF4.MFDataset(sam_filename)
sam_z = nl_nc.variables['z']
sam_time = nl_nc.variables['time']
indx_sz = (np.abs(clb_z[:] - analysis_height)).argmin()
indx_st = (np.abs(sam_time[:] - analysis_time)).argmin()

# Now extract SAM variables
nl_nc = netCDF4.MFDataset(sam_micro_filename)
x1 = nl_nc.variables[SAM_VARx]
x1 = np.swapaxes(x1,1,3)
x1 = np.squeeze(x1)

nl_nc = netCDF4.MFDataset(sam_filename)
y1 = nl_nc.variables[SAM_VARy]
y1 = np.swapaxes(y1,1,3)
y1 = np.squeeze(y1)

# Find grid dimensions so we can reshape SAM's output into one column
sam_gridx = np.size(x1,axis=1)
sam_gridy = np.size(x1,axis=2)
x1 = x1[indx_st,:,:,indx_sz]
y1 = y1[indx_st,:,:,indx_sz]

# Now extract clubb variables of interest
nl_nc = netCDF4.Dataset(clb_silhs_filename)
x2 = nl_nc.variables[CLUBB_VARx]
y2 = nl_nc.variables[CLUBB_VARy]
x2 = x2[indx_ct,indx_cz,:,0]
y2 = y2[indx_ct,indx_cz,:,0]


x1 = np.squeeze(np.reshape(x1,(sam_gridx*sam_gridy,1)))
y1 = np.squeeze(np.reshape(y1,(sam_gridx*sam_gridy,1)))
x2 = np.squeeze(x2)
y2 = np.squeeze(y2)

plot_joint_pdf(x1,x2,SAM_VARx,x_unit, y1,y2,SAM_VARy,y_unit, model1,model2,Title)
