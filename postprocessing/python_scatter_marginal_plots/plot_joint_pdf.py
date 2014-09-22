#
# plot_joint_pdf.py
# 
# Description: Includes a definition that will read in 1D data from two
#               different models and plot the scatter-plot and joint pdfs
#               of the two variates. This is callable from outside the file.
#               Example code for extracting data is provided, however. 
# An example call to the definition:

# plot_joint_pdf(sam_chi,clb_chi,'chi','[kg kg^-1]',sam_rr,clb_rr,'rrm','[kg s^-1]','Default','Double gamma_coef',Title)


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