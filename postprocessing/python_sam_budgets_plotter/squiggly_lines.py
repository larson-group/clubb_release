#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import help.plot_budgets as pb
from help.plot_defs import *
import cases.bomex_case as bc
from matplotlib import pyplot as mpp
from matplotlib.transforms import Affine2D
from matplotlib.backends.backend_pdf import PdfPages
from netCDF4 import Dataset
import numpy as np
import matplotlib.path as mpath
import matplotlib.patches as mpatches
import matplotlib as mpl

#wave_col = 'cornflowerblue'
wave_col = 'purple'
env_col = '#377eb8'

def draw_vertical_wavy_arrow(ax,xcoord,start,end,nwaves,amplitude):
    l = end-start
    #amplitude = l/nwaves/2
    x = np.linspace(start,end,1000)
    y = xcoord+amplitude*np.sin((x-start)/l*2*np.pi*nwaves)
    ax.plot(y,x,color='purple',ls='-', lw=5)
    ax.arrow(end, xcoord, (end-start)*.1,0, width=1,length_includes_head=True)

def curly_arrow(start, end, head_len = 1, head_width = 1, n = 5, col='gray', lw=1., width = 0.1):
    xmin, ymin = start
    xmax, ymax = end
    dist = np.sqrt((xmin - xmax)**2 + (ymin - ymax)**2)
    n0 = dist / (2 * np.pi)
    
    x = np.linspace(0, dist, 151) + xmin
    y = width * np.sin(n*1.005 * (x-xmin) / n0) + ymin
    line = mpp.Line2D(x,y, color=col, lw=lw)
    
    del_x = xmax - xmin
    del_y = ymax - ymin
    ang = np.arctan2(del_y, del_x)
    
    line.set_transform(mpl.transforms.Affine2D().rotate_around(xmin, ymin, ang) + ax.transData)
    ax.add_line(line)
    
    
    verts = np.array([[0,head_width],[0,-head_width],[head_len,0],[0,head_width]]).astype(float)
    verts[:,1] += ymax
    verts[:,0] += xmax
    path = mpath.Path(verts)
    patch = mpatches.PathPatch(path, fc=col, ec=col)
    
    patch.set_transform(mpl.transforms.Affine2D().rotate_around(xmax, ymax, ang) + ax.transData)
    return patch

nc = Dataset(bc.sam_file)
#TODO: incorporate cld for environment profile
u = np.squeeze(nc.variables['U'])
ucld = np.squeeze(nc.variables['UCLD'])
cld = np.squeeze(nc.variables['CLD'])
uenv = (u-cld*ucld)/(1-cld)
#print(np.nanmean(1-cld,axis=0))
uenv = np.nanmean(uenv[181:],axis=0)
print(uenv)
#uw = np.squeeze(nc.variables['UW'])
#uwcld = np.squeeze(nc.variables['UWCLD'])
ucld = np.where(ucld<-1000,np.nan,ucld)
#uwcld = np.where(uwcld<-1000,np.nan,uwcld)
h = np.squeeze(nc.variables['z'])
u = np.nanmean(u[181:,:], axis=0)
#uw = np.nanmean(uw[181:,:], axis=0)
ucld = np.nanmean(ucld[181:,:], axis=0)
#uwcld = np.nanmean(uwcld[181:,:], axis=0)

fontsizes = {
    'labels' : 20,
    'ticks' : 15,
    'title' : 20,
    'legend' : 13.5,
    }

figsize=np.array(mpp.rcParams['figure.figsize'])
scale = .85
figsize=(7.7,6)
fig, ax = mpp.subplots(1,1,figsize=figsize)

fig.suptitle(r'Eastward wind, BOMEX, t=181-360 min', y=.95, fontsize=fontsizes['title'])
#u_line = ax.plot(uenv,h, color='k',ls='-', lw=3, label=r'Layer avg. wind, $\overline{u}$')
u_line = ax.plot(uenv,h, color=env_col,ls='--', lw=3, label=r'Env. avg. wind, $\overline{u}^\mathrm{env}$')
#u_line = ax.plot(u,h, color='k',ls='-', lw=3, label=r'$\overline{u}$')
ucld_line = ax.plot(ucld,h, color='darkorange',ls=':',lw=5, label=r'Cloud avg. wind, $\overline{u}^\mathrm{{cld}}$')
#u_total = ax.plot(u,h, color='r',ls=':', lw=3, label=r'Layer avg. wind, $\overline{u}$')
#ucld_line = ax.plot(ucld,h, color='darkorange',ls=':',lw=5, label=r'$\overline{u}^\mathrm{{cld}}$')
ax.set_xlim(-8.2,-6.25)
ax.set_ylim(0,2000)
xmin,xmax = ax.get_xlim()
ax.set_xlabel(r'$\overline{u}^\mathrm{{cld}},\ \overline{u}^\mathrm{env}\ \mathrm{\left[m\,s^{-1}\right]}$', fontsize=fontsizes['labels'])
ax.set_ylabel('Height [m]', fontsize=fontsizes['labels'],x=10)
jet_max = h[u.argmin()]
idx = np.nanargmin(np.abs(u-ucld))
intersect = h[idx]
print(jet_max)
print(intersect)

# horizontal line
ax.axhline(intersect,color='k',ls='--')
#mpp.text(-7.75,920,"Region of upgradient",color='black')
#mpp.text(-7.75,880,"momentum flux",color='black')
ax.legend(bbox_to_anchor=(.53,.73),fontsize=fontsizes['legend'])
ax.tick_params(axis='both', labelsize=fontsizes['ticks'])

#draw_vertical_wavy_arrow(ax,-8,250,1000,5,.1)
ax.add_patch(curly_arrow((-7.2,550),(-7.2,1380),100,.02,width=.02,lw=2,col=wave_col))
#mpp.annotate('',xy=(-7.2,1480),xytext=(-7.2,600),arrowprops=dict(width=5,headwidth=10,fc='purple',ec='purple'))
mpp.text(-7.17,750,r"$(\overline{u}^\mathrm{env}-c)<0$",color=wave_col,fontsize=14)
ax.add_patch(curly_arrow((-7.4,1750),(-7.4,1420),60,.015,n=3,width=.013,lw=1.8,col=wave_col))
mpp.text(-7.38,1700,r"$(\overline{u}^\mathrm{env}-c)>0$",color=wave_col,fontsize=14)
#mpp.annotate('',xy=(-7.2,1550),xytext=(-7.2,1900),arrowprops=dict(width=5,headwidth=10,fc='purple',ec='purple'))
ax.add_patch(curly_arrow((-7.7,1250),(-7.7,1170),30,.008,n=1.5,width=.007,lw=1.5,col=wave_col))
#mpp.annotate('',xy=(-7.7,1150),xytext=(-7.7,1280),arrowprops=dict(width=4,headwidth=8,fc='purple',ec='purple'))
#mpp.text(-6.45,100,r"$w'>0$",fontsize=20)
#mpp.text(-5,1800,r"a)",fontsize=20)
#mpp.annotate('',xy=(-6.5,250),xytext=(-7.15,250),arrowprops=dict(arrowstyle='->',lw=2))
#mpp.text(-7.05,300,r"$u'>0$",fontsize=20)
##axes[i].xaxis.get_offset_text().set_size(fontsizes['ticks'])
fig.savefig('/home/sdomke/workspace/plotgen_out/publishing_runs/squiggly.png')
pdf = PdfPages('/home/sdomke/workspace/plotgen_out/publishing_runs/squoggly.pdf')
pdf.savefig(fig)
pdf.close()