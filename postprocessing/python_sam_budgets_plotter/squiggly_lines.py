#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import help.plot_budgets as pb
from help.plot_defs import *
import cases.bomex_case as bc
from matplotlib import pyplot as mpp
from matplotlib.backends.backend_pdf import PdfPages
from netCDF4 import Dataset
import numpy as np
from matplotlib.patches import Rectangle, Patch


def draw_wavy_arrow(ax,start,end,amplitude,):
    

nc = Dataset(bc.sam_file)

u = np.squeeze(nc.variables['U'])
ucld = np.squeeze(nc.variables['UCLD'])
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
u_line = ax.plot(u,h, color='k',ls='-', lw=3, label=r'Layer avg. wind, $\overline{u}$')
#u_line = ax.plot(u,h, color='k',ls='-', lw=3, label=r'$\overline{u}$')
ucld_line = ax.plot(ucld,h, color='darkorange',ls=':',lw=5, label=r'Cloud avg. wind, $\overline{u}^\mathrm{{cld}}$')
#ucld_line = ax.plot(ucld,h, color='darkorange',ls=':',lw=5, label=r'$\overline{u}^\mathrm{{cld}}$')
ax.set_xlim(-8.2,-6.25)
ax.set_ylim(0,2000)
xmin,xmax = ax.get_xlim()
ax.set_xlabel(r'$\overline{u}^\mathrm{{cld}},\ \overline{u}\ \mathrm{\left[m\,s^{-1}\right]}$', fontsize=fontsizes['labels'])
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

mpp.annotate('',xy=(-7.2,1480),xytext=(-7.2,600),arrowprops=dict(width=5,headwidth=10,fc='purple',ec='purple'))
mpp.text(-7.15,750,r"$(\overline{u}-c)<0$",color='purple',fontsize=15)
mpp.annotate('',xy=(-7.2,1550),xytext=(-7.2,1900),arrowprops=dict(width=5,headwidth=10,fc='purple',ec='purple'))
mpp.annotate('',xy=(-7.7,1150),xytext=(-7.7,1280),arrowprops=dict(width=4,headwidth=8,fc='purple',ec='purple'))
#mpp.text(-6.45,100,r"$w'>0$",fontsize=20)
#mpp.text(-5,1800,r"a)",fontsize=20)
#mpp.annotate('',xy=(-6.5,250),xytext=(-7.15,250),arrowprops=dict(arrowstyle='->',lw=2))
#mpp.text(-7.05,300,r"$u'>0$",fontsize=20)
##axes[i].xaxis.get_offset_text().set_size(fontsizes['ticks'])
fig.savefig('/home/sdomke/workspace/plotgen_out/publishing_runs/squiggly.png')
pdf = PdfPages('/home/sdomke/workspace/plotgen_out/publishing_runs/squoggly.pdf')
pdf.savefig(fig)
pdf.close()