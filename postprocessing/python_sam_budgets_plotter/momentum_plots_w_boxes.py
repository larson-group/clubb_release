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


rect_col = ['grey', 'lightcoral', 'cornflowerblue']
rect_lab = ['Downgradient flux', 'Upgradient flux', 'Weak flux']

nc = Dataset(bc.sam_file)

u = np.squeeze(nc.variables['U'])
ucld = np.squeeze(nc.variables['UCLD'])
uw = np.squeeze(nc.variables['UW'])
uwcld = np.squeeze(nc.variables['UWCLD'])
ucld = np.where(ucld<-1000,np.nan,ucld)
uwcld = np.where(uwcld<-1000,np.nan,uwcld)
h = np.squeeze(nc.variables['z'])
u = np.nanmean(u[181:,:], axis=0)
uw = np.nanmean(uw[181:,:], axis=0)
ucld = np.nanmean(ucld[181:,:], axis=0)
uwcld = np.nanmean(uwcld[181:,:], axis=0)


fontsizes = {
    'labels' : 20,
    'ticks' : 15,
    'title' : 20,
    'legend' : 16,
    }

figsize=np.array(mpp.rcParams['figure.figsize'])
scale = .85
figsize=(7.7,6)
fig, ax = mpp.subplots(1,1,figsize=figsize)

fig.suptitle(r'Eastward wind, BOMEX, t=181-360 min', y=.95, fontsize=fontsizes['title'])
u_line = ax.plot(u,h, color='k',ls='-', lw=3, label=r'Layer-averaged wind, $\overline{u}$')
#u_line = ax.plot(u,h, color='k',ls='-', lw=3, label=r'$\overline{u}$')
ucld_line = ax.plot(ucld,h, color='darkorange',ls=':',lw=5, label=r'Cloud-averaged wind, $\overline{u}^\mathrm{{cld}}$')
#ucld_line = ax.plot(ucld,h, color='darkorange',ls=':',lw=5, label=r'$\overline{u}^\mathrm{{cld}}$')
#ax.set_xlim(-8.2,-7.5)
ax.set_ylim(0,2000)
xmin,xmax = ax.get_xlim()
ax.set_xlabel(r'$\overline{u}^\mathrm{{cld}},\ \overline{u}\ \mathrm{\left[m\,s^{-1}\right]}$', fontsize=fontsizes['labels'])
ax.set_ylabel('Height [m]', fontsize=fontsizes['labels'],x=10)
jet_max = h[u.argmin()]
tmp = np.roll(uwcld,1)*uwcld
idx = np.where(tmp<0)[0][0]
intersect = (h[idx] + h[idx-1])/2.
print(jet_max)
print(intersect)


# Plot rectangles
rect1 = Rectangle((xmin,0),xmax-xmin,jet_max,color='grey',alpha=.5)
rect2 = Rectangle((xmin,jet_max),xmax-xmin,intersect-jet_max,color='lightcoral',alpha=.5)
rect3 = Rectangle((xmin,intersect),xmax-xmin,2000-intersect,color='cornflowerblue',alpha=.5)
ax.add_patch(rect1)
ax.add_patch(rect2)
ax.add_patch(rect3)
#mpp.text(-7.75,920,"Region of upgradient",color='black')
#mpp.text(-7.75,880,"momentum flux",color='black')
lines = u_line+ucld_line
rectangles = [Patch(color=rect_col[i],alpha=.5,label=rect_lab[i]) for i in range(len(rect_col))][::-1]
ax.legend(handles = lines+rectangles, fontsize=fontsizes['legend'], loc=7)
ax.tick_params(axis='both', labelsize=fontsizes['ticks'])

mpp.annotate('',xy=(-6.5,250),xytext=(-6.5,40),arrowprops=dict(arrowstyle='->',lw=2))
mpp.text(-6.45,100,r"$w'>0$",fontsize=20)
mpp.text(-5,1800,r"a)",fontsize=20)
mpp.annotate('',xy=(-6.5,250),xytext=(-7.15,250),arrowprops=dict(arrowstyle='->',lw=2))
mpp.text(-7.05,300,r"$u'>0$",fontsize=20)
#axes[i].xaxis.get_offset_text().set_size(fontsizes['ticks'])
fig.savefig('/home/sdomke/workspace/plotgen_out/publishing_runs/ucld_w_regions.png')
pdf = PdfPages('/home/sdomke/workspace/plotgen_out/publishing_runs/u_w_boxes.pdf')
pdf.savefig(fig)
pdf.close()

ax.clear()
fig.suptitle(r"Momentum flux $\overline{u'w'}$, BOMEX, t=181-360 min", y=.96, fontsize=fontsizes['title'])
ax.axvline(0,c='k',ls='--')
uw_line = ax.plot(uw,h, color='k',ls='-', lw=3, label=r"Layer-averaged flux, $\overline{u'w'}$")
#uw_line = ax.plot(uw,h, color='k',ls='-', lw=3, label=r"$\overline{u'w'}$")
#uwcld_line = ax.plot(uwcld,h, color='k',ls='-', lw=3, label=r"Cloud averaged, $\overline{u'w'}^\mathrm{cld}$")
ax.set_ylim(0,2000)
xmin,xmax = ax.get_xlim()
ax.set_xlabel(r"Momentum flux, $\overline{u'w'}\ \mathrm{\left[m^{2}\,s^{-2}\right]}$", fontsize=fontsizes['labels'])
ax.set_ylabel('Height [m]', fontsize=fontsizes['labels'])

# Plot rectangles
rect1 = Rectangle((xmin,0),xmax-xmin,jet_max,color='grey',alpha=.5)
rect2 = Rectangle((xmin,jet_max),xmax-xmin,intersect-jet_max,color='lightcoral',alpha=.5)
rect3 = Rectangle((xmin,intersect),xmax-xmin,2000-intersect,color='cornflowerblue',alpha=.5)
ax.add_patch(rect1)
ax.add_patch(rect2)
ax.add_patch(rect3)
#mpp.text(-7.75,920,"Region of upgradient",color='black')
#mpp.text(-7.75,880,"momentum flux",color='black')
#ax.legend(handles = uw_line+uwcld_line+rectangles)
ax.legend(handles = uw_line+rectangles, fontsize=fontsizes['legend'], loc=7)
ax.tick_params(axis='both', labelsize=fontsizes['ticks'])
mpp.text(.07,1800,r"b)",fontsize=20)

fig.savefig('/home/sdomke/workspace/plotgen_out/publishing_runs/uw_w_regions.png')
pdf = PdfPages('/home/sdomke/workspace/plotgen_out/publishing_runs/uw_w_boxes.pdf')
pdf.savefig(fig)
pdf.close()