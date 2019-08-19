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
from matplotlib.patches import Rectangle

nc = Dataset(bc.sam_file)

pdf = PdfPages('/home/sdomke/workspace/thesis/masters/data/upgradient.pdf')
u = np.squeeze(nc.variables['U'])
ucld = np.squeeze(nc.variables['UCLD'])
ucld = np.where(ucld<-1000,np.nan,ucld)
h = np.squeeze(nc.variables['z'])
u = np.nanmean(u[180:,:], axis=0)
ucld = np.nanmean(ucld[180:,:], axis=0)

mpp.suptitle('Upgradient momentum flux within clouds (BOMEX)', y=.93)
mpp.plot(ucld,h, color='#e41a1c',ls='--',lw=3, label=r'Cloud averaged, $\bar{u}^\mathrm{cld}$')
mpp.plot(u,h, color='k',ls='-', lw=3, label=r'Layer averaged, $\bar{u}$')
#mpp.plot(u,h, color='#377eb8',ls='-', lw=3, label=r'Layer averaged, $\bar{u}$')
mpp.xlim(-8.2,-7.5)
mpp.ylim(500,1250)
mpp.legend()
mpp.xlabel(r'$\mathrm{\bar{u}\ \left[\frac{m}{s}\right]}$')
mpp.ylabel('Height [m]')
# w' arrow DONE
mpp.annotate('',xy=(-8.1,600),xytext=(-8.1,500),arrowprops=dict(arrowstyle='->'))
#mpp.arrow(-8.1,500,0,100,color='k', head_width=.01, head_length=25,length_includes_head=True)
mpp.text(-8.09,550,r"$w'$")
# u' arrow 1 DONE
mpp.annotate('',xy=(-8+.27,750),xytext=(-8,750),arrowprops=dict(arrowstyle='->'))
#mpp.arrow(-8.0,750,0.26,0,color='k', head_width=20, head_length=.02,length_includes_head=True)
mpp.text(-7.9,710,r"$u'>0$")
# u' arrow 2 DONE
mpp.annotate('',xy=(-7.968+.12,900),xytext=(-7.968,900),arrowprops=dict(arrowstyle='->'))
#mpp.arrow(-7.965,900,0.115,0,color='k', head_width=20, head_length=.02,length_includes_head=True)
mpp.text(-7.94,860,r"$u'>0$")
# gradient 1
mpp.annotate('',xy=(-7.995-.03,790),xytext=(-7.995,710),arrowprops=dict(arrowstyle='->'))
#mpp.arrow(-7.995,710,-.03,80,color='k', head_width=.01, head_length=25,length_includes_head=True)
mpp.text(-8.07,730,r"$\frac{d \bar{u}}{dz} < 0$")
# gradient 2
mpp.annotate('',xy=(-7.995+0.04,880+60),xytext=(-7.995,880),arrowprops=dict(arrowstyle='->'))
#mpp.arrow(-7.99,880,.04,60,color='k', head_width=.01, head_length=25,length_includes_head=True)
mpp.text(-8.04,910,r"$\frac{d \bar{u}}{dz} > 0$")

rect = Rectangle((-8.2,800),.7,200,color='grey',alpha=.5)
mpp.gca().add_patch(rect)
mpp.text(-7.75,920,"Region of upgradient",color='black')
mpp.text(-7.75,880,"momentum flux",color='black')


pdf.savefig()

pdf.close()