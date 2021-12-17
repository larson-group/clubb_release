## python script to generate figures for Stephens, Larson, and Mironov 2022:  ``A parameterization 
## of non-hydrostatic pressure transport terms in second- and third-order turbulence equations''


import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset



#folders---USER MODIFIES THIS
clubb_folder  = '/path/' #specify location of CLUBB "new" and "dflt" files 
sam_folder = '/path/' #specify location of SAM files
output_folder = '/path/' #specify output folder for figures



#load files---need SAM and CLUBB files for BOMEX, RF01, and WANGARA
cbztd=Dataset(clubb_folder+'bomex_zt_dflt.nc')
cbzmd=Dataset(clubb_folder+'bomex_zm_dflt.nc')
cbzt=Dataset(clubb_folder+'bomex_zt_new.nc')
cbzm=Dataset(clubb_folder+'bomex_zm_new.nc')
sb=Dataset(sam_folder+'BOMEX_64x64x75_100m_40m_1s.nc')

crztd=Dataset(clubb_folder+'dycoms2_rf01_zt_dflt.nc')
crzmd=Dataset(clubb_folder+'dycoms2_rf01_zm_dflt.nc')
crzt=Dataset(clubb_folder+'dycoms2_rf01_zt_new.nc')
crzm=Dataset(clubb_folder+'dycoms2_rf01_zm_new.nc')
sr=Dataset(sam_folder+'DYCOMS_RF01_96x96x320.nc')

cwztd=Dataset(clubb_folder+'wangara_zt_dflt.nc')
cwzmd=Dataset(clubb_folder+'wangara_zm_dflt.nc')
cwzt=Dataset(clubb_folder+'wangara_zt_new.nc')
cwzm=Dataset(clubb_folder+'wangara_zm_new.nc')
sw=Dataset(sam_folder+'WANGARA_64x64x80_100m_40m_LES.nc')


#define variables

###clubb
cbhztd=np.array(cbztd.variables['altitude'])
cbhzmd=np.array(cbzmd.variables['altitude'])
crhztd=np.array(crztd.variables['altitude'])
crhzmd=np.array(crzmd.variables['altitude'])
cwhztd=np.array(cwztd.variables['altitude'])
cwhzmd=np.array(cwzmd.variables['altitude'])

cbhzt=np.array(cbzt.variables['altitude'])
cbhzm=np.array(cbzm.variables['altitude'])
crhzt=np.array(crzt.variables['altitude'])
crhzm=np.array(crzm.variables['altitude'])
cwhzt=np.array(cwzt.variables['altitude'])
cwhzm=np.array(cwzm.variables['altitude'])

cbwp2d=np.mean(np.array(cbzmd.variables['wp2']),(2,3))
crwp2d=np.mean(np.array(crzmd.variables['wp2']),(2,3))
cwwp2d=np.mean(np.array(cwzmd.variables['wp2']),(2,3))

cbwp2=np.mean(np.array(cbzm.variables['wp2']),(2,3))
crwp2=np.mean(np.array(crzm.variables['wp2']),(2,3))
cwwp2=np.mean(np.array(cwzm.variables['wp2']),(2,3))

cbwp3d=np.mean(np.array(cbztd.variables['wp3']),(2,3))
crwp3d=np.mean(np.array(crztd.variables['wp3']),(2,3))
cwwp3d=np.mean(np.array(cwztd.variables['wp3']),(2,3))

cbwp3=np.mean(np.array(cbzt.variables['wp3']),(2,3))
crwp3=np.mean(np.array(crzt.variables['wp3']),(2,3))
cwwp3=np.mean(np.array(cwzt.variables['wp3']),(2,3))

#wp2
cbwp2_bt=np.mean(np.array(cbzm.variables['wp2_bt']),(2,3))
cbwp2_ma=np.mean(np.array(cbzm.variables['wp2_ma']),(2,3))
cbwp2_ta=np.mean(np.array(cbzm.variables['wp2_ta']),(2,3))
cbwp2_ac=np.mean(np.array(cbzm.variables['wp2_ac']),(2,3))
cbwp2_bp=np.mean(np.array(cbzm.variables['wp2_bp']),(2,3))
cbwp2_pr1=np.mean(np.array(cbzm.variables['wp2_pr1']),(2,3))
cbwp2_pr2=np.mean(np.array(cbzm.variables['wp2_pr2']),(2,3))
cbwp2_pr3=np.mean(np.array(cbzm.variables['wp2_pr3']),(2,3))
cbwp2_pr_dfsn=np.mean(np.array(cbzm.variables['wp2_pr_dfsn']),(2,3))
cbwp2_dp1=np.mean(np.array(cbzm.variables['wp2_dp1']),(2,3))
cbwp2_dp2=np.mean(np.array(cbzm.variables['wp2_dp2']),(2,3))
cbwp2_cl=np.mean(np.array(cbzm.variables['wp2_cl']),(2,3))
cbwp2_pd=np.mean(np.array(cbzm.variables['wp2_pd']),(2,3))
cbwp2_splat=np.mean(np.array(cbzm.variables['wp2_splat']),(2,3))
cbwp2_sf=np.mean(np.array(cbzm.variables['wp2_sf']),(2,3))
cbwp2_pr_dfsn_tot = cbwp2_pr_dfsn + cbwp2_dp2
cbwp2_pres_scram=cbwp2_pr1+cbwp2_pr2+cbwp2_pr3+cbwp2_splat
cbwp2_adv=cbwp2_ma+cbwp2_ta+cbwp2_ac
cbwp2_buoy=cbwp2_bp
cbwp2_diss=cbwp2_dp1
cbwp2_tot=cbwp2_bt
cbwp2_other=cbwp2_cl+cbwp2_pd+cbwp2_sf
cbwp2_res=cbwp2_bt-(cbwp2_ma+cbwp2_ta+cbwp2_ac+cbwp2_bp+cbwp2_pr1+cbwp2_pr2+cbwp2_pr3+cbwp2_pr_dfsn+cbwp2_dp1+cbwp2_dp2+cbwp2_cl+cbwp2_pd+cbwp2_splat+cbwp2_sf)

crwp2_bt=np.mean(np.array(crzm.variables['wp2_bt']),(2,3))
crwp2_ma=np.mean(np.array(crzm.variables['wp2_ma']),(2,3))
crwp2_ta=np.mean(np.array(crzm.variables['wp2_ta']),(2,3))
crwp2_ac=np.mean(np.array(crzm.variables['wp2_ac']),(2,3))
crwp2_bp=np.mean(np.array(crzm.variables['wp2_bp']),(2,3))
crwp2_pr1=np.mean(np.array(crzm.variables['wp2_pr1']),(2,3))
crwp2_pr2=np.mean(np.array(crzm.variables['wp2_pr2']),(2,3))
crwp2_pr3=np.mean(np.array(crzm.variables['wp2_pr3']),(2,3))
crwp2_pr_dfsn=np.mean(np.array(crzm.variables['wp2_pr_dfsn']),(2,3))
crwp2_dp1=np.mean(np.array(crzm.variables['wp2_dp1']),(2,3))
crwp2_dp2=np.mean(np.array(crzm.variables['wp2_dp2']),(2,3))
crwp2_cl=np.mean(np.array(crzm.variables['wp2_cl']),(2,3))
crwp2_pd=np.mean(np.array(crzm.variables['wp2_pd']),(2,3))
crwp2_splat=np.mean(np.array(crzm.variables['wp2_splat']),(2,3))
crwp2_sf=np.mean(np.array(crzm.variables['wp2_sf']),(2,3))
crwp2_pr_dfsn_tot = crwp2_pr_dfsn + crwp2_dp2
crwp2_pres_scram=crwp2_pr1+crwp2_pr2+crwp2_pr3+crwp2_splat
crwp2_adv=crwp2_ma+crwp2_ta+crwp2_ac
crwp2_buoy=crwp2_bp
crwp2_diss=crwp2_dp1
crwp2_tot=crwp2_bt
crwp2_other=crwp2_cl+crwp2_pd+crwp2_sf
crwp2_res=crwp2_bt-(crwp2_ma+crwp2_ta+crwp2_ac+crwp2_bp+crwp2_pr1+crwp2_pr2+crwp2_pr3+crwp2_pr_dfsn+crwp2_dp1+crwp2_dp2+crwp2_cl+crwp2_pd+crwp2_splat+crwp2_sf)

cwwp2_bt=np.mean(np.array(cwzm.variables['wp2_bt']),(2,3))
cwwp2_ma=np.mean(np.array(cwzm.variables['wp2_ma']),(2,3))
cwwp2_ta=np.mean(np.array(cwzm.variables['wp2_ta']),(2,3))
cwwp2_ac=np.mean(np.array(cwzm.variables['wp2_ac']),(2,3))
cwwp2_bp=np.mean(np.array(cwzm.variables['wp2_bp']),(2,3))
cwwp2_pr1=np.mean(np.array(cwzm.variables['wp2_pr1']),(2,3))
cwwp2_pr2=np.mean(np.array(cwzm.variables['wp2_pr2']),(2,3))
cwwp2_pr3=np.mean(np.array(cwzm.variables['wp2_pr3']),(2,3))
cwwp2_pr_dfsn=np.mean(np.array(cwzm.variables['wp2_pr_dfsn']),(2,3))
cwwp2_dp1=np.mean(np.array(cwzm.variables['wp2_dp1']),(2,3))
cwwp2_dp2=np.mean(np.array(cwzm.variables['wp2_dp2']),(2,3))
cwwp2_cl=np.mean(np.array(cwzm.variables['wp2_cl']),(2,3))
cwwp2_pd=np.mean(np.array(cwzm.variables['wp2_pd']),(2,3))
cwwp2_splat=np.mean(np.array(cwzm.variables['wp2_splat']),(2,3))
cwwp2_sf=np.mean(np.array(cwzm.variables['wp2_sf']),(2,3))
cwwp2_pr_dfsn_tot = cwwp2_pr_dfsn + cwwp2_dp2
cwwp2_pres_scram=cwwp2_pr1+cwwp2_pr2+cwwp2_pr3+cwwp2_splat
cwwp2_adv=cwwp2_ma+cwwp2_ta+cwwp2_ac
cwwp2_buoy=cwwp2_bp
cwwp2_diss=cwwp2_dp1
cwwp2_tot=cwwp2_bt
cwwp2_other=cwwp2_cl+cwwp2_pd+cwwp2_sf
cwwp2_res=cwwp2_bt-(cwwp2_ma+cwwp2_ta+cwwp2_ac+cwwp2_bp+cwwp2_pr1+cwwp2_pr2+cwwp2_pr3+cwwp2_pr_dfsn+cwwp2_dp1+cwwp2_dp2+cwwp2_cl+cwwp2_pd+cwwp2_splat+cwwp2_sf)

#wp3
cbwp3_bt=np.mean(np.array(cbzt.variables['wp3_bt']),(2,3))
cbwp3_ma=np.mean(np.array(cbzt.variables['wp3_ma']),(2,3))
cbwp3_ta=np.mean(np.array(cbzt.variables['wp3_ta']),(2,3))
cbwp3_ac=np.mean(np.array(cbzt.variables['wp3_ac']),(2,3))
cbwp3_bp1=np.mean(np.array(cbzt.variables['wp3_bp1']),(2,3))
cbwp3_pr1=np.mean(np.array(cbzt.variables['wp3_pr1']),(2,3))
cbwp3_pr2=np.mean(np.array(cbzt.variables['wp3_pr2']),(2,3))
cbwp3_pr3=np.mean(np.array(cbzt.variables['wp3_pr3']),(2,3))
cbwp3_pr_tp=np.mean(np.array(cbzt.variables['wp3_pr_tp']),(2,3))
cbwp3_pr_dfsn=np.mean(np.array(cbzt.variables['wp3_pr_dfsn']),(2,3))
cbwp3_dp1=np.mean(np.array(cbzt.variables['wp3_dp1']),(2,3))
cbwp3_tp=np.mean(np.array(cbzt.variables['wp3_tp']),(2,3))
cbwp3_cl=np.mean(np.array(cbzt.variables['wp3_cl']),(2,3))
cbwp3_splat=np.mean(np.array(cbzt.variables['wp3_splat']),(2,3))
cbwp3_pr_dfsn_tot = cbwp3_pr_dfsn + cbwp3_dp1
cbwp3_pres_scram=cbwp3_pr1+cbwp3_pr2+cbwp3_pr3+cbwp3_pr_tp+cbwp3_splat
cbwp3_adv=cbwp3_ma+cbwp3_ta+cbwp3_ac
cbwp3_buoy=cbwp3_bp1
cbwp3_diss=cbwp3_splat
cbwp3_tot=cbwp3_bt
cbwp3_res=cbwp3_bt-(cbwp3_ma+cbwp3_ta+cbwp3_ac+cbwp3_bp1+cbwp3_pr1+cbwp3_pr2+cbwp3_pr3+cbwp3_pr_tp+cbwp3_pr_dfsn+cbwp3_dp1+cbwp3_tp+cbwp3_cl+cbwp3_splat)

crwp3_bt=np.mean(np.array(crzt.variables['wp3_bt']),(2,3))
crwp3_ma=np.mean(np.array(crzt.variables['wp3_ma']),(2,3))
crwp3_ta=np.mean(np.array(crzt.variables['wp3_ta']),(2,3))
crwp3_ac=np.mean(np.array(crzt.variables['wp3_ac']),(2,3))
crwp3_bp1=np.mean(np.array(crzt.variables['wp3_bp1']),(2,3))
crwp3_pr1=np.mean(np.array(crzt.variables['wp3_pr1']),(2,3))
crwp3_pr2=np.mean(np.array(crzt.variables['wp3_pr2']),(2,3))
crwp3_pr3=np.mean(np.array(crzt.variables['wp3_pr3']),(2,3))
crwp3_pr_tp=np.mean(np.array(crzt.variables['wp3_pr_tp']),(2,3))
crwp3_pr_dfsn=np.mean(np.array(crzt.variables['wp3_pr_dfsn']),(2,3))
crwp3_dp1=np.mean(np.array(crzt.variables['wp3_dp1']),(2,3))
crwp3_tp=np.mean(np.array(crzt.variables['wp3_tp']),(2,3))
crwp3_cl=np.mean(np.array(crzt.variables['wp3_cl']),(2,3))
crwp3_splat=np.mean(np.array(crzt.variables['wp3_splat']),(2,3))
crwp3_pr_dfsn_tot = crwp3_pr_dfsn + crwp3_dp1
crwp3_pres_scram=crwp3_pr1+crwp3_pr2+crwp3_pr3+crwp3_pr_tp+crwp3_splat
crwp3_adv=crwp3_ma+crwp3_ta+crwp3_ac
crwp3_buoy=crwp3_bp1
crwp3_diss=crwp3_splat
crwp3_tot=crwp3_bt
crwp3_res=crwp3_bt-(crwp3_ma+crwp3_ta+crwp3_ac+crwp3_bp1+crwp3_pr1+crwp3_pr2+crwp3_pr3+crwp3_pr_tp+crwp3_pr_dfsn+crwp3_dp1+crwp3_tp+crwp3_cl+crwp3_splat)

cwwp3_bt=np.mean(np.array(cwzt.variables['wp3_bt']),(2,3))
cwwp3_ma=np.mean(np.array(cwzt.variables['wp3_ma']),(2,3))
cwwp3_ta=np.mean(np.array(cwzt.variables['wp3_ta']),(2,3))
cwwp3_ac=np.mean(np.array(cwzt.variables['wp3_ac']),(2,3))
cwwp3_bp1=np.mean(np.array(cwzt.variables['wp3_bp1']),(2,3))
cwwp3_pr1=np.mean(np.array(cwzt.variables['wp3_pr1']),(2,3))
cwwp3_pr2=np.mean(np.array(cwzt.variables['wp3_pr2']),(2,3))
cwwp3_pr3=np.mean(np.array(cwzt.variables['wp3_pr3']),(2,3))
cwwp3_pr_tp=np.mean(np.array(cwzt.variables['wp3_pr_tp']),(2,3))
cwwp3_pr_dfsn=np.mean(np.array(cwzt.variables['wp3_pr_dfsn']),(2,3))
cwwp3_dp1=np.mean(np.array(cwzt.variables['wp3_dp1']),(2,3))
cwwp3_tp=np.mean(np.array(cwzt.variables['wp3_tp']),(2,3))
cwwp3_cl=np.mean(np.array(cwzt.variables['wp3_cl']),(2,3))
cwwp3_splat=np.mean(np.array(cwzt.variables['wp3_splat']),(2,3))
cwwp3_pr_dfsn_tot = cwwp3_pr_dfsn + cwwp3_dp1
cwwp3_pres_scram=cwwp3_pr1+cwwp3_pr2+cwwp3_pr3+cwwp3_pr_tp+cwwp3_splat
cwwp3_adv=cwwp3_ma+cwwp3_ta+cwwp3_ac
cwwp3_buoy=cwwp3_bp1
cwwp3_diss=cwwp3_splat
cwwp3_tot=cwwp3_bt
cwwp3_res=cwwp3_bt-(cwwp3_ma+cwwp3_ta+cwwp3_ac+cwwp3_bp1+cwwp3_pr1+cwwp3_pr2+cwwp3_pr3+cwwp3_pr_tp+cwwp3_pr_dfsn+cwwp3_dp1+cwwp3_tp+cwwp3_cl+cwwp3_splat)


###sam-les
sbh=np.array(sb.variables['z'])
srh=np.array(sr.variables['z'])
swh=np.array(sw.variables['z'])

sbup2=np.mean(np.array(sb.variables['U2']),(2,3))
sbvp2=np.mean(np.array(sb.variables['V2']),(2,3))

sbwp2=np.mean(np.array(sb.variables['W2']),(2,3))
sbwp2adv=np.mean(np.array(sb.variables['W2ADV']),(2,3))
sbwp2redis=np.mean(np.array(sb.variables['W2REDIS']),(2,3))
sbwp2buoy=np.mean(np.array(sb.variables['W2BUOY']),(2,3))
sbwp2diff=np.mean(np.array(sb.variables['W2DFSN']),(2,3))
sbwp2bt=np.mean(np.array(sb.variables['W2BT']),(2,3))
sbwp2pres=np.mean(np.array(sb.variables['W2PRES']),(2,3))
sbwp2presdfsn=np.mean(np.array(sb.variables['W2PRESDFSN']),(2,3))
sbwp2res = sbwp2bt - (sbwp2adv + sbwp2pres + sbwp2redis + sbwp2buoy + sbwp2diff)

sbwppp=np.mean(np.array(sb.variables['WPPP']),(2,3))
sbwptke=np.mean(np.array(sb.variables['WPTKE']),(2,3))
sbwp2pp=np.mean(np.array(sb.variables['WP2PP']),(2,3))
sbwp2tke=np.mean(np.array(sb.variables['WP2TKE']),(2,3))
sbwp2tke = sbwp2tke - 0.5*sbup2*sbwp2 - 0.5*sbvp2*sbwp2 - 0.5*sbwp2**2

sbtkeadv=np.mean(np.array(sb.variables['ADVTR']),(2,3))
sbtkeshear=np.mean(np.array(sb.variables['SHEAR']),(2,3))
sbtkebuoy=np.mean(np.array(sb.variables['BUOYA']),(2,3))
sbtkepres=np.mean(np.array(sb.variables['PRESSTR']),(2,3))
sbtkediftr=np.mean(np.array(sb.variables['DIFTR']),(2,3))
sbtkedissip=np.mean(np.array(sb.variables['DISSIP']),(2,3))
sbtkediftrplusdissip=sbtkediftr+sbtkedissip
sbtkebt=np.mean(np.array(sb.variables['BT']),(2,3))
sbtkeres=sbtkebt-(sbtkeshear+sbtkebuoy+sbtkeadv+sbtkepres+sbtkediftrplusdissip)

srup2=np.mean(np.array(sr.variables['U2']),(2,3))
srvp2=np.mean(np.array(sr.variables['V2']),(2,3))

srwp2=np.mean(np.array(sr.variables['W2']),(2,3))
srwp2adv=np.mean(np.array(sr.variables['W2ADV']),(2,3))
srwp2redis=np.mean(np.array(sr.variables['W2REDIS']),(2,3))
srwp2buoy=np.mean(np.array(sr.variables['W2BUOY']),(2,3))
srwp2diff=np.mean(np.array(sr.variables['W2DFSN']),(2,3))
srwp2bt=np.mean(np.array(sr.variables['W2BT']),(2,3))
srwp2pres=np.mean(np.array(sr.variables['W2PRES']),(2,3))
srwp2presdfsn=np.mean(np.array(sr.variables['W2PRESDFSN']),(2,3))
srwp2res = srwp2bt - (srwp2adv + srwp2pres + srwp2redis + srwp2buoy + srwp2diff)

srwppp=np.mean(np.array(sr.variables['WPPP']),(2,3))
srwptke=np.mean(np.array(sr.variables['WPTKE']),(2,3))
srwp2pp=np.mean(np.array(sr.variables['WP2PP']),(2,3))
srwp2tke=np.mean(np.array(sr.variables['WP2TKE']),(2,3))
srwp2tke = srwp2tke - 0.5*srup2*srwp2 - 0.5*srvp2*srwp2 - 0.5*srwp2**2

srtkeadv=np.mean(np.array(sr.variables['ADVTR']),(2,3))
srtkeshear=np.mean(np.array(sr.variables['SHEAR']),(2,3))
srtkebuoy=np.mean(np.array(sr.variables['BUOYA']),(2,3))
srtkepres=np.mean(np.array(sr.variables['PRESSTR']),(2,3))
srtkediftr=np.mean(np.array(sr.variables['DIFTR']),(2,3))
srtkedissip=np.mean(np.array(sr.variables['DISSIP']),(2,3))
srtkediftrplusdissip=srtkediftr+srtkedissip
srtkebt=np.mean(np.array(sr.variables['BT']),(2,3))
srtkeres=srtkebt-(srtkeshear+srtkebuoy+srtkeadv+srtkepres+srtkediftrplusdissip)

swup2=np.mean(np.array(sw.variables['U2']),(2,3))
swvp2=np.mean(np.array(sw.variables['V2']),(2,3))

swwp2=np.mean(np.array(sw.variables['W2']),(2,3))
swwp2adv=np.mean(np.array(sw.variables['W2ADV']),(2,3))
swwp2redis=np.mean(np.array(sw.variables['W2REDIS']),(2,3))
swwp2buoy=np.mean(np.array(sw.variables['W2BUOY']),(2,3))
swwp2diff=np.mean(np.array(sw.variables['W2DFSN']),(2,3))
swwp2bt=np.mean(np.array(sw.variables['W2BT']),(2,3))
swwp2pres=np.mean(np.array(sw.variables['W2PRES']),(2,3))
swwp2presdfsn=np.mean(np.array(sw.variables['W2PRESDFSN']),(2,3))
swwp2res = swwp2bt - (swwp2adv + swwp2pres + swwp2redis + swwp2buoy + swwp2diff)

swwppp=np.mean(np.array(sw.variables['WPPP']),(2,3))
swwptke=np.mean(np.array(sw.variables['WPTKE']),(2,3))
swwp2pp=np.mean(np.array(sw.variables['WP2PP']),(2,3))
swwp2tke=np.mean(np.array(sw.variables['WP2TKE']),(2,3))
swwp2tke = swwp2tke - 0.5*swup2*swwp2 - 0.5*swvp2*swwp2 - 0.5*swwp2**2

swtkeadv=np.mean(np.array(sw.variables['ADVTR']),(2,3))
swtkeshear=np.mean(np.array(sw.variables['SHEAR']),(2,3))
swtkebuoy=np.mean(np.array(sw.variables['BUOYA']),(2,3))
swtkepres=np.mean(np.array(sw.variables['PRESSTR']),(2,3))
swtkediftr=np.mean(np.array(sw.variables['DIFTR']),(2,3))
swtkedissip=np.mean(np.array(sw.variables['DISSIP']),(2,3))
swtkediftrplusdissip=swtkediftr+swtkedissip
swtkebt=np.mean(np.array(sw.variables['BT']),(2,3))
swtkeres=swtkebt-(swtkeshear+swtkebuoy+swtkeadv+swtkepres+swtkediftrplusdissip)

sbwp3=np.mean(np.array(sb.variables['W3']),(2,3))
sbwp3tp=np.mean(np.array(sb.variables['W3TP']),(2,3))
sbwp3adv=np.mean(np.array(sb.variables['W3ADV']),(2,3))
sbwp3pres=np.mean(np.array(sb.variables['W3PRES']),(2,3))
sbwp3buoy=np.mean(np.array(sb.variables['W3BUOY']),(2,3))
sbwp3diff=np.mean(np.array(sb.variables['W3DFSN']),(2,3))
sbwp3bt=np.mean(np.array(sb.variables['W3BT']),(2,3))
sbwp3presdfsn=np.mean(np.array(sb.variables['W3PRESDFSN']),(2,3))
sbwp3presscr=np.mean(np.array(sb.variables['W3PRESSCR']),(2,3))
sbwp3pres=sbwp3pres-sbwp3tp
sbwp3res = sbwp3bt - (sbwp3adv + sbwp3pres + sbwp3buoy + sbwp3diff + sbwp3tp)

srwp3=np.mean(np.array(sr.variables['W3']),(2,3))
srwp3tp=np.mean(np.array(sr.variables['W3TP']),(2,3))
srwp3adv=np.mean(np.array(sr.variables['W3ADV']),(2,3))
srwp3pres=np.mean(np.array(sr.variables['W3PRES']),(2,3))
srwp3buoy=np.mean(np.array(sr.variables['W3BUOY']),(2,3))
srwp3diff=np.mean(np.array(sr.variables['W3DFSN']),(2,3))
srwp3bt=np.mean(np.array(sr.variables['W3BT']),(2,3))
srwp3presdfsn=np.mean(np.array(sr.variables['W3PRESDFSN']),(2,3))
srwp3presscr=np.mean(np.array(sr.variables['W3PRESSCR']),(2,3))
srwp3pres=srwp3pres-srwp3tp
srwp3res = srwp3bt - (srwp3adv + srwp3pres + srwp3buoy + srwp3diff + srwp3tp)

swwp3=np.mean(np.array(sw.variables['W3']),(2,3))
swwp3tp=np.mean(np.array(sw.variables['W3TP']),(2,3))
swwp3adv=np.mean(np.array(sw.variables['W3ADV']),(2,3))
swwp3pres=np.mean(np.array(sw.variables['W3PRES']),(2,3))
swwp3buoy=np.mean(np.array(sw.variables['W3BUOY']),(2,3))
swwp3diff=np.mean(np.array(sw.variables['W3DFSN']),(2,3))
swwp3bt=np.mean(np.array(sw.variables['W3BT']),(2,3))
swwp3presdfsn=np.mean(np.array(sw.variables['W3PRESDFSN']),(2,3))
swwp3presscr=np.mean(np.array(sw.variables['W3PRESSCR']),(2,3))
swwp3pres=swwp3pres-swwp3tp
swwp3res = swwp3bt - (swwp3adv + swwp3pres + swwp3buoy + swwp3diff + swwp3tp)

#define averaging times & heights
bs=180
be=359
bh=70

rs=180
re=239
rhc=100
rhs=200

ws=180
we=239
wh=50


# colors
cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']


# any "label" specifications below are only relevant if ax.legend() is used for that subplot.

##########################################################
##########################################################
#figure testing wp2 and wp3 sensitivity
fig=plt.figure(figsize=(10.5,7.))
#
ax=fig.add_subplot(2,3,1)
ax.plot(np.mean(sbwp2[bs:be,0:bh],0),sbh[0:bh],'k',label='SAM-LES',linewidth=2.)
ax.plot(np.mean(cbwp2[bs:be,0:bh],0),cbhzt[0:bh],'c',label='CLUBB')
ax.plot(np.mean(cbwp2d[bs:be,0:bh],0),cbhztd[0:bh],color='orange',label='CLUBB')
xmin,xmax = ax.get_xlim()
ymin,ymax = ax.get_ylim()
ax.set_xlabel(r"$\overline{w'^2}$ [m$^2$s$^{-2}$]")
ax.set_ylabel('Height [m]')
ax.set_title('BOMEX shallow cumulus')
#
ax=fig.add_subplot(2,3,2)
ax.plot(np.mean(srwp2[rs:re,0:rhs],0),srh[0:rhs],'k',label='SAM-LES',linewidth=2.)
ax.plot(np.mean(crwp2[rs:re,0:rhc],0),crhzt[0:rhc],'c',label='CLUBB')
ax.plot(np.mean(crwp2d[rs:re,0:rhc],0),crhztd[0:rhc],color='orange',label='CLUBB')
xmin,xmax = ax.get_xlim()
ymin,ymax = ax.get_ylim()
ax.set_xlabel(r"$\overline{w'^2}$ [m$^2$s$^{-2}$]")
ax.set_title('DYCOMS2_RF01 stratocumulus')
#
ax=fig.add_subplot(2,3,3)
ax.plot(np.mean(swwp2[ws:we,0:wh],0),swh[0:wh],'k',label='SAM-LES',linewidth=2.)
ax.plot(np.mean(cwwp2[ws:we,0:wh],0),cwhzt[0:wh],'c',label='CLUBB new')
ax.plot(np.mean(cwwp2d[ws:we,0:wh],0),cwhztd[0:wh],color='orange',label='CLUBB dflt.')
xmin,xmax = ax.get_xlim()
ymin,ymax = ax.get_ylim()
ax.set_xlabel(r"$\overline{w'^2}$ [m$^2$s$^{-2}$]")
ax.set_title('Wangara dry convection')
ax.legend(frameon=False)
#
ax=fig.add_subplot(2,3,4)
ax.plot(np.mean(sbwp3[bs:be,0:bh],0),sbh[0:bh],'k',label='SAM-LES',linewidth=2.)
ax.plot(np.mean(cbwp3[bs:be,0:bh],0),cbhzt[0:bh],'c',label='CLUBB')
ax.plot(np.mean(cbwp3d[bs:be,0:bh],0),cbhztd[0:bh],color='orange',label='CLUBB')
xmin,xmax = ax.get_xlim()
ymin,ymax = ax.get_ylim()
ax.set_xlabel(r"$\overline{w'^3}$ [m$^3$s$^{-3}$]")
ax.set_ylabel('Height [m]')
#
ax=fig.add_subplot(2,3,5)
ax.plot(np.mean(srwp3[rs:re,0:rhs],0),srh[0:rhs],'k',label='SAM-LES',linewidth=2.)
ax.plot(np.mean(crwp3[rs:re,0:rhc],0),crhzt[0:rhc],'c',label='CLUBB')
ax.plot(np.mean(crwp3d[rs:re,0:rhc],0),crhztd[0:rhc],color='orange',label='CLUBB')
xmin,xmax = ax.get_xlim()
ymin,ymax = ax.get_ylim()
ax.set_xlabel(r"$\overline{w'^3}$ [m$^3$s$^{-3}$]")
#
ax=fig.add_subplot(2,3,6)
ax.plot(np.mean(swwp3[ws:we,0:wh],0),swh[0:wh],'k',label='SAM-LES',linewidth=2.)
ax.plot(np.mean(cwwp3[ws:we,0:wh],0),cwhzt[0:wh],'c',label='CLUBB')
ax.plot(np.mean(cwwp3d[ws:we,0:wh],0),cwhztd[0:wh],color='orange',label='CLUBB')
xmin,xmax = ax.get_xlim()
ymin,ymax = ax.get_ylim()
ax.set_xlabel(r"$\overline{w'^3}$ [m$^3$s$^{-3}$]")
#

fig.subplots_adjust(hspace=0.325)
plt.savefig(output_folder+'sensitivity.png',dpi=300,bbox_inches="tight")
plt.show()


##########################################################
##########################################################
#figure showing SAM wp2 and wp3 budgets
fig=plt.figure(figsize=(10.5,7.))

ax=fig.add_subplot(2,3,1)
ax.plot(np.mean(sbwp2pres[bs:be,0:bh],0),sbh[0:bh],label='pres. trans.')
ax.plot(np.mean(sbwp2redis[bs:be,0:bh],0),sbh[0:bh],label='pres. scram.')
ax.plot(np.mean(sbwp2adv[bs:be,0:bh],0),sbh[0:bh],label='adv.')
ax.plot(np.mean(sbwp2buoy[bs:be,0:bh],0),sbh[0:bh],label='buoy.')
ax.plot(np.mean(sbwp2diff[bs:be,0:bh],0),sbh[0:bh],label='diss.')
ax.plot(np.mean(sbwp2bt[bs:be,0:bh],0),sbh[0:bh],label='tot.')
ax.plot(np.mean(sbwp2res[bs:be,0:bh],0),sbh[0:bh],label='resid.')
ax.set_ylabel('Height [m]')
ax.set_xlabel(r"$\overline{w'^2}$ budgets [m$^2$s$^{-3}$]")
xmin,xmax = ax.get_xlim()
ymin,ymax = ax.get_ylim()
ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
ax.set_title('BOMEX shallow cumulus')
#
ax=fig.add_subplot(2,3,2)
ax.plot(np.mean(srwp2pres[rs:re,0:rhs],0),srh[0:rhs],label='W2PRES')
ax.plot(np.mean(srwp2redis[rs:re,0:rhs],0),srh[0:rhs],label='W2REDIS')
ax.plot(np.mean(srwp2adv[rs:re,0:rhs],0),srh[0:rhs],label='W2ADV')
ax.plot(np.mean(srwp2buoy[rs:re,0:rhs],0),srh[0:rhs],label='W2BUOY')
ax.plot(np.mean(srwp2diff[rs:re,0:rhs],0),srh[0:rhs],label='W2DFSN')
ax.plot(np.mean(srwp2bt[rs:re,0:rhs],0),srh[0:rhs],label='W2BT')
ax.plot(np.mean(srwp2res[rs:re,0:rhs],0),srh[0:rhs],label='W2_RES')
ax.set_xlabel(r"$\overline{w'^2}$ budgets [m$^2$s$^{-3}$]")
ax.set_xlim(-0.006,0.006)
xmin,xmax = ax.get_xlim()
ymin,ymax = ax.get_ylim()
ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
ax.set_title('DYCOMS2_RF01 stratocumulus')
#
ax=fig.add_subplot(2,3,3)
ax.plot(np.mean(swwp2pres[ws:we,0:wh],0),swh[0:wh],label='pres. trans.')
ax.plot(np.mean(swwp2redis[ws:we,0:wh],0),swh[0:wh],label='pres. scram.')
ax.plot(np.mean(swwp2adv[ws:we,0:wh],0),swh[0:wh],label='3rd-ord. trans.')
ax.plot(np.mean(swwp2buoy[ws:we,0:wh],0),swh[0:wh],label='buoy. prod.')
ax.plot(np.mean(swwp2diff[ws:we,0:wh],0),swh[0:wh],label='diss.')
ax.plot(np.mean(swwp2bt[ws:we,0:wh],0),swh[0:wh],label='tot. tend.')
ax.plot(np.mean(swwp2res[ws:we,0:wh],0),swh[0:wh],label='resid.')
ax.set_xlabel(r"$\overline{w'^2}$ budgets [m$^2$s$^{-3}$]")
ax.set_xlim(-0.0125,0.0125)
xmin,xmax = ax.get_xlim()
ymin,ymax = ax.get_ylim()
ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
ax.set_title('Wangara dry convection')
ax.legend(frameon=False,loc='upper left',prop={'size': 7})

ax=fig.add_subplot(2,3,4)
ax.plot(np.mean(sbwp3presdfsn[bs:be,0:bh],0),sbh[0:bh],label='W3PRESDFSN')
ax.plot(np.mean(sbwp3presscr[bs:be,0:bh],0),sbh[0:bh],label='W3PRESSCR')
ax.plot(np.mean(sbwp3tp[bs:be,0:bh],0),sbh[0:bh],color=cycle[9],label='W3TP')
ax.plot(np.mean(sbwp3adv[bs:be,0:bh],0),sbh[0:bh],color=cycle[2],label='W3ADV')
ax.plot(np.mean(sbwp3buoy[bs:be,0:bh],0),sbh[0:bh],color=cycle[3],label='W3BUOY')
ax.plot(np.mean(sbwp3diff[bs:be,0:bh],0),sbh[0:bh],color=cycle[4],label='W3DFSN')
ax.plot(np.mean(sbwp3bt[bs:be,0:bh],0),sbh[0:bh],color=cycle[5],label='W3BT')
ax.plot(np.mean(sbwp3res[bs:be,0:bh],0),sbh[0:bh],color=cycle[6],label='W3_RES')
ax.set_ylabel('Height [m]')
ax.set_xlabel(r"$\overline{w'^3}$ budgets [m$^3$s$^{-4}$]")
xmin,xmax = ax.get_xlim()
ymin,ymax = ax.get_ylim()
ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
#
ax=fig.add_subplot(2,3,5)
ax.plot(np.mean(srwp3presdfsn[rs:re,0:rhs],0),srh[0:rhs],label='W3PRESDFSN')
ax.plot(np.mean(srwp3presscr[rs:re,0:rhs],0),srh[0:rhs],label='W3PRESSCR')
ax.plot(np.mean(srwp3tp[rs:re,0:rhs],0),srh[0:rhs],color=cycle[9],label='W3TP')
ax.plot(np.mean(srwp3adv[rs:re,0:rhs],0),srh[0:rhs],color=cycle[2],label='W3ADV')
ax.plot(np.mean(srwp3buoy[rs:re,0:rhs],0),srh[0:rhs],color=cycle[3],label='W3BUOY')
ax.plot(np.mean(srwp3diff[rs:re,0:rhs],0),srh[0:rhs],color=cycle[4],label='W3DFSN')
ax.plot(np.mean(srwp3bt[rs:re,0:rhs],0),srh[0:rhs],color=cycle[5],label='W3BT')
ax.plot(np.mean(srwp3res[rs:re,0:rhs],0),srh[0:rhs],color=cycle[6],label='W3_RES')
ax.set_xlabel(r"$\overline{w'^3}$ budgets [m$^3$s$^{-4}$]")
xmin,xmax = ax.get_xlim()
ymin,ymax = ax.get_ylim()
ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
#
ax=fig.add_subplot(2,3,6)
ax.plot(np.mean(swwp3presdfsn[ws:we,0:wh],0),swh[0:wh],label="\"pres. trans.\"")
ax.plot(np.mean(swwp3presscr[ws:we,0:wh],0),swh[0:wh],label="\"pres. scram.\"")
ax.plot(np.mean(swwp3tp[ws:we,0:wh],0),swh[0:wh],color=cycle[9],label='turb. prod.')
ax.plot(np.mean(swwp3adv[ws:we,0:wh],0),swh[0:wh],color=cycle[2],label='4th-ord. trans.')
ax.plot(np.mean(swwp3buoy[ws:we,0:wh],0),swh[0:wh],color=cycle[3],label='buoy. prod.')
ax.plot(np.mean(swwp3diff[ws:we,0:wh],0),swh[0:wh],color=cycle[4],label='diss.')
ax.plot(np.mean(swwp3bt[ws:we,0:wh],0),swh[0:wh],color=cycle[5],label='tot. tend.')
ax.plot(np.mean(swwp3res[ws:we,0:wh],0),swh[0:wh],color=cycle[6],label='resid.')
ax.set_xlabel(r"$\overline{w'^3}$ budgets [m$^3$s$^{-4}$]")
xmin,xmax = ax.get_xlim()
ymin,ymax = ax.get_ylim()
ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
ax.legend(frameon=False,loc='upper left',prop={'size': 7})

fig.subplots_adjust(hspace=0.325)

plt.savefig(output_folder+'samwp23budgets.png',dpi=300,bbox_inches="tight")
plt.show()



##########################################################
##########################################################
#figure showing CLUBB wp2 and wp3 budgets
fig=plt.figure(figsize=(10.5,7.))

ax=fig.add_subplot(2,3,1)
ax.plot(np.mean(cbwp2_pr_dfsn_tot[bs:be,0:bh],0),cbhzm[0:bh],label='pres. trans.')
ax.plot(np.mean(cbwp2_pres_scram[bs:be,0:bh],0),cbhzm[0:bh],label='pres. scram.')
ax.plot(np.mean(cbwp2_adv[bs:be,0:bh],0),cbhzm[0:bh],label='adv.')
ax.plot(np.mean(cbwp2_buoy[bs:be,0:bh],0),cbhzm[0:bh],label='buoy.')
ax.plot(np.mean(cbwp2_diss[bs:be,0:bh],0),cbhzm[0:bh],label='diss.')
ax.plot(np.mean(cbwp2_tot[bs:be,0:bh],0),cbhzm[0:bh],label='tot.')
ax.plot(np.mean(cbwp2_res[bs:be,0:bh],0),cbhzm[0:bh],label='resid.')
ax.plot(np.mean(cbwp2_other[bs:be,0:bh],0),cbhzm[0:bh],label='other')
ax.set_ylabel('Height [m]')
ax.set_xlabel(r"$\overline{w'^2}$ budgets [m$^2$s$^{-3}$]")
xmin,xmax = ax.get_xlim()
ymin,ymax = ax.get_ylim()
ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
ax.set_title('BOMEX shallow cumulus')
#
ax=fig.add_subplot(2,3,2)
ax.plot(np.mean(crwp2_pr_dfsn_tot[rs:re,0:rhc],0),crhzm[0:rhc],label='W2PRES')
ax.plot(np.mean(crwp2_pres_scram[rs:re,0:rhc],0),crhzm[0:rhc],label='W2REDIS')
ax.plot(np.mean(crwp2_adv[rs:re,0:rhc],0),crhzm[0:rhc],label='W2ADV')
ax.plot(np.mean(crwp2_buoy[rs:re,0:rhc],0),crhzm[0:rhc],label='W2BUOY')
ax.plot(np.mean(crwp2_diss[rs:re,0:rhc],0),crhzm[0:rhc],label='W2DFSN')
ax.plot(np.mean(crwp2_tot[rs:re,0:rhc],0),crhzm[0:rhc],label='W2BT')
ax.plot(np.mean(crwp2_res[rs:re,0:rhc],0),crhzm[0:rhc],label='W2_RES')
ax.plot(np.mean(crwp2_other[rs:re,0:rhc],0),crhzm[0:rhc],label='W2_RES')
ax.set_xlabel(r"$\overline{w'^2}$ budgets [m$^2$s$^{-3}$]")
xmin,xmax = ax.get_xlim()
ymin,ymax = ax.get_ylim()
ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
ax.set_title('DYCOMS2_RF01 stratocumulus')
#
ax=fig.add_subplot(2,3,3)
ax.plot(np.mean(cwwp2_pr_dfsn_tot[ws:we,0:wh],0),cwhzm[0:wh],label='pres. trans.')
ax.plot(np.mean(cwwp2_pres_scram[ws:we,0:wh],0),cwhzm[0:wh],label='pres. scram.')
ax.plot(np.mean(cwwp2_adv[ws:we,0:wh],0),cwhzm[0:wh],label='3rd-ord. trans.')
ax.plot(np.mean(cwwp2_buoy[ws:we,0:wh],0),cwhzm[0:wh],label='buoy. prod.')
ax.plot(np.mean(cwwp2_diss[ws:we,0:wh],0),cwhzm[0:wh],label='diss.')
ax.plot(np.mean(cwwp2_tot[ws:we,0:wh],0),cwhzm[0:wh],label='tot. tend.')
ax.plot(np.mean(cwwp2_res[ws:we,0:wh],0),cwhzm[0:wh],label='resid.')
ax.plot(np.mean(cwwp2_other[ws:we,0:wh],0),cwhzm[0:wh],label='other')
ax.set_xlabel(r"$\overline{w'^2}$ budgets [m$^2$s$^{-3}$]")
xmin,xmax = ax.get_xlim()
ymin,ymax = ax.get_ylim()
ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
ax.set_title('Wangara dry convection')
ax.legend(frameon=False,loc='upper right',prop={'size': 7})

ax=fig.add_subplot(2,3,4)
ax.plot(np.mean(cbwp3_pr_dfsn_tot[bs:be,0:bh],0),cbhzt[0:bh],label='W3PRESDFSN')
ax.plot(np.mean(cbwp3_pres_scram[bs:be,0:bh],0),cbhzt[0:bh],label='W3PRESSCR')
ax.plot(np.mean(cbwp3_tp[bs:be,0:bh],0),cbhzt[0:bh],color=cycle[9],label='W3TP')
ax.plot(np.mean(cbwp3_adv[bs:be,0:bh],0),cbhzt[0:bh],color=cycle[2],label='W3ADV')
ax.plot(np.mean(cbwp3_buoy[bs:be,0:bh],0),cbhzt[0:bh],color=cycle[3],label='W3BUOY')
#ax.plot(np.mean(cbwp3_diss[bs:be,0:bh],0),cbhzt[0:bh],color=cycle[4],label='W3DFSN')
ax.plot(np.mean(cbwp3_tot[bs:be,0:bh],0),cbhzt[0:bh],color=cycle[5],label='W3BT')
ax.plot(np.mean(cbwp3_res[bs:be,0:bh],0),cbhzt[0:bh],color=cycle[6],label='W3_RES')
ax.plot(np.mean(cbwp3_cl[bs:be,0:bh],0),cbhzt[0:bh],color=cycle[7],label='W3_RES')
ax.set_ylabel('Height [m]')
ax.set_xlabel(r"$\overline{w'^3}$ budgets [m$^3$s$^{-4}$]")
xmin,xmax = ax.get_xlim()
ymin,ymax = ax.get_ylim()
ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
#
ax=fig.add_subplot(2,3,5)
ax.plot(np.mean(crwp3_pr_dfsn_tot[rs:re,0:rhc],0),crhzt[0:rhc],label='W3PRESDFSN')
ax.plot(np.mean(crwp3_pres_scram[rs:re,0:rhc],0),crhzt[0:rhc],label='W3PRESSCR')
ax.plot(np.mean(crwp3_tp[rs:re,0:rhc],0),crhzt[0:rhc],color=cycle[9],label='W3TP')
ax.plot(np.mean(crwp3_adv[rs:re,0:rhc],0),crhzt[0:rhc],color=cycle[2],label='W3ADV')
ax.plot(np.mean(crwp3_buoy[rs:re,0:rhc],0),crhzt[0:rhc],color=cycle[3],label='W3BUOY')
#ax.plot(np.mean(crwp3_diss[rs:re,0:rhc],0),crhzt[0:rhc],color=cycle[4],label='W3DFSN')
ax.plot(np.mean(crwp3_tot[rs:re,0:rhc],0),crhzt[0:rhc],color=cycle[5],label='W3BT')
ax.plot(np.mean(crwp3_res[rs:re,0:rhc],0),crhzt[0:rhc],color=cycle[6],label='W3_RES')
ax.plot(np.mean(crwp3_cl[rs:re,0:rhc],0),crhzt[0:rhc],color=cycle[7],label='W3_RES')
ax.set_xlabel(r"$\overline{w'^3}$ budgets [m$^3$s$^{-4}$]")
xmin,xmax = ax.get_xlim()
ymin,ymax = ax.get_ylim()
ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
#
ax=fig.add_subplot(2,3,6)
ax.plot(np.mean(cwwp3_pr_dfsn_tot[ws:we,0:wh],0),cwhzt[0:wh],label="\"pres. trans.\"")
ax.plot(np.mean(cwwp3_pres_scram[ws:we,0:wh],0),cwhzt[0:wh],label="\"pres. scram.\"")
ax.plot(np.mean(cwwp3_tp[ws:we,0:wh],0),cwhzt[0:wh],color=cycle[9],label="turb. prod.")
ax.plot(np.mean(cwwp3_adv[ws:we,0:wh],0),cwhzt[0:wh],color=cycle[2],label='4th-ord. trans.')
ax.plot(np.mean(cwwp3_buoy[ws:we,0:wh],0),cwhzt[0:wh],color=cycle[3],label='buoy. prod.')
#ax.plot(np.mean(cwwp3_diss[ws:we,0:wh],0),cwhzt[0:wh],color=cycle[4],label='diss.')
ax.plot(np.mean(cwwp3_tot[ws:we,0:wh],0),cwhzt[0:wh],color=cycle[5],label='tot. tend.')
ax.plot(np.mean(cwwp3_res[ws:we,0:wh],0),cwhzt[0:wh],color=cycle[6],label='resid.')
ax.plot(np.mean(cwwp3_cl[ws:we,0:wh],0),cwhzt[0:wh],color=cycle[7],label='other')
ax.set_xlabel(r"$\overline{w'^3}$ budgets [m$^3$s$^{-4}$]")
xmin,xmax = ax.get_xlim()
ymin,ymax = ax.get_ylim()
ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
ax.legend(frameon=False,loc='upper right',prop={'size': 7})

fig.subplots_adjust(hspace=0.325)

plt.savefig(output_folder+'clubbwp23budgets.png',dpi=300,bbox_inches="tight")
plt.show()



##########################################################
##########################################################
#figure showing SAM tke budgets
fig=plt.figure(figsize=(10.5,3.5))

ax=fig.add_subplot(1,3,1)
ax.plot(np.mean(sbtkepres[bs:be,0:bh],0),sbh[0:bh],label='PRES')
ax.plot(np.mean(sbtkeshear[bs:be,0:bh],0),sbh[0:bh],color=cycle[8],label='SHEAR')
ax.plot(np.mean(sbtkeadv[bs:be,0:bh],0),sbh[0:bh],color=cycle[2],label='ADV')
ax.plot(np.mean(sbtkebuoy[bs:be,0:bh],0),sbh[0:bh],color=cycle[3],label='BUOY')
ax.plot(np.mean(sbtkediftrplusdissip[bs:be,0:bh],0),sbh[0:bh],color=cycle[4],label='DSSP')
ax.plot(np.mean(sbtkebt[bs:be,0:bh],0),sbh[0:bh],color=cycle[5],label='BT')
ax.plot(np.mean(sbtkeres[bs:be,0:bh],0),sbh[0:bh],color=cycle[6],label='RES')
ax.set_ylabel('Height [m]')
ax.set_xlabel(r"TKE budgets [m$^2$s$^{-3}$]")
xmin,xmax = ax.get_xlim()
ymin,ymax = ax.get_ylim()
ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
ax.set_title('BOMEX shallow cumulus')
plt.text((xmax-xmin)*0.85+xmin,(ymax-ymin)*0.90+ymin,'(a)')
#
ax=fig.add_subplot(1,3,2)
ax.plot(np.mean(srtkepres[rs:re,0:rhs],0),srh[0:rhs],label='W2PRES')
ax.plot(np.mean(srtkeshear[rs:re,0:rhs],0),srh[0:rhs],color=cycle[8],label='W2REDIS')
ax.plot(np.mean(srtkeadv[rs:re,0:rhs],0),srh[0:rhs],color=cycle[2],label='W2ADV')
ax.plot(np.mean(srtkebuoy[rs:re,0:rhs],0),srh[0:rhs],color=cycle[3],label='W2BUOY')
ax.plot(np.mean(srtkediftrplusdissip[rs:re,0:rhs],0),srh[0:rhs],color=cycle[4],label='W2DFSN')
ax.plot(np.mean(srtkebt[rs:re,0:rhs],0),srh[0:rhs],color=cycle[5],label='W2BT')
ax.plot(np.mean(srtkeres[rs:re,0:rhs],0),srh[0:rhs],color=cycle[6],label='W2_RES')
ax.set_xlabel(r"TKE budgets [m$^2$s$^{-3}$]")
ax.set_xlim(-0.003,0.003)
xmin,xmax = ax.get_xlim()
ax.set_title('DYCOMS2_RF01 stratocumulus')
ymin,ymax = ax.get_ylim()
ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.text((xmax-xmin)*0.85+xmin,(ymax-ymin)*0.90+ymin,'(b)')
#
ax=fig.add_subplot(1,3,3)
ax.plot(np.mean(swtkepres[ws:we,0:wh],0),swh[0:wh],label='pres. trans.')
ax.plot(np.mean(swtkeshear[ws:we,0:wh],0),swh[0:wh],color=cycle[8],label='shear prod.')
ax.plot(np.mean(swtkeadv[ws:we,0:wh],0),swh[0:wh],color=cycle[2],label='3rd-ord. trans.')
ax.plot(np.mean(swtkebuoy[ws:we,0:wh],0),swh[0:wh],color=cycle[3],label='buoy. prod.')
ax.plot(np.mean(swtkediftrplusdissip[ws:we,0:wh],0),swh[0:wh],color=cycle[4],label='diss.')
ax.plot(np.mean(swtkebt[ws:we,0:wh],0),swh[0:wh],color=cycle[5],label='tot. tend.')
ax.plot(np.mean(swtkeres[ws:we,0:wh],0),swh[0:wh],color=cycle[6],label='resid.')
ax.set_xlabel(r"TKE budgets [m$^2$s$^{-3}$]")
xmin,xmax = ax.get_xlim()
ymin,ymax = ax.get_ylim()
ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
ax.set_title('Wangara dry convection')
plt.text((xmax-xmin)*0.85+xmin,(ymax-ymin)*0.90+ymin,'(c)')
ax.legend(frameon=False,loc='upper left',prop={'size': 8})

plt.subplots_adjust(bottom=0.2)

plt.savefig(output_folder+'samtkebudgets.png',dpi=300,bbox_inches="tight")
plt.show()

##########################################################
##########################################################
#figure showing SAM breakdown of wp2 presure terms
fig=plt.figure(figsize=(10.5,7.))

ax=fig.add_subplot(2,3,1)
ax.plot(np.mean(sbwp2pres[bs:be,0:bh]+sbwp2redis[bs:be,0:bh],0),sbh[0:bh],'k',linewidth=2.,label='W2PRESTOT')
ax.plot(np.mean(sbwp2pres[bs:be,0:bh],0),sbh[0:bh],label='W2PRESDFSN')
ax.plot(np.mean(sbwp2redis[bs:be,0:bh],0),sbh[0:bh],label='W2PRESSCR')
ax.set_ylabel('Height [m]')
ax.set_xlabel(r"[m$^2$s$^{-3}$]")
xmin,xmax = ax.get_xlim()
ymin,ymax = ax.get_ylim()
ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
ax.set_title('BOMEX shallow cumulus')
#
ax=fig.add_subplot(2,3,2)
ax.plot(np.mean(srwp2pres[rs:re,0:rhs]+srwp2redis[rs:re,0:rhs],0),srh[0:rhs],'k',linewidth=2.,label='W2PRESTOT')
ax.plot(np.mean(srwp2pres[rs:re,0:rhs],0),srh[0:rhs],label='W2PRESDFSN')
ax.plot(np.mean(srwp2redis[rs:re,0:rhs],0),srh[0:rhs],label='W2PRESSCR')
ax.set_xlabel(r"[m$^2$s$^{-3}$]")
xmin,xmax = ax.get_xlim()
ax.set_title('DYCOMS2_RF01 stratocumulus')
ymin,ymax = ax.get_ylim()
ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
#
ax=fig.add_subplot(2,3,3)
ax.plot(np.mean(swwp2pres[ws:we,0:wh]+swwp2redis[ws:we,0:wh],0),swh[0:wh],'k',linewidth=2.,label='total pres. term')
ax.plot(np.mean(swwp2pres[ws:we,0:wh],0),swh[0:wh],label='pres. trans.')
ax.plot(np.mean(swwp2redis[ws:we,0:wh],0),swh[0:wh],label='pres. scram.')
ax.set_xlabel(r"[m$^2$s$^{-3}$]")
xmin,xmax = ax.get_xlim()
ymin,ymax = ax.get_ylim()
ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
ax.set_title('Wangara dry convection')
ax.legend(frameon=False,prop={'size': 8},loc='upper left')

ax=fig.add_subplot(2,3,4)
ax.plot(np.mean(sbwp3pres[bs:be,0:bh],0),sbh[0:bh],'k',linewidth=2.,label='W3PRES')
ax.plot(np.mean(sbwp3presdfsn[bs:be,0:bh],0),sbh[0:bh],label='W3PRESDFSN')
ax.plot(np.mean(sbwp3presscr[bs:be,0:bh],0),sbh[0:bh],label='W3PRESSCR')
ax.set_ylabel('Height [m]')
ax.set_xlabel(r"[m$^3$s$^{-4}$]")
xmin,xmax = ax.get_xlim()
ymin,ymax = ax.get_ylim()
#
ax=fig.add_subplot(2,3,5)
ax.plot(np.mean(srwp3pres[rs:re,0:rhs],0),srh[0:rhs],'k',linewidth=2.,label='W3PRES')
ax.plot(np.mean(srwp3presdfsn[rs:re,0:rhs],0),srh[0:rhs],label='W3PRESDFSN')
ax.plot(np.mean(srwp3presscr[rs:re,0:rhs],0),srh[0:rhs],label='W3PRESSCR')
ax.set_xlabel(r"[m$^3$s$^{-4}$]")
xmin,xmax = ax.get_xlim()
ymin,ymax = ax.get_ylim()
#
ax=fig.add_subplot(2,3,6)
ax.plot(np.mean(swwp3pres[ws:we,0:wh],0),swh[0:wh],'k',linewidth=2.,label='total pres. term')
ax.plot(np.mean(swwp3presdfsn[ws:we,0:wh],0),swh[0:wh],label="\"pres. trans.\"")
ax.plot(np.mean(swwp3presscr[ws:we,0:wh],0),swh[0:wh],label="\"pres. scram.\"")
ax.set_xlabel(r"[m$^3$s$^{-4}$]")
xmin,xmax = ax.get_xlim()
ymin,ymax = ax.get_ylim()
ax.legend(frameon=False,loc='upper right',prop={'size': 8})

fig.subplots_adjust(hspace=0.325)

plt.savefig(output_folder+'samwp23presdecomp.png',dpi=300,bbox_inches="tight")
plt.show()


##########################################################
##########################################################
#figure showing SAM w'p' and w'tke etc. for BOMEX
coef1=-0.20
coef2=-0.15

fig=plt.figure(figsize=(7.,3.5))

ax=fig.add_subplot(1,2,1)
ax.plot(np.mean(sbwppp[bs:be,0:bh],0),sbh[0:bh],label=r"$\overline{w'p'}/\rho$")
ax.plot(coef1*np.mean(sbwptke[bs:be,0:bh],0),sbh[0:bh],label=r"$c_{wp2}\overline{w'e}$")
ax.set_ylabel('Height [m]')
ax.set_xlabel(r"[m$^3$s$^{-3}$]")
xmin,xmax = ax.get_xlim()
ymin,ymax = ax.get_ylim()
ax.legend(frameon=False)

ax=fig.add_subplot(1,2,2)
ax.plot(np.mean(sbwp2pp[bs:be,0:bh],0),sbh[0:bh],label=r"$\overline{w'^2p'}/\rho$")
ax.plot(coef2*np.mean(sbwp2tke[bs:be,0:bh],0),sbh[0:bh],label=r"$c_{wp3}(\overline{w'^2e}-\overline{w'^2}\overline{e})$")
ax.set_xlabel(r"[m$^4$s$^{-4}$]")
xmin,xmax = ax.get_xlim()
ymin,ymax = ax.get_ylim()
ax.legend(frameon=False)

fig.subplots_adjust(hspace=0.325)

plt.savefig(output_folder+'sam_wp2pp_param.png',dpi=300,bbox_inches="tight")
plt.show()


##########################################################
##########################################################
#figure combining pr_dfsn and eddy diffusivity
fig=plt.figure(figsize=(10.5,7.))

ax=fig.add_subplot(2,3,1)
ax.plot(np.mean(sbwp2pres[bs:be,0:bh],0),sbh[0:bh],'k',linewidth=2.,label='W2PRES')
ax.plot(np.mean(cbwp2_pr_dfsn_tot[bs:be,0:bh],0),cbhzm[0:bh],label='tot')
ax.set_ylabel('Height [m]')
ax.set_xlabel(r"[m$^2$s$^{-3}$]")
ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
ax.set_title('BOMEX shallow cumulus')

ax=fig.add_subplot(2,3,2)
ax.plot(np.mean(srwp2pres[rs:re,0:rhs],0),srh[0:rhs],'k',linewidth=2.,label='W2PRES')
ax.plot(np.mean(crwp2_pr_dfsn_tot[rs:re,0:rhc],0),crhzm[0:rhc],label='CLUBB')
ax.set_xlabel(r"[m$^2$s$^{-3}$]")
ax.set_title('DYCOMS2_RF01 stratocumulus')
ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

ax=fig.add_subplot(2,3,3)
ax.plot(np.mean(swwp2pres[ws:we,0:wh],0),swh[0:wh],'k',linewidth=2.,label='SAM pres. trans.')
ax.plot(np.mean(cwwp2_pr_dfsn_tot[ws:we,0:wh],0),cwhzm[0:wh],label='CLUBB pres. trans.')
ax.set_xlabel(r"[m$^2$s$^{-3}$]")
ax.set_title('Wangara dry convection')
ax.legend(frameon=False,prop={'size': 8})
ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

ax=fig.add_subplot(2,3,4)
ax.plot(np.mean(sbwp3presdfsn[bs:be,0:bh],0),sbh[0:bh],'k',linewidth=2.,label='SAM')
ax.plot(np.mean(cbwp3_pr_dfsn_tot[bs:be,0:bh],0),cbhzt[0:bh],label='CLUBB before')
ax.set_ylabel('Height [m]')
ax.set_xlabel(r"[m$^3$s$^{-4}$]")
ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

ax=fig.add_subplot(2,3,5)
ax.plot(np.mean(srwp3presdfsn[rs:re,0:rhs],0),srh[0:rhs],'k',linewidth=2.,label='SAM')
ax.plot(np.mean(crwp3_pr_dfsn_tot[rs:re,0:rhc],0),crhzt[0:rhc],label='CLUBB before')
ax.set_xlabel(r"[m$^3$s$^{-4}$]")
ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

ax=fig.add_subplot(2,3,6)
ax.plot(np.mean(swwp3presdfsn[ws:we,0:wh],0),swh[0:wh],'k',linewidth=2.,label="SAM \"pres. trans.\"")
ax.plot(np.mean(cwwp3_pr_dfsn_tot[ws:we,0:wh],0),cwhzt[0:wh],label="CLUBB \"pres. trans.\"")
ax.set_xlabel(r"[m$^3$s$^{-4}$]")
ax.legend(frameon=False,prop={'size': 8})
ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

fig.subplots_adjust(hspace=0.325)

plt.savefig(output_folder+'samVclubb_pr_dfsn.png',dpi=300,bbox_inches="tight")
plt.show()



##########################################################
##########################################################
#figure with CLUBB pressure transport decomposition into pr_dfsn and eddy diffusivity
fig=plt.figure(figsize=(10.5,7.))

ax=fig.add_subplot(2,3,1)
ax.plot(np.mean(cbwp2_pr_dfsn_tot[bs:be,0:bh],0),cbhzm[0:bh],linewidth=3.,label='total')
ax.plot(np.mean(cbwp2_pr_dfsn[bs:be,0:bh],0),cbhzm[0:bh],label='pr. dfsn.')
ax.plot(np.mean(cbwp2_dp2[bs:be,0:bh],0),cbhzm[0:bh],label='eddy dfsn.')
ax.set_ylabel('Height [m]')
ax.set_xlabel(r"[m$^2$s$^{-3}$]")
ax.set_title('BOMEX shallow cumulus')
ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

ax=fig.add_subplot(2,3,2)
ax.plot(np.mean(crwp2_pr_dfsn_tot[rs:re,0:rhc],0),crhzm[0:rhc],linewidth=3.,label='total')
ax.plot(np.mean(crwp2_pr_dfsn[rs:re,0:rhc],0),crhzm[0:rhc],label='pr. dfsn.')
ax.plot(np.mean(crwp2_dp2[rs:re,0:rhc],0),crhzm[0:rhc],label='eddy dfsn.')
ax.set_xlabel(r"[m$^2$s$^{-3}$]")
ax.set_title('DYCOMS2_RF01 stratocumulus')
ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

ax=fig.add_subplot(2,3,3)
ax.plot(np.mean(cwwp2_pr_dfsn_tot[ws:we,0:wh],0),cwhzm[0:wh],linewidth=3.,label='total')
ax.plot(np.mean(cwwp2_pr_dfsn[ws:we,0:wh],0),cwhzm[0:wh],label='pr. dfsn.')
ax.plot(np.mean(cwwp2_dp2[ws:we,0:wh],0),cwhzm[0:wh],label='eddy dfsn.')
ax.set_xlabel(r"[m$^2$s$^{-3}$]")
ax.set_title('Wangara dry convection')
ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

ax=fig.add_subplot(2,3,4)
ax.plot(np.mean(cbwp3_pr_dfsn_tot[bs:be,0:bh],0),cbhzt[0:bh],linewidth=3.,label='CLUBB before')
ax.plot(np.mean(cbwp3_pr_dfsn[bs:be,0:bh],0),cbhzt[0:bh],label='CLUBB before')
ax.plot(np.mean(cbwp3_dp1[bs:be,0:bh],0),cbhzt[0:bh],label='CLUBB before')
ax.set_ylabel('Height [m]')
ax.set_xlabel(r"[m$^3$s$^{-4}$]")
ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

ax=fig.add_subplot(2,3,5)
ax.plot(np.mean(crwp3_pr_dfsn_tot[rs:re,0:rhc],0),crhzt[0:rhc],linewidth=3.,label='CLUBB before')
ax.plot(np.mean(crwp3_pr_dfsn[rs:re,0:rhc],0),crhzt[0:rhc],label='CLUBB before')
ax.plot(np.mean(crwp3_dp1[rs:re,0:rhc],0),crhzt[0:rhc],label='CLUBB before')
ax.set_xlabel(r"[m$^3$s$^{-4}$]")
ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

ax=fig.add_subplot(2,3,6)
ax.plot(np.mean(cwwp3_pr_dfsn_tot[ws:we,0:wh],0),cwhzt[0:wh],linewidth=3.,label='total pres. trans.')
ax.plot(np.mean(cwwp3_pr_dfsn[ws:we,0:wh],0),cwhzt[0:wh],label='Lumley')
ax.plot(np.mean(cwwp3_dp1[ws:we,0:wh],0),cwhzt[0:wh],label='new eddy dfsn.')
ax.set_xlabel(r"[m$^3$s$^{-4}$]")
ax.legend(frameon=False)
ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

fig.subplots_adjust(hspace=0.325)

plt.savefig(output_folder+'clubb_pr_dfsn_breakdown.png',dpi=300,bbox_inches="tight")
plt.show()
