import netCDF4
import numpy as np
import pylab as pl

clubb_nc       = netCDF4.Dataset('rico_lh_zt.nc')
silhs_sfc_nc   = netCDF4.Dataset('rico_lh_lh_sfc.nc')
silhs_2D_u_nc  = netCDF4.Dataset('rico_lh_u_lh_sample_points_2D.nc')
silhs_2D_l_nc  = netCDF4.Dataset('rico_lh_nl_lh_sample_points_2D.nc')

cf1 = clubb_nc.variables['cloud_frac_1']
cf2 = clubb_nc.variables['cloud_frac_2']
pf1 = clubb_nc.variables['precip_frac_1']
pf2 = clubb_nc.variables['precip_frac_2']
mxf = clubb_nc.variables['mixt_frac']

k_lh_start = silhs_sfc_nc.variables['k_lh_start'][:]

X_mixt_comp = silhs_2D_u_nc.variables['X_mixt_comp'][:]
weights = silhs_2D_u_nc.variables['lh_sample_point_weights'][:]
chi_n = silhs_2D_l_nc.variables['chi'][:]

###################################################
clubb_var = clubb_nc.variables['rrm_auto'][:]
rho = clubb_nc.variables['rho'][:]
l_time_shift = False
################
rr_ln   = silhs_2D_l_nc.variables['rr'][:]
Nr_ln   = silhs_2D_l_nc.variables['Nr'][:]
Ncn_ln   = silhs_2D_l_nc.variables['Ncn'][:]
###################################################

num_samples = rr_ln.shape[2]

time1 = 3000
time2 = 4320

clubb_var_plt = np.empty(time2-time1)

vars_plt = list()
for i in range(0,8):
    vars_plt.append(np.empty(time2-time1))

silhs_var_plt = np.empty(time2-time1)

category_labels = ['c_p_1', 'c_p_2', 'c_np_1', 'c_np_2', 'nc_p_1', 'nc_p_2', 'nc_np_1', 'nc_np_2']

for t in range(time1,time2):
    k = int(round(k_lh_start[t,0,0,0])) - 1
    cf1_t = cf1[t,k,0,0]
    cf2_t = cf2[t,k,0,0]
    pf1_t = pf1[t,k,0,0]
    pf2_t = pf2[t,k,0,0]
    mxf_t = mxf[t,k,0,0]
    if l_time_shift:
        clubb_var_plt[t-time1] = clubb_var[t-1,k,0,0]
    else:
        clubb_var_plt[t-time1] = clubb_var[t,k,0,0]
    for i in range(0,8):
        vars_plt[i][t-time1] = 0.
    silhs_var_plt[t-time1] = 0.
    for i in range(0,num_samples):
        # Determine category elements
        if chi_n[t,k,i,0] > 0.0:
            l_cld = True
        else:
            l_cld = False
        comp = X_mixt_comp[t,k,i,0]
        if comp == 1:
            l_comp_1 = True
        else:
            l_comp_1 = False
        if rr_ln[t,k,i,0] > 0.0:
            l_prc = True
        else:
            l_prc = False
    
        if l_cld and l_prc and l_comp_1:
            j = 0
        elif l_cld and l_prc and not l_comp_1:
            j = 1
        elif l_cld and not l_prc and l_comp_1:
            j = 2
        elif l_cld and not l_prc and not l_comp_1:
            j = 3
        elif not l_cld and l_prc and l_comp_1:
            j = 4
        elif not l_cld and l_prc and not l_comp_1:
            j = 5
        elif not l_cld and not l_prc and l_comp_1:
            j = 6
        elif not l_cld and not l_prc and not l_comp_1:
            j = 7

        # Autoconversion (KK)
        if chi_n[t,k,i,0] > 0.0:
            vars_plt[j][t-time1] = vars_plt[j][t-time1] + (1350.0*(rho[t,k,0,0]/1.0e6)**-1.79) * chi_n[t,k,i,0]**2.47 * Ncn_ln[t,k,i,0]**-1.79 * weights[t,k,i,0]

        # Count
        #vars_plt[j][t-time1] = vars_plt[j][t-time1] + 1

    for i in range(0,8):
        vars_plt[i][t-time1] = vars_plt[i][t-time1] / num_samples
        silhs_var_plt[t-time1] = silhs_var_plt[t-time1] + vars_plt[i][t-time1]

colors = ['red', 'cyan', 'magenta', 'yellow', 'black', 'orange', 'brown', 'pink']

pl.figure()
pl.plot(range(time1,time2), clubb_var_plt[:], label='analytic')
pl.plot(range(time1,time2), silhs_var_plt[:], label='SILHS-sum')
for u in range(0,len(vars_plt)):
    if (u==7):
        pl.plot(range(time1,time2), vars_plt[u][:], colors[u], label=category_labels[u])
    else:
        pl.plot(range(time1,time2), vars_plt[u][:], colors[u], label=category_labels[u])

#clusters = [(0,1,2,3),(4,5,6,7)]
#for u in range(0,len(clusters)):
#    pl.plot(range(time1,time2), vars_plt[clusters[u][0]] + vars_plt[clusters[u][1]] + vars_plt[clusters[u][2]] + vars_plt[clusters[u][3]], label="Cluster "+str(u))

#clusters = [(0,2),(1,3),(4,6),(5,7)]
#for u in range(0,len(clusters)):
#    pl.plot(range(time1,time2), vars_plt[clusters[u][0]] + vars_plt[clusters[u][1]], label="Cluster "+str(u))

pl.legend()
pl.show()
