#
# cam_hole_filling_analysis.py
# 
# Description: Under construction


import netCDF4
import numpy as np
import matplotlib.pyplot as plt

outdir = "./"
out_filename_tavg = "clip_tend_tavg"
out_filename_zavg = "clip_tend_zavg"
out_filetype =".png"


l_show_plots = False

### Function Definitions ###

#----------------------------------------------------------
def get_2d_profile(nc_file, var_name):
#
# Reads a 2d profile (time,lev) from a Scam data file 
#
#----------------------------------------------------------
    import numpy as np

    var = nc_file.variables[var_name]
    var = np.squeeze(var)

    return var
#----------------------------------------------------------

#----------------------------------------------------------
def plot_z_profile(ax, data_label, data, lev,):
#
# Creates a time average of a 2 dimensional (time, lev) array
# and plots the resulting profile. 
#
#----------------------------------------------------------
    import numpy as np

    ax.plot(data,lev,'b-',[0,0],[0,np.max(lev)],'r--')
    ax.set_xlabel(data_label)
    ax.set_ylabel('Altitude [m]')
#----------------------------------------------------------

#----------------------------------------------------------
def plot_t_profile(title, data_label, data, out_filepath):
#
# Creates a time average of a 2 dimensional (time, lev) array
# and plots the resulting profile. 
#
#----------------------------------------------------------
    import numpy as np
    import matplotlib.pyplot as plt

    fig = plt.figure()
    fig.text(.5,.95,title,horizontalalignment='center',)
    plt.plot(range(2001), data)
    plt.ylabel(data_label)
    plt.xlabel('Time [days]')
    fig.savefig(out_filepath)
    plt.close()
#----------------------------------------------------------

#----------------------------------------------------------
def calc_invers_dzt(zm):
#
#
#----------------------------------------------------------
    invers_dzt = [[0 for x in range(nlev)] for x in range(ntimes)]

    for i in range(ntimes):

        for j in range(nlev):
            invers_dzt[i][j] = 1./(zm[i,j+1]-zm[i,j])
        #end for j in range(1,len(zm))

    #end for i in range(1,len(times))

    invers_dzt = np.array(invers_dzt)
    return invers_dzt
#----------------------------------------------------------

#----------------------------------------------------------
def vertical_integral(field, rho_ds, invrs_dz):
#
#
#----------------------------------------------------------
    import numpy as np

    vertical_integral = range(len(field))

    vertical_integral = np.sum( field * rho_ds / invrs_dz )

    return vertical_integral

#----------------------------------------------------------

#----------------------------------------------------------
def vertical_average(field, rho_ds, invrs_dz):
#
#
#----------------------------------------------------------
    import numpy as np

    field_zavg = [0 for x in range(ntimes)]
    dummy_one = [1 for x in range(nlev)]
    
    for i in range(1,ntimes):
    
       numer = vertical_integral(field[i,:], rho, invers_dzt)
       denom = vertical_integral(dummy_one, rho, invers_dzt)
       field_zavg[i] = numer/denom
    # end for i in range(1,ntimes)

    field_zavg = np.array(field_zavg)
    return field_zavg
#----------------------------------------------------------


### Main script starts here ###

# CAM data file
nc_file_path = './camclubb707_L30_T1200.cam.h0.0001-01-01-00000.nc'

nc = netCDF4.Dataset(nc_file_path)

# Grid cell altitudes
lev = nc.variables['lev']
nlev = len(lev)

times = nc.variables['time']
ntimes = len(times)

# Calculate the inverse of each grid cell height (for averaging over z)
rho = get_2d_profile(nc,'RHO_DS_HF')
tmp = [[0 for x in range(nlev)] for x in range(ntimes)]
for i in range(ntimes):
    for j in range(nlev):
        tmp[i][j] = rho[i,j+1]
    #end for
#end for

rho = np.array(tmp)

zm = get_2d_profile(nc,'ZM_HF')
invers_dzt = calc_invers_dzt(zm)

ice_clip_tend = get_2d_profile(nc,'INEGCLPTEND')
ice_clip_tend_tavg = np.average(ice_clip_tend, axis=0)

liq_clip_tend = get_2d_profile(nc,'LNEGCLPTEND')
liq_clip_tend_tavg = np.average(liq_clip_tend, axis=0)

vap_clip_tend = get_2d_profile(nc,'VNEGCLPTEND')
vap_clip_tend_tavg = np.average(vap_clip_tend, axis=0)

# Plot profiles
title='CAM-CLUBB-SILHS clipping budgets'

fig = plt.figure()
fig.text(.5,.95,title,horizontalalignment='center',)

ax1 = plt.subplot2grid((2, 2), (0, 0))
ax2 = plt.subplot2grid((2, 2), (0, 1))
ax3 = plt.subplot2grid((2, 2), (1, 0))
ax4 = plt.subplot2grid((2, 2), (1, 1))

plot_z_profile(ax1,'INEGCLPTEND', ice_clip_tend_tavg,lev)
plot_z_profile(ax2,'LNEGCLPTEND', liq_clip_tend_tavg,lev)
plot_z_profile(ax3,'VNEGCLPTEND', vap_clip_tend_tavg,lev)
plot_z_profile(ax4,'Sum', ice_clip_tend_tavg+liq_clip_tend_tavg+vap_clip_tend_tavg,lev)

fig.savefig(outdir+out_filename_tavg+out_filetype)
plt.close()

# Plot timelines of height averages
levavg_sum = vertical_average(ice_clip_tend+liq_clip_tend+vap_clip_tend, rho, invers_dzt)
levavg_vap = vertical_average(vap_clip_tend, rho, invers_dzt)
levavg_liq = vertical_average(liq_clip_tend, rho, invers_dzt)
levavg_ice = vertical_average(ice_clip_tend, rho, invers_dzt)

print("---- Total ----")
print("Sum: "+str(np.sum(levavg_sum)))
print("Vap: "+str(np.sum(levavg_vap)))
print("Liq: "+str(np.sum(levavg_liq)))
print("Ice: "+str(np.sum(levavg_ice)))

title = 'Sum of Vapor/Ice/Liquid mixing ratio tendencies due to hole filling (total)'
label = 'Mixing ratios (Sum)'
out_filepath = outdir+out_filename_zavg+"_sum"+out_filetype
plot_t_profile(title, label, levavg_sum, out_filepath)

title = 'Vapor mixing ratio tendencies due to hole filling (total)'
label = 'Vapor mixing ratio'
out_filepath = outdir+out_filename_zavg+"_vap"+out_filetype
plot_t_profile(title, label, levavg_vap, out_filepath)

title = 'Cloud liquid mixing ratio tendencies due to hole filling (total)'
label = 'Cloud liquid mixing ratio'
out_filepath = outdir+out_filename_zavg+"_liq"+out_filetype
plot_t_profile(title, label, levavg_liq, out_filepath)

title = 'Cloud ice mixing ratio tendencies due to hole filling (total)'
label = 'Cloud ice mixing ratio'
out_filepath = outdir+out_filename_zavg+"_ice"+out_filetype
plot_t_profile(title, label, levavg_ice, out_filepath)


# Plot profiles of clipping tendencies due to vertical hole filling
ice_clip_tend_vhf = get_2d_profile(nc,'INEGCLPTEND_VHF')
ice_clip_tend_vhf_tavg = np.average(ice_clip_tend_vhf, axis=0)

liq_clip_tend_vhf = get_2d_profile(nc,'LNEGCLPTEND_VHF')
liq_clip_tend_vhf_tavg = np.average(liq_clip_tend_vhf, axis=0)

fig = plt.figure()
fig.text(.5,.95,title,horizontalalignment='center',)

ax1 = plt.subplot2grid((2, 2), (0, 0))
ax2 = plt.subplot2grid((2, 2), (0, 1))
ax3 = plt.subplot2grid((2, 2), (1, 0))

plot_z_profile(ax1,'INEGCLPTEND_VHF', ice_clip_tend_vhf_tavg,lev)
plot_z_profile(ax2,'LNEGCLPTEND_VHF', liq_clip_tend_vhf_tavg,lev)
plot_z_profile(ax3,'Ice+Liq VHF', ice_clip_tend_vhf_tavg+liq_clip_tend_vhf_tavg,lev)

fig.savefig(outdir+out_filename_tavg+"_vhf"+out_filetype)
plt.close()

levavg_sum_vhf = vertical_average(ice_clip_tend_vhf+liq_clip_tend_vhf, rho, invers_dzt)
levavg_liq_vhf = vertical_average(liq_clip_tend_vhf, rho, invers_dzt)
levavg_ice_vhf = vertical_average(ice_clip_tend_vhf, rho, invers_dzt)

print("---- VertHF ----")
print("Sum: "+str(np.sum(levavg_sum_vhf)))
print("Liq: "+str(np.sum(levavg_liq_vhf)))
print("Ice: "+str(np.sum(levavg_ice_vhf)))

title = 'Sum of Vapor/Ice/Liquid mixing ratio tendencies due to vertical hole filling'
label = 'Mixing ratios (Sum)'
out_filepath = outdir+out_filename_zavg+"_sum_vhf"+out_filetype
plot_t_profile(title, label, levavg_sum_vhf, out_filepath)

title = 'Cloud liquid mixing ratio tendencies due to vertical hole filling'
label = 'Cloud liquid mixing ratio'
out_filepath = outdir+out_filename_zavg+"_liq_vhf"+out_filetype
plot_t_profile(title, label, levavg_liq_vhf, out_filepath)

title = 'Cloud ice mixing ratio tendencies due to vertical hole filling'
label = 'Cloud ice mixing ratio'
out_filepath = outdir+out_filename_zavg+"_ice_vhf"+out_filetype
plot_t_profile(title, label, levavg_ice_vhf, out_filepath)


# Plot profiles of clipping tendencies due to water vapor hole filling
ice_clip_tend_whf = get_2d_profile(nc,'INEGCLPTEND_WHF')
ice_clip_tend_whf_tavg = np.average(ice_clip_tend_whf, axis=0)

liq_clip_tend_whf = get_2d_profile(nc,'LNEGCLPTEND_WHF')
liq_clip_tend_whf_tavg = np.average(liq_clip_tend_whf, axis=0)

fig = plt.figure()
fig.text(.5,.95,title,horizontalalignment='center',)

ax1 = plt.subplot2grid((2, 2), (0, 0))
ax2 = plt.subplot2grid((2, 2), (0, 1))
ax3 = plt.subplot2grid((2, 2), (1, 0))

plot_z_profile(ax1,'INEGCLPTEND_WHF', ice_clip_tend_whf_tavg,lev)
plot_z_profile(ax2,'LNEGCLPTEND_WHF', liq_clip_tend_whf_tavg,lev)
plot_z_profile(ax3,'Ice+Liq WHF', ice_clip_tend_whf_tavg+liq_clip_tend_whf_tavg,lev)

fig.savefig(outdir+out_filename_tavg+"_whf"+out_filetype)
plt.close()

levavg_sum_whf = vertical_average(ice_clip_tend_whf+liq_clip_tend_whf, rho, invers_dzt)
levavg_liq_whf = vertical_average(liq_clip_tend_whf, rho, invers_dzt)
levavg_ice_whf = vertical_average(ice_clip_tend_whf, rho, invers_dzt)

print("---- WVHF ----")
print("Sum: "+str(np.sum(levavg_sum_whf)))
print("Liq: "+str(np.sum(levavg_liq_whf)))
print("Ice: "+str(np.sum(levavg_ice_whf)))

title = 'Sum of Vapor/Ice/Liquid mixing ratio tendencies due to water vapor hole filling'
label = 'Mixing ratios (Sum)'
out_filepath = outdir+out_filename_zavg+"_sum_whf"+out_filetype
plot_t_profile(title, label, levavg_sum_whf, out_filepath)

title = 'Cloud liquid mixing ratio tendencies due to water vapor hole filling'
label = 'Cloud liquid mixing ratio'
out_filepath = outdir+out_filename_zavg+"_liq_whf"+out_filetype
plot_t_profile(title, label, levavg_liq_whf, out_filepath)

title = 'Cloud ice mixing ratio tendencies due to water vapor hole filling'
label = 'Cloud ice mixing ratio'
out_filepath = outdir+out_filename_zavg+"_ice_whf"+out_filetype
plot_t_profile(title, label, levavg_ice_whf, out_filepath)

# Plot profiles of clipping tendencies due to clipping
ice_clip_tend_clp = get_2d_profile(nc,'INEGCLPTEND_CLP')
ice_clip_tend_clp_tavg = np.average(ice_clip_tend_clp, axis=0)

liq_clip_tend_clp = get_2d_profile(nc,'LNEGCLPTEND_CLP')
liq_clip_tend_clp_tavg = np.average(liq_clip_tend_clp, axis=0)

fig = plt.figure()
fig.text(.5,.95,title,horizontalalignment='center',)

ax1 = plt.subplot2grid((2, 2), (0, 0))
ax2 = plt.subplot2grid((2, 2), (0, 1))
ax3 = plt.subplot2grid((2, 2), (1, 0))

plot_z_profile(ax1,'INEGCLPTEND_CLP', ice_clip_tend_clp_tavg,lev)
plot_z_profile(ax2,'LNEGCLPTEND_CLP', liq_clip_tend_clp_tavg,lev)
plot_z_profile(ax3,'Ice+Liq CLP', ice_clip_tend_clp_tavg+liq_clip_tend_clp_tavg,lev)

fig.savefig(outdir+out_filename_tavg+"_clp"+out_filetype)
plt.close()

levavg_sum_clp = vertical_average(ice_clip_tend_clp+liq_clip_tend_clp, rho, invers_dzt)
levavg_liq_clp = vertical_average(liq_clip_tend_clp, rho, invers_dzt)
levavg_ice_clp = vertical_average(ice_clip_tend_clp, rho, invers_dzt)

print("---- CLIP ----")
print("Sum: "+str(np.sum(levavg_sum_clp)))
print("Liq: "+str(np.sum(levavg_liq_clp)))
print("Ice: "+str(np.sum(levavg_ice_clp)))

title = 'Sum of Vapor/Ice/Liquid mixing ratio tendencies due to clipping'
label = 'Mixing ratios (Sum)'
out_filepath = outdir+out_filename_zavg+"_sum_clp"+out_filetype
plot_t_profile(title, label, levavg_sum_clp, out_filepath)

title = 'Cloud liquid mixing ratio tendencies due to clipping'
label = 'Cloud liquid mixing ratio'
out_filepath = outdir+out_filename_zavg+"_liq_clp"+out_filetype
plot_t_profile(title, label, levavg_liq_clp, out_filepath)

title = 'Cloud ice mixing ratio tendencies due to clipping'
label = 'Cloud ice mixing ratio'
out_filepath = outdir+out_filename_zavg+"_ice_clp"+out_filetype
plot_t_profile(title, label, levavg_ice_clp, out_filepath)

# Plot profiles of clipping tendencies due to rest
ice_clip_tend_rst = get_2d_profile(nc,'INEGCLPTEND_RST')
ice_clip_tend_rst_tavg = np.average(ice_clip_tend_rst, axis=0)

liq_clip_tend_rst = get_2d_profile(nc,'LNEGCLPTEND_RST')
liq_clip_tend_rst_tavg = np.average(liq_clip_tend_rst, axis=0)

fig = plt.figure()
fig.text(.5,.95,title,horizontalalignment='center',)

ax1 = plt.subplot2grid((2, 2), (0, 0))
ax2 = plt.subplot2grid((2, 2), (0, 1))
ax3 = plt.subplot2grid((2, 2), (1, 0))

plot_z_profile(ax1,'INEGCLPTEND_RST', ice_clip_tend_rst_tavg,lev)
plot_z_profile(ax2,'LNEGCLPTEND_RST', liq_clip_tend_rst_tavg,lev)
plot_z_profile(ax3,'Ice+Liq RST', ice_clip_tend_rst_tavg+liq_clip_tend_rst_tavg,lev)

fig.savefig(outdir+out_filename_tavg+"_rst"+out_filetype)
plt.close()

levavg_sum_rst = vertical_average(ice_clip_tend_rst+liq_clip_tend_rst, rho, invers_dzt)
levavg_liq_rst = vertical_average(liq_clip_tend_rst, rho, invers_dzt)
levavg_ice_rst = vertical_average(ice_clip_tend_rst, rho, invers_dzt)

print("---- REST ----")
print("Sum: "+str(np.sum(levavg_sum_rst)))
print("Liq: "+str(np.sum(levavg_liq_rst)))
print("Ice: "+str(np.sum(levavg_ice_rst)))

title = 'Sum of Vapor/Ice/Liquid mixing ratio tendencies due to clipping'
label = 'Mixing ratios (Sum)'
out_filepath = outdir+out_filename_zavg+"_sum_rst"+out_filetype
plot_t_profile(title, label, levavg_sum_rst, out_filepath)

title = 'Cloud liquid mixing ratio tendencies due to clipping'
label = 'Cloud liquid mixing ratio'
out_filepath = outdir+out_filename_zavg+"_liq_rst"+out_filetype
plot_t_profile(title, label, levavg_liq_rst, out_filepath)

title = 'Cloud ice mixing ratio tendencies due to clipping'
label = 'Cloud ice mixing ratio'
out_filepath = outdir+out_filename_zavg+"_ice_rst"+out_filetype
plot_t_profile(title, label, levavg_ice_rst, out_filepath)

# Plot profiles of clipping tendencies due to qneg3 in physics_update
ice_clip_tend_nohf_pu = get_2d_profile(nc,'INEGCLPTEND_NOHF_PU')
ice_clip_tend_nohf_pu_tavg = np.average(ice_clip_tend_nohf_pu, axis=0)

liq_clip_tend_nohf_pu = get_2d_profile(nc,'LNEGCLPTEND_NOHF_PU')
liq_clip_tend_nohf_pu_tavg = np.average(liq_clip_tend_nohf_pu, axis=0)

vap_clip_tend_nohf_pu = get_2d_profile(nc,'VNEGCLPTEND_NOHF_PU')
vap_clip_tend_nohf_pu_tavg = np.average(liq_clip_tend_nohf_pu, axis=0)

fig = plt.figure()
fig.text(.5,.95,title,horizontalalignment='center',)

ax1 = plt.subplot2grid((2, 2), (0, 0))
ax2 = plt.subplot2grid((2, 2), (0, 1))
ax3 = plt.subplot2grid((2, 2), (1, 0))
ax4 = plt.subplot2grid((2, 2), (1, 1))

plot_z_profile(ax1,'INEGCLPTEND_NOHF_PU', ice_clip_tend_nohf_pu_tavg,lev)
plot_z_profile(ax2,'LNEGCLPTEND_NOHF_PU', liq_clip_tend_nohf_pu_tavg,lev)
plot_z_profile(ax3,'VNEGCLPTEND_NOHF_PU', vap_clip_tend_nohf_pu_tavg,lev)
plot_z_profile(ax4,'Ice+Liq+Vap NOHF PU', ice_clip_tend_nohf_pu_tavg+liq_clip_tend_nohf_pu_tavg+vap_clip_tend_nohf_pu_tavg,lev)

fig.savefig(outdir+out_filename_tavg+"_nohf_pu"+out_filetype)
plt.close()

levavg_sum_nohf_pu = vertical_average(ice_clip_tend_nohf_pu+liq_clip_tend_nohf_pu+vap_clip_tend_nohf_pu, rho, invers_dzt)
levavg_liq_nohf_pu = vertical_average(liq_clip_tend_nohf_pu, rho, invers_dzt)
levavg_ice_nohf_pu = vertical_average(ice_clip_tend_nohf_pu, rho, invers_dzt)
levavg_vap_nohf_pu = vertical_average(vap_clip_tend_nohf_pu, rho, invers_dzt)

print("---- Phys update no hf ----")
print("Sum: "+str(np.sum(levavg_sum_nohf_pu)))
print("Liq: "+str(np.sum(levavg_liq_nohf_pu)))
print("Ice: "+str(np.sum(levavg_ice_nohf_pu)))

title = 'Sum of Vapor/Ice/Liquid mixing ratio tendencies due to clipping'
label = 'Mixing ratios (Sum)'
out_filepath = outdir+out_filename_zavg+"_sum_nohf_pu"+out_filetype
plot_t_profile(title, label, levavg_sum_nohf_pu, out_filepath)

title = 'Cloud liquid mixing ratio tendencies due to clipping'
label = 'Cloud liquid mixing ratio'
out_filepath = outdir+out_filename_zavg+"_liq_nohf_pu"+out_filetype
plot_t_profile(title, label, levavg_liq_nohf_pu, out_filepath)

title = 'Cloud ice mixing ratio tendencies due to clipping'
label = 'Cloud ice mixing ratio'
out_filepath = outdir+out_filename_zavg+"_ice_nohf_pu"+out_filetype
plot_t_profile(title, label, levavg_ice_nohf_pu, out_filepath)

title = 'Water vapor mixing ratio tendencies due to clipping'
label = 'Water vapor mixing ratio'
out_filepath = outdir+out_filename_zavg+"_vap_nohf_pu"+out_filetype
plot_t_profile(title, label, levavg_vap_nohf_pu, out_filepath)

# Plot profiles of clipping tendencies due to qneg3 (TPHYSBCB)
ice_clip_tend_nohf_b = get_2d_profile(nc,'INEGCLPTEND_NOHF_B')
ice_clip_tend_nohf_b_tavg = np.average(ice_clip_tend_nohf_b, axis=0)

liq_clip_tend_nohf_b = get_2d_profile(nc,'LNEGCLPTEND_NOHF_B')
liq_clip_tend_nohf_b_tavg = np.average(liq_clip_tend_nohf_b, axis=0)

vap_clip_tend_nohf_b = get_2d_profile(nc,'VNEGCLPTEND_NOHF_B')
vap_clip_tend_nohf_b_tavg = np.average(liq_clip_tend_nohf_b, axis=0)

fig = plt.figure()
fig.text(.5,.95,title,horizontalalignment='center',)

ax1 = plt.subplot2grid((2, 2), (0, 0))
ax2 = plt.subplot2grid((2, 2), (0, 1))
ax3 = plt.subplot2grid((2, 2), (1, 0))
ax4 = plt.subplot2grid((2, 2), (1, 1))

plot_z_profile(ax1,'INEGCLPTEND_NOHF_B', ice_clip_tend_nohf_b_tavg,lev)
plot_z_profile(ax2,'LNEGCLPTEND_NOHF_B', liq_clip_tend_nohf_b_tavg,lev)
plot_z_profile(ax3,'VNEGCLPTEND_NOHF_B', vap_clip_tend_nohf_b_tavg,lev)
plot_z_profile(ax4,'Ice+Liq+Vap NOHF TPHYSBCB', ice_clip_tend_nohf_b_tavg+liq_clip_tend_nohf_b_tavg+vap_clip_tend_nohf_b_tavg,lev)

fig.savefig(outdir+out_filename_tavg+"_nohf_b"+out_filetype)
plt.close()

levavg_sum_nohf_b = vertical_average(ice_clip_tend_nohf_b+liq_clip_tend_nohf_b+vap_clip_tend_nohf_b, rho, invers_dzt)
levavg_liq_nohf_b = vertical_average(liq_clip_tend_nohf_b, rho, invers_dzt)
levavg_ice_nohf_b = vertical_average(ice_clip_tend_nohf_b, rho, invers_dzt)
levavg_vap_nohf_b = vertical_average(vap_clip_tend_nohf_b, rho, invers_dzt)

print("---- TPHYSBCB no hf ----")
print("Sum: "+str(np.sum(levavg_sum_nohf_b)))
print("Liq: "+str(np.sum(levavg_liq_nohf_b)))
print("Ice: "+str(np.sum(levavg_ice_nohf_b)))

title = 'Sum of Vapor/Ice/Liquid mixing ratio tendencies due to clipping'
label = 'Mixing ratios (Sum)'
out_filepath = outdir+out_filename_zavg+"_sum_nohf_b"+out_filetype
plot_t_profile(title, label, levavg_sum_nohf_b, out_filepath)

title = 'Cloud liquid mixing ratio tendencies due to clipping'
label = 'Cloud liquid mixing ratio'
out_filepath = outdir+out_filename_zavg+"_liq_nohf_b"+out_filetype
plot_t_profile(title, label, levavg_liq_nohf_b, out_filepath)

title = 'Cloud ice mixing ratio tendencies due to clipping'
label = 'Cloud ice mixing ratio'
out_filepath = outdir+out_filename_zavg+"_ice_nohf_b"+out_filetype
plot_t_profile(title, label, levavg_ice_nohf_b, out_filepath)

title = 'Water vapor mixing ratio tendencies due to clipping'
label = 'Water vapor mixing ratio'
out_filepath = outdir+out_filename_zavg+"_vap_nohf_b"+out_filetype
plot_t_profile(title, label, levavg_vap_nohf_b, out_filepath)

# Plot profiles of clipping tendencies due to qneg3 (TPHYSBCC)
ice_clip_tend_nohf_c = get_2d_profile(nc,'INEGCLPTEND_NOHF_C')
ice_clip_tend_nohf_c_tavg = np.average(ice_clip_tend_nohf_c, axis=0)

liq_clip_tend_nohf_c = get_2d_profile(nc,'LNEGCLPTEND_NOHF_C')
liq_clip_tend_nohf_c_tavg = np.average(liq_clip_tend_nohf_c, axis=0)

vap_clip_tend_nohf_c = get_2d_profile(nc,'VNEGCLPTEND_NOHF_C')
vap_clip_tend_nohf_c_tavg = np.average(liq_clip_tend_nohf_c, axis=0)

fig = plt.figure()
fig.text(.5,.95,title,horizontalalignment='center',)

ax1 = plt.subplot2grid((2, 2), (0, 0))
ax2 = plt.subplot2grid((2, 2), (0, 1))
ax3 = plt.subplot2grid((2, 2), (1, 0))
ax4 = plt.subplot2grid((2, 2), (1, 1))

plot_z_profile(ax1,'INEGCLPTEND_NOHF_C', ice_clip_tend_nohf_c_tavg,lev)
plot_z_profile(ax2,'LNEGCLPTEND_NOHF_C', liq_clip_tend_nohf_c_tavg,lev)
plot_z_profile(ax3,'VNEGCLPTEND_NOHF_C', vap_clip_tend_nohf_c_tavg,lev)
plot_z_profile(ax4,'Ice+Liq+Vap NOHF TPHYSBCC', ice_clip_tend_nohf_c_tavg+liq_clip_tend_nohf_c_tavg+vap_clip_tend_nohf_c_tavg,lev)

fig.savefig(outdir+out_filename_tavg+"_nohf_c"+out_filetype)
plt.close()

levavg_sum_nohf_c = vertical_average(ice_clip_tend_nohf_c+liq_clip_tend_nohf_c+vap_clip_tend_nohf_c, rho, invers_dzt)
levavg_liq_nohf_c = vertical_average(liq_clip_tend_nohf_c, rho, invers_dzt)
levavg_ice_nohf_c = vertical_average(ice_clip_tend_nohf_c, rho, invers_dzt)
levavg_vap_nohf_c = vertical_average(vap_clip_tend_nohf_c, rho, invers_dzt)

print("---- TPHYSBCB no hf ----")
print("Sum: "+str(np.sum(levavg_sum_nohf_c)))
print("Liq: "+str(np.sum(levavg_liq_nohf_c)))
print("Ice: "+str(np.sum(levavg_ice_nohf_c)))

title = 'Sum of Vapor/Ice/Liquid mixing ratio tendencies due to clipping'
label = 'Mixing ratios (Sum)'
out_filepath = outdir+out_filename_zavg+"_sum_nohf_c"+out_filetype
plot_t_profile(title, label, levavg_sum_nohf_c, out_filepath)

title = 'Cloud liquid mixing ratio tendencies due to clipping'
label = 'Cloud liquid mixing ratio'
out_filepath = outdir+out_filename_zavg+"_liq_nohf_c"+out_filetype
plot_t_profile(title, label, levavg_liq_nohf_c, out_filepath)

title = 'Cloud ice mixing ratio tendencies due to clipping'
label = 'Cloud ice mixing ratio'
out_filepath = outdir+out_filename_zavg+"_ice_nohf_c"+out_filetype
plot_t_profile(title, label, levavg_ice_nohf_c, out_filepath)

title = 'Water vapor mixing ratio tendencies due to clipping'
label = 'Water vapor mixing ratio'
out_filepath = outdir+out_filename_zavg+"_vap_nohf_c"+out_filetype
plot_t_profile(title, label, levavg_vap_nohf_c, out_filepath)


print('---- TEST (should be all zero) ----')
print(np.sum(levavg_sum_vhf)+np.sum(levavg_sum_whf)+np.sum(levavg_sum_clp)+np.sum(levavg_sum_rst)-np.sum(levavg_sum)+np.sum(levavg_vap))
print(np.sum(levavg_liq_vhf)+np.sum(levavg_liq_whf)+np.sum(levavg_liq_clp)+np.sum(levavg_liq_rst)-np.sum(levavg_liq))
print(np.sum(levavg_ice_vhf)+np.sum(levavg_ice_whf)+np.sum(levavg_ice_clp)+np.sum(levavg_ice_rst)-np.sum(levavg_ice))

# Those profiles should match

fig = plt.figure()
fig.text(.5,.95,title,horizontalalignment='center',)
plt.plot(range(2001), levavg_vap)
plt.plot(range(2001), -levavg_sum_whf, 'ro')
fig.savefig(outdir+out_filename_zavg+"_test"+out_filetype)


fig = plt.figure()
fig.text(.5,.95,title,horizontalalignment='center',)
plt.plot(range(2001), levavg_sum-levavg_vap)
plt.plot(range(2001), levavg_sum_vhf+levavg_sum_whf+levavg_sum_clp+levavg_sum_rst, 'ro')
fig.savefig(outdir+out_filename_zavg+"_test2"+out_filetype)


fig = plt.figure()
fig.text(.5,.95,title,horizontalalignment='center',)

liq_clip_tend_vhf_tavg = np.average(liq_clip_tend_vhf, axis=0)
liq_clip_tend_whf_tavg = np.average(liq_clip_tend_whf, axis=0)
liq_clip_tend_clp_tavg = np.average(liq_clip_tend_clp, axis=0)
liq_clip_tend_rst_tavg = np.average(liq_clip_tend_rst, axis=0)

plt.plot(liq_clip_tend_tavg, range(nlev))
plt.plot(liq_clip_tend_vhf_tavg+liq_clip_tend_whf_tavg+liq_clip_tend_clp_tavg+liq_clip_tend_rst_tavg, range(nlev),'r--')

fig.savefig(outdir+out_filename_tavg+"_compare_liq"+out_filetype)

fig = plt.figure()
fig.text(.5,.95,title,horizontalalignment='center',)

ice_clip_tend_vhf_tavg = np.average(ice_clip_tend_vhf, axis=0)
ice_clip_tend_whf_tavg = np.average(ice_clip_tend_whf, axis=0)
ice_clip_tend_clp_tavg = np.average(ice_clip_tend_clp, axis=0)
ice_clip_tend_rst_tavg = np.average(ice_clip_tend_rst, axis=0)

plt.plot(ice_clip_tend_tavg, range(nlev))
plt.plot(ice_clip_tend_vhf_tavg+ice_clip_tend_whf_tavg+ice_clip_tend_clp_tavg+ice_clip_tend_rst_tavg, range(nlev),'r--')


# Overplot profiles for comparison

fig.savefig(outdir+out_filename_tavg+"_compare_ice"+out_filetype)
plt.close()

fig = plt.figure()
title = 'Liquid time-averaged profiles'
fig.text(.5,.95,title,horizontalalignment='center',)

plt_tot, = plt.plot(liq_clip_tend_tavg, zm[1,0:nlev], 'b')
plt_vhf, = plt.plot(liq_clip_tend_vhf_tavg, zm[1,0:nlev],'r')
plt_whf, = plt.plot(liq_clip_tend_whf_tavg, zm[1,0:nlev], 'g')
plt_clp, = plt.plot(liq_clip_tend_clp_tavg, zm[1,0:nlev], 'm')
plt_rst, = plt.plot(liq_clip_tend_rst_tavg, zm[1,0:nlev], 'c')
plt.legend([plt_tot, plt_vhf, plt_whf, plt_clp, plt_rst], ['Total', 'Vert. HF', 'WV HF', 'Clipping', 'Rest'], prop={'size':'small'})

fig.savefig(outdir+out_filename_tavg+"_liq_profile_cmp"+out_filetype)
plt.close()

fig = plt.figure()
title = 'Ice time-averaged profiles'
fig.text(.5,.95,title,horizontalalignment='center',)

plt_tot, = plt.plot(ice_clip_tend_tavg, zm[1,0:nlev], 'b')
plt_vhf, = plt.plot(ice_clip_tend_vhf_tavg, zm[1,0:nlev],'r')
plt_whf, = plt.plot(ice_clip_tend_whf_tavg, zm[1,0:nlev], 'g')
plt_clp, = plt.plot(ice_clip_tend_clp_tavg, zm[1,0:nlev], 'm')
plt_rst, = plt.plot(ice_clip_tend_rst_tavg, zm[1,0:nlev], 'c')
plt.legend([plt_tot, plt_vhf, plt_whf, plt_clp, plt_rst], ['Total', 'Vert. HF', 'WV HF', 'Clipping', 'Rest'], prop={'size':'small'})

fig.savefig(outdir+out_filename_tavg+"_ice_profile_cmp"+out_filetype)
plt.close()

fig = plt.figure()
title = 'Liquid height-averaged timelines'
fig.text(.5,.95,title,horizontalalignment='center',)

plt_tot, = plt.plot(range(ntimes), levavg_liq, 'b', linewidth=2.0)
plt_vhf, = plt.plot(range(ntimes), levavg_liq_vhf,'r-.', linewidth=2.0)
plt_whf, = plt.plot(range(ntimes), levavg_liq_whf, 'g:', linewidth=2.0)
plt_clp, = plt.plot(range(ntimes), levavg_liq_clp, 'm--', linewidth=2.0)
plt_rst, = plt.plot(range(ntimes), levavg_liq_rst, 'c--', linewidth=2.0)
plt.legend([plt_tot, plt_vhf, plt_whf, plt_clp, plt_rst], ['Total', 'Vert. HF', 'WV HF', 'Clipping', 'Rest'], loc=2, prop={'size':'small'})

fig.savefig(outdir+out_filename_tavg+"_liq_timeln_cmp"+out_filetype)
plt.close()

fig = plt.figure()
title = 'Ice height-averaged timelines'
fig.text(.5,.95,title,horizontalalignment='center',)

plt_tot, = plt.plot(range(ntimes), levavg_ice, 'b', linewidth=2.0)
plt_vhf, = plt.plot(range(ntimes), levavg_ice_vhf,'r-.', linewidth=2.0)
plt_whf, = plt.plot(range(ntimes), levavg_ice_whf, 'g:', linewidth=2.0)
plt_clp, = plt.plot(range(ntimes), levavg_ice_clp, 'm--', linewidth=2.0)
plt_rst, = plt.plot(range(ntimes), levavg_ice_rst, 'c--', linewidth=2.0)
plt.legend([plt_tot, plt_vhf, plt_whf, plt_clp, plt_rst], ['Total', 'Vert. HF', 'WV HF', 'Clipping', 'Rest'], loc=2, prop={'size':'small'})

fig.savefig(outdir+out_filename_tavg+"_ice_timeln_cmp"+out_filetype)
plt.close()

if ( l_show_plots ):
   plt.show()

