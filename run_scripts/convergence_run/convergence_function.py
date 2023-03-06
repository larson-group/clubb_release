import argparse
import numpy as np
import os
import sys
import shutil

def modify_ic_profile(clubb_dir, case_dir, grid_dir, case, dz, refine): 

  print('construct initial condition file for ' + case)

  if (case != "dycoms2_rf02"): 
    # read the initial sounding profile and increase height levels 
    # note: the data in new levels are all set to missing values as 
    #       we do not attempt to change the original sounding profiles 
    #       the purpose is to have denser height levels so that the cubic 
    #       spline interpolation does not generate strange profiles 
    input_file_name = case + '_sounding.in'
    input_file  = open(os.path.join(case_dir, input_file_name), 'r')
    input_lines = input_file.readlines()
    input_file.close()
    line_collections = []
    line_data = []
    for line in input_lines:
      if (line.strip().startswith('!')):
        line_collections.append(line)
      elif (line.strip().startswith('z[m]')):
        line_collections.append(line)
        refline = line
      else: 
        line_data.append(line)
    zdat = []
    znew = []
    zmin = 0 
    zmax = 0
    for line in line_data: 
      zdat.append(float(line.split()[0]))
      if (float(line.split()[0]) < zmin): 
        zmin = float(line.split()[0])
      if (float(line.split()[0]) > zmax):
        zmax = float(line.split()[0])
    for val in np.arange(zmin,zmax,100):
      znew.append(val)
    zlev = znew + list(set(zdat)- set(znew))
    zlev.sort()
    nlev = len(zlev)
    for i in np.arange(len(zlev)):
      line = refline
      for item in refline.split():
        if (item == "z[m]"):
          line = line.replace(item, str("{:.1f}".format(zlev[i])))
        else: 
          line = line.replace(item, '-999.9')
      for j in np.arange(len(zdat)): 
        if (zdat[j] == zlev[i]):
          line = line_data[j]
      line_collections.append(line)

    out_file_name = case + '_sounding.in'
    output_file = open(os.path.join(case_dir, out_file_name), 'w')
    for line in line_collections:
      output_file.write(line)
    output_file.close()

  else: # case == "dycoms2_rf02"
    # read namelist file to obtain model height in case setup 
    model_file_name = case + '_nd_model.in'
    input_file = open(os.path.join(case_dir, model_file_name), 'r')
    input_lines = input_file.readlines()
    input_file.close()
    for line in input_lines:
      if ('zm_top' in line):
        height =  float(line.split()[2])
        htop   = height + 150.0 
      if ('zm_init' in line):
        hgtlow =  float(line.split()[2])

    # if refinement specified, create grid file and set appropriate file name
    if (refine != -9999):
      #for stretched grid 
      gridtmp = np.loadtxt(os.path.join(grid_dir,'deep_convection_128lev_27km_zt_grid.grd'))
      grid128 = 0.5*(gridtmp[0:len(gridtmp)-2] + gridtmp[1:len(gridtmp)-1])
      ind = np.where(grid128 > htop)[0][0]
      grid128 = grid128[:ind]
      grid128[0] = 0.0
      nlev = (len(grid128)-1)*2**refine + 1
      grid = np.empty((nlev))
      coarse_ind = np.arange(len(grid128))
      refine_ind = np.arange(0,nlev,2**refine)
      grid[refine_ind] = grid128[coarse_ind]
      for level in range(refine):
        coarse_ind = np.arange(0, nlev, 2**(refine-level))
        refine_ind = coarse_ind[:-1] + 2**(refine-level-1)
        grid[refine_ind] = 0.5*(grid[coarse_ind[1:]] + grid[coarse_ind[:-1]])
    elif (dz != -9999): 
      # for grid with fixed spacing  
      grid = np.arange(hgtlow, htop, dz, dtype=float)
      nlev = len(grid)
    else:  
      sys.exit('No refinement parameter specified, no need to construct the profile')
   
    #construct initial condition profile (Wyant, et al. 2007, eq 1--4) 
    #rtm and thlm profile 
    zthres   = 795.0
    g_per_kg = 1000.0
    thlm = 295.0 + (np.sign(grid - zthres)) * (np.abs(grid - zthres)) ** (1.0 / 3.0)
    rtm = (5.0 - 3.0 * (1.0 - np.exp((zthres - grid)/500.0))) / g_per_kg
    for i in np.arange(nlev):
      if (grid[i] < zthres):
        rtm[i] = 9.45/g_per_kg 
        thlm[i] = 288.3
    #um and vm profile 
    um = 3.0 + (4.3*grid/1000.0)
    vm = -9.0 + (5.6*grid/1000.0)
    ug = um
    vg = vm
    #wm profile 
    z0 = [0, height, htop]
    w0 = [0.,-0.006, -0.006]
    wm = np.interp(grid, z0, w0)
    #apply smoothing on the profiles of thlm and rtm 
    #smoothed heaviside function with dzt = 30 m 
    dzt   = 30
    thlm2 = 295.0 + (dzt)**(1.0/3.0)
    thlm1 = 288.3
    rtm2  = (5.0 - 3.0 * (1.0 - np.exp(-dzt/500.0))) / g_per_kg
    rtm1  = 9.45/g_per_kg
    frac  = (grid - zthres) / dzt
    for i in np.arange(nlev):
     if( frac[i] >= -1.0 and frac[i] <= 1.0 ):
        H0 = 0.5*(1.0+frac[i] + (1.0/np.pi)*np.sin(np.pi*frac[i]))
        thlm[i] = thlm1 + H0*(thlm2 - thlm1)
        rtm[i]  = rtm1  + H0*(rtm2 - rtm1)
  
    #write out the sounding profile 
    input_file_name = case + '_sounding.in'
    input_file = open(os.path.join(case_dir, input_file_name), 'r')
    input_lines = input_file.readlines()
    input_file.close()
    line_collections = []
    for line in input_lines:
      if (line.strip().startswith('!')):
        line_collections.append(line)
      if (line.strip().startswith('z[m]')):
        line_collections.append(line)
        refline = line
    
    for i in np.arange(nlev):
      line = refline 
      for item in refline.split():
        if(item == "z[m]"): 
           line = line.replace(item, str("{:.3f}".format(grid[i])))
        if(item == "thlm[K]"):
           line = line.replace(item, str("{:.5f}".format(thlm[i])))
        if(item == "rt[kg\\kg]"):
           line = line.replace(item, str("{:.5f}".format(rtm[i])))
        if(item == "u[m\\s]"):
           line = line.replace(item, str("{:.5f}".format(um[i])))
        if(item == "v[m\\s]"):
           line = line.replace(item, str("{:.5f}".format(vm[i])))
        if(item == "w[m\\s]"):
           line = line.replace(item, str("{:.6f}".format(wm[i])))     
        if(item == "ug[m\\s]"):
           line = line.replace(item, str("{:.5f}".format(ug[i])))
        if(item == "vg[m\\s]"):
           line = line.replace(item, str("{:.5f}".format(vg[i])))
      line_collections.append(line)
     
    out_file_name = case + '_sounding.in' 
    output_file = open(os.path.join(case_dir, out_file_name), 'w')
    for line in line_collections:
      output_file.write(line)
    output_file.close()
  
  return 
