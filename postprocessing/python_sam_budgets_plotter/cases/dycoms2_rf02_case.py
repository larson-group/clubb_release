"""
-------------------------------------------------------------------------------
G E N E R A L   I N F O R M A T I O N
-------------------------------------------------------------------------------
This file contains the case specific input needed for plotgen.py.
to plot height profiles and budgets for the DYCOMS2_RF02 case
"""

#-------------------------------------------------------------------------------
#   G E N E R A L   D E F A U L T   S E T T I N G S
#-------------------------------------------------------------------------------

#case = 'DYCOMS2_RF02'
case = 'DYCOMS_RF02'
full_name = 'DYCOMS-II RF02'
case_folder = '/home/sdomke/workspace/clubb/sam_clubb/{case}'.format(case=case)
enabled = True 
#type = 'budget' # not used in plotgen
nx = 128
ny = 128
nz = 96
dxy = 50                # [m]
dz = '5-80'             # [m]
dt = .5                 # [s]
startTime = 121.0       # [minutes]
endTime = 240.0         # [minutes]
startHeight = 20.0      # [m]
endHeight = 1200.0      # [m]
time_3d = 43200.0       # [s]

# run entry for header in html file
run = '{case}_{nx}x{ny}x{nz}_{dxy}m_{dz}m_{dt}s'.format(case=case, nx=nx, ny=ny, nz=nz, dxy=dxy, dz=dz, dt=dt)
# jpg file names
plot_case_name = '{case}_{dx}x{dx}_{{type}}_{{date}}_{{plot}}'.format(case=case.lower(), dx=nx)
## use absolute paths or relative paths originating from the directory containing plotgen.py
# directory for output files
out_dir = '/home/sdomke/workspace/plotgen_out/{case}_{{date}}/'.format(case=case.lower())
# pdf output name
out_pdf = '{case}_{dx}x{dx}_{{type}}_{{date}}.pdf'.format(case=case.lower(),dx=nx)

## input .nc file
## SAM
# nc file generated from .stat output
#sam_file = '/home/sdomke/workspace/clubb/avi_out/DYCOMS_RF02_128x128x96_dr_nosed_190206.nc'
sam_file = '/home/sdomke/workspace/clubb/avi_out/grid_change/DYCOMS_RF02_128x128x96_dr_nosed.nc'
#sam_file = '/home/sdomke/workspace/clubb/sam_out/outstat/DYCOMS_RF02_64x64x96_7200_dr_nosed_181205.nc'
# nc file generated from .bin3D output
sam_3d_file = '/home/sdomke/workspace/clubb/avi_out/out3d/DYCOMS_RF02_128x128x96_dr_nosed_64_0000043200.nc'

# nc files for publishing runs with bigger horizontal grid (256x256):
out_dir = '/home/sdomke/workspace/plotgen_out/publishing_runs/{case}_{{date}}/'.format(case=case.lower())
sam_file = '/home/sdomke/workspace/clubb/avi_out/publishing_runs/DYCOMS_RF02_256x256x96_dr_nosed.nc'
sam_3d_file = '/home/sdomke/workspace/clubb/avi_out/publishing_runs/3d/DYCOMS_RF02_256x256x96_dr_nosed_256_0000043200.nc'

## CLUBB
clubb_zm_file = '/home/sdomke/workspace/clubb/clubb_out/dycoms2_rf02_do_zm.nc'
clubb_zt_file = '/home/sdomke/workspace/clubb/clubb_out/dycoms2_rf02_do_zt.nc'
## old CLUBB
old_clubb_zm_file = '/home/sdomke/workspace/clubb/clubb_out/dycoms2_rf02_do_zm_old.nc'
old_clubb_zt_file = '/home/sdomke/workspace/clubb/clubb_out/dycoms2_rf02_do_zt_old.nc'

## case setup files
sam_prm = case_folder+'/prm.les2'
sam_grd = case_folder+'/grd.les'

# header in html file
headerText = '{run} {{type}} Minutes {start}-{end}, {bottom}m-{top}m'.format(run=run, start=startTime, end=endTime, bottom=startHeight, top=endHeight)

#-------------------------------------------------------------------------------
#   G E N E R A L   P L O T   S E T T I N G S
#-------------------------------------------------------------------------------
lw = 5
color = 'nipy_spectral'
yLabel = 'Height [m]'
