"""
-------------------------------------------------------------------------------
G E N E R A L   I N F O R M A T I O N
-------------------------------------------------------------------------------
This file contains the case specific input needed for plotgen.py.
to plot height profiles and budgets for the DYCOMS2_RF01 case
"""

#-------------------------------------------------------------------------------
#   G E N E R A L   D E F A U L T   S E T T I N G S
#-------------------------------------------------------------------------------

#case = 'DYCOMS2_RF01'
case = 'DYCOMS_RF01'
full_name = 'DYCOMS-II RF01'
case_folder = '/home/sdomke/workspace/clubb/sam_clubb/{case}'.format(case=case)
enabled = True
#type = 'budget' # not used in plotgen
nx = 256
ny = 256
nz = 320
dxy = 35                # [m]
dz = 5                  # [m]
dt = .5                 # [s]
startTime = 121.0       # [minutes]
endTime = 240.0         # [minutes]
startHeight = 20.0      # [m]
endHeight = 1100.0      # [m]
time_3d = 28800.0       # [s]

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
#sam_file = '/home/sdomke/workspace/clubb/avi_out/DYCOMS_RF01_96x96x320_190207.nc'
sam_file = '/home/sdomke/workspace/clubb/avi_out/grid_change/DYCOMS_RF01_96x96x320.nc'
# nc file generated from .bin3D output
sam_3d_file = '/home/sdomke/workspace/clubb/avi_out/out3d/DYCOMS_RF01_96x96x320_96_0000028800.nc'

# nc files for publishing runs with bigger horizontal grid (256x256):
out_dir = '/home/sdomke/workspace/plotgen_out/publishing_runs/{case}_{{date}}/'.format(case=case.lower())
sam_file = '/home/sdomke/workspace/clubb/avi_out/publishing_runs/DYCOMS_RF01_256x256x320.nc'
sam_3d_file = '/home/sdomke/workspace/clubb/avi_out/publishing_runs/3d/DYCOMS_RF01_256x256x320_256_0000028800.nc'

## CLUBB
clubb_zm_file = '/home/sdomke/workspace/clubb/clubb_out/dycoms2_rf01_zm.nc'
clubb_zt_file = '/home/sdomke/workspace/clubb/clubb_out/dycoms2_rf01_zt.nc'
## old CLUBB
old_clubb_zm_file = '/home/sdomke/workspace/clubb/clubb_out/dycoms2_rf01_zm.nc'
old_clubb_zt_file = '/home/sdomke/workspace/clubb/clubb_out/dycoms2_rf01_zt.nc'

## case setup files
sam_prm = case_folder+'/prm'
sam_grd = case_folder+'/grd'

# header in html file
headerText = '{run} {{type}} Minutes {start}-{end}, {bottom}m-{top}m'.format(run=run, start=startTime, end=endTime, bottom=startHeight, top=endHeight)

#-------------------------------------------------------------------------------
#   G E N E R A L   P L O T   S E T T I N G S
#-------------------------------------------------------------------------------
lw = 5
color = 'nipy_spectral'
yLabel = 'Height [m]'
