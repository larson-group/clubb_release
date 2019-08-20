"""
-------------------------------------------------------------------------------
G E N E R A L   I N F O R M A T I O N
-------------------------------------------------------------------------------
This file contains the case specific input needed for plotgen.py.
to plot height profiles and budgets for the LBA case
"""

#-------------------------------------------------------------------------------
#   G E N E R A L   D E F A U L T   S E T T I N G S
#-------------------------------------------------------------------------------

case = 'LBA'
full_name = case
case_folder = '/home/sdomke/workspace/clubb/sam_clubb/{case}'.format(case=case)
enabled = True
#type = 'budget' # not used in plotgen
nx = 128
ny = 128
nz = 128
dxy = 1000            # [m]
dz = '25-250'         # [m]
dt = 6                # [s]
startTime = 301.0     # [minutes]
endTime = 360.0       # [minutes]
startHeight = 0.0     # [m]
endHeight = 15000.0   # [m]
time_3d = 3600.0      # dt

# run entry for header in html file
run = '{case}_{nx}x{ny}x{nz}_{dxy}m_{dz}m_{dt}s'.format(case=case, nx=nx, ny=ny, nz=nz, dxy=dxy, dz=dz, dt=dt)
# jpg file names
plot_case_name = '{case}_{dx}x{dx}_{{type}}_{{date}}_{{plot}}'.format(case=case.lower(), dx=nx)
## use absolute paths or relative paths originating from the directory containing plotgen.py
# directory for output files
out_dir = '../../output/plotgen_out/{case}_{{date}}/'.format(case=case.lower())
# pdf output name
out_pdf = '{case}_{dx}x{dx}_{{type}}_{{date}}.pdf'.format(case=case.lower(),dx=nx)

## input .nc file
## SAM
# nc file generated from .stat output
#sam_file = '/home/sdomke/workspace/clubb/avi_out/LBA_128kmx128kmx128_1km_Morrison_190205.nc'
sam_file = '../../output/LBA_128kmx128kmx128_1km_Morrison.nc'
# nc file generated from .bin3D output
sam_3d_file = '../../output/sam_3d/LBA_128kmx128kmx128_1km_Morrison_64_0000003600.nc'

## CLUBB (case not included)
clubb_zm_file = '../../output/prog_mom_flux_paper/new/lba_zm.nc'
clubb_zt_file = '../../output/prog_mom_flux_paper/new/lba_zt.nc'
## old CLUBB
old_clubb_zm_file = '../../output/prog_mom_flux_paper/old/lba_zm.nc'
old_clubb_zt_file = '../../output/prog_mom_flux_paper/old/lba_zt.nc'
## case setup files
sam_prm = case_folder+'/prm.les'
sam_grd = case_folder+'/grd.les'

# header in html file
headerText = '{run} {{type}} Minutes {start}-{end}, {bottom}m-{top}m'.format(run=run, start=startTime, end=endTime, bottom=startHeight, top=endHeight)

#-------------------------------------------------------------------------------
#   G E N E R A L   P L O T   S E T T I N G S
#-------------------------------------------------------------------------------
lw = 5
color = 'nipy_spectral'
yLabel = 'Height [m]'