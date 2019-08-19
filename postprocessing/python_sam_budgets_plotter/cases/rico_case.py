"""
-------------------------------------------------------------------------------
G E N E R A L   I N F O R M A T I O N
-------------------------------------------------------------------------------
This file contains the case specific input needed for plotgen.py.
to plot height profiles and budgets for the RICO case
"""

#-------------------------------------------------------------------------------
#   G E N E R A L   D E F A U L T   S E T T I N G S
#-------------------------------------------------------------------------------

case = 'RICO'
full_name = case
case_folder = '/home/sdomke/workspace/clubb/sam_clubb/{case}'.format(case=case)
enabled = True
#type = 'budget' # not used in plotgen
nx = 256
ny = 256
nz = 100
dxy = 100             # [m]
dz = 40               # [m]
dt = 1                # [s]
startTime = 721.0     # [minutes]
endTime = 1440.0       # [minutes]
startHeight = 0.0     # [m]
endHeight = 4000.0    # [m]
time_3d = 259200.0    # dt

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
sam_file = '../../output/RICO_512x512x100_drizzle.nc'
# nc file generated from .bin3D output
sam_3d_file = '../../output/sam_3d/RICO_256x256x100_drizzle_256_0000259200.nc'

## CLUBB
clubb_zm_file = '../../output/prog_mom_flux_paper/new/rico_zm.nc'
clubb_zt_file = '../../output/prog_mom_flux_paper/new/rico_zt.nc'
## old CLUBB
old_clubb_zm_file = '../../output/prog_mom_flux_paper/old/rico_zm.nc'
old_clubb_zt_file = '../../output/prog_mom_flux_paper/old/rico_zt.nc'

## case setup files
sam_prm = case_folder+'/prm.les'
sam_grd = case_folder+'/grd'

# header in html file
headerText = '{run} {{type}} Minutes {start}-{end}, {bottom}m-{top}m'.format(run=run, start=startTime, end=endTime, bottom=startHeight, top=endHeight)

#-------------------------------------------------------------------------------
#   G E N E R A L   P L O T   S E T T I N G S
#-------------------------------------------------------------------------------
lw = 5
color = 'nipy_spectral'
yLabel = 'Height [m]'
