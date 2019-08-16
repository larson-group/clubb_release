"""
-------------------------------------------------------------------------------
G E N E R A L   I N F O R M A T I O N
-------------------------------------------------------------------------------
This file contains the case specific input needed for plotgen.py.
to plot height profiles and budgets for the GCSSARM case
"""

#-------------------------------------------------------------------------------
#   G E N E R A L   D E F A U L T   S E T T I N G S
#-------------------------------------------------------------------------------

case = 'GCSSARM'
full_name = case
case_folder = '/home/sdomke/workspace/clubb/sam_clubb/{case}'.format(case=case)
enabled = True # not used in plotgen, no idea what this should do
#type = 'budget' # not used in plotgen
nx = 96
ny = 96
nz = 110
dxy = 66.7            # [m]
dz = 80               # [m]
dt = 1                # [s]
startTime = 301.0     # [minutes]
endTime = 840.0       # [minutes]
startHeight = 0.0   # [m]
endHeight = 3500.0    # [m]
time_3d = 21600.0     # dt

# run entry for header in html file
run = '{case}_{nx}x{ny}x{nz}_{dxy}m_{dz}m_{dt}s'.format(case=case, nx=nx, ny=ny, nz=nz, dxy=dxy, dz=dz, dt=dt)
# jpg file names
plot_case_name = '{case}_{dx}x{dx}_{{type}}_{{date}}_{{plot}}'.format(case=case.lower(), dx=nx)
## use absolute paths or relative paths originating from the directory containing plotgen.py
# directory for output files
out_dir = '../../output/plotgen_out/{case}_{{date}}/'.format(case=case.lower())
# pdf output name
out_pdf = '{case}_{dx}x{dx}_{{type}}_{{date}}.pdf'.format(case=case.lower(),dx=nx)

## input .nc files
## SAM
# nc file generated from .stat output
sam_file = '../../output/GCSSARM_96x96x110_67m_40m_1s.nc'
# nc file generated from .bin3D output
sam_3d_file = '/home/sdomke/workspace/clubb/avi_out/out3d/GCSSARM_96x96x110_67m_40m_1s_96_0000052200.nc'

## CLUBB
clubb_zm_file = '../../output/prog_mom_flux_paper/new/arm_zm.nc'
clubb_zt_file = '../../output/prog_mom_flux_paper/new/arm_zt.nc'
## old CLUBB
old_clubb_zm_file = '../../output/prog_mom_flux_paper/old/arm_zm.nc'
old_clubb_zt_file = '../../output/prog_mom_flux_paper/old/arm_zt.nc'

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