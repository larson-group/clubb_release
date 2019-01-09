"""
-------------------------------------------------------------------------------
G E N E R A L   I N F O R M A T I O N
-------------------------------------------------------------------------------
This file contains the case specific input needed for plotgen.py.
to plot height profiles and budgets for the BOMEX case
"""

#-------------------------------------------------------------------------------
#   G E N E R A L   D E F A U L T   S E T T I N G S
#-------------------------------------------------------------------------------
# TODO: check parameters in Heinze/Siebesma and (auto-)rename name,case,out_dir etc

case = 'BOMEX'
enabled = True # not used in plotgen, no idea what this should do
#type = 'budget' # not used in plotgen
nx = 64
ny = 64
nz = 75
dxy = 100             # [m]
dz = 40               # [m]
dt = 1                # [s]
startTime = 181.0     # [minutes]
endTime = 360.0       # [minutes]
startHeight = 0.0     # [m]
endHeight = 2500.0    # [m]

# run entry for header in html file
run = '{case}_{nx}x{ny}x{nz}_{dxy}m_{dz}m_{dt}s'.format(case=case, nx=nx, ny=ny, nz=nz, dxy=dxy, dz=dz, dt=dt)
# jpg file names
plot_case_name = '{case}_{{type}}_{{date}}_{{plot}}'.format(case=case.lower())
## use absolute paths or relative paths originating from the directory containing plotgen.py
# directory for output files
out_dir = '/home/sdomke/workspace/plotgen_out/{case}_{{date}}/'.format(case=case.lower())
# pdf output name
pdf_out = '{case}_{{type}}_{{date}}'.format(case=case.lower())

# input .nc file
sam_file = '/home/sdomke/workspace/clubb/avi_out/BOMEX_64x64x75_100m_40m_1s_181115.nc'
# header in html file
headerText = '{run} {{type}} Minutes {start}-{end}, {bottom}m-{top}m'.format(run=run, start=startTime, end=endTime, bottom=startHeight, top=endHeight)

#-------------------------------------------------------------------------------
#   G E N E R A L   P L O T   S E T T I N G S
#-------------------------------------------------------------------------------
lineWidth = 2
color = 'nipy_spectral'
yLabel = 'Height [m]'