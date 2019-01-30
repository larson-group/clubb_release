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
# TODO: check parameters in Heinze/Siebesma and (auto-)rename name,case,out_dir etc

case = 'DYCOMS2_RF02'
enabled = True # not used in plotgen, no idea what this should do
#type = 'budget' # not used in plotgen
nx = 128
ny = 128
nz = 96
dxy = 50              # [m]
dz = '5-80'           # [m]
dt = .5               # [s]
startTime = 0.0      # [minutes]
endTime = 60.0       # [minutes]
startHeight = 20.0    # [m]
endHeight = 1200.0    # [m]

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
sam_file = '/home/sdomke/workspace/clubb/avi_out/dycoms2_dbg_out_181218.nc'
#sam_file = '/home/sdomke/workspace/clubb/sam_out/outstat/DYCOMS_RF02_64x64x96_7200_dr_nosed_181205.nc'
# header in html file
headerText = '{run} {{type}} Minutes {start}-{end}, {bottom}m-{top}m'.format(run=run, start=startTime, end=endTime, bottom=startHeight, top=endHeight)

#-------------------------------------------------------------------------------
#   G E N E R A L   P L O T   S E T T I N G S
#-------------------------------------------------------------------------------
lineWidth = 2
color = 'nipy_spectral'
yLabel = 'Height [m]'