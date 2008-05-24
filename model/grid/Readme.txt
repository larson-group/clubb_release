$Id: Readme.txt,v 1.1 2008-05-24 22:10:57 griffinb Exp $

MODEL/GRID DIRECTORY OVERVIEW
=============================

This directory contains files that can be read in as grid altitude inputs on
either momentum levels (*zm_grid.dat files) or thermodynamic levels 
(*zt_grid.dat files), depending if either grid_type = 3 or grid_type = 2,
respectively, was selected in the appropriate model.in file.

Note:  To use an evenly-spaced grid, simply select grid_type = 1 in the 
       appropriate model.in file.  Also enter the grid spacing (deltaz) and
       the lowest altitude level (zm_init) (found at momentum level 1) in the 
       appropriate model.in file.

The number of entries in the *zm_grid.dat or *zt_grid.dat file needs to match
the number of vertical levels (nzmax) declared in the appropriate model.in 
file.

All *zm_grid.dat or *zt_grid.dat files need to be in written in the format of 
one altitude entry per line in ascending order, with the lowest grid altitude
listed as the first entry in the file.

The following is a list of cases and files that have stretched grids (please
update this list as stretched grids are changed or added):

--------------------------------------------------------------------------------

ARM 3-year (arm_3year):

arm_3year_zt_grid.dat file:

The ARM 3-year zt_grid.dat file is based on the stretched grid that is found 
in SAM for the ARM 97 case.  There are 128 grid levels that cover the 27000 m. 
vertical domain.  The lowest grid levels are at -12 m. and 12 m.  Grid levels 
are 38 m. to 65 m. apart near the surface and gradually increase in thickness 
until they reach 250 m. apart at an altitude of 5750 m.  Above that level, the 
grid levels remain 250 m. apart.

--------------------------------------------------------------------------------

ARM 0003 (arm_0003):

arm_0003_zt_grid.dat file:

The ARM 0003 zt_grid.dat file is based on the stretched grid that is found 
in SAM for the ARM 97 case.  There are 128 grid levels that cover the 27000 m. 
vertical domain.  The lowest grid levels are at -12 m. and 12 m.  Grid levels 
are 38 m. to 65 m. apart near the surface and gradually increase in thickness 
until they reach 250 m. apart at an altitude of 5750 m.  Above that level, the 
grid levels remain 250 m. apart.

--------------------------------------------------------------------------------

ARM 97 (arm_97):

arm_97_zt_grid.dat file:

The ARM 97 zt_grid.dat file is based on the stretched grid that is found in SAM
for that case.  There are 128 grid levels that cover the 27000 m. vertical 
domain.  The lowest grid levels are at -12 m. and 12 m.  Grid levels are 38 m. 
to 65 m. apart near the surface and gradually increase in thickness until they 
reach 250 m. apart at an altitude of 5750 m.  Above that level, the grid levels 
remain 250 m. apart.

--------------------------------------------------------------------------------

DYCOMS-II RF02 
(dycoms2_rf02_do, dycoms2_rf02_ds, dycoms2_rf02_nd, and dycoms2_rf02_so):

dycoms2_rf02_*_zm_grid.dat files:

The DYCOMS-II RF02 zm_grid.dat files are taken from the results of running the 
grid-generating source code found on the DYCOMS-II RF02 website.  The DYCOMS-II 
RF02 intercomparison specifies a stretched grid that is based on 97 momentum 
levels covering a 1500 m. vertical domain.  The thermodynamic levels are placed 
halfway between the momentum levels.


dycoms2_rf02_*_zt_grid.dat files:

The DYCOMS-II RF02 zt_grid.dat files are based on the stretched grid that is 
found in SAM for that case.  The stretched grid found in SAM is generated from 
the same source code found at the DYCOMS-II RF02 website as mentioned above, but
altered so that it takes 148 momentum levels to cover the 1500 m. vertical 
domain.  Again, the source code places the thermodynamic levels halfway between 
the momentum levels.  However, for purposes of SAM, the thermodynamic levels are
output from the source code.  SAM reads the thermodynamic levels as input, and 
then places it's momentum levels halfway between the thermodynamic levels.  
Thus, the resulting momentum levels in SAM are at slightly different altitudes 
than the momentum levels output by the original source code.

--------------------------------------------------------------------------------

LBA (lba):

lba_zt_grid.dat file:

The LBA zt_grid.dat file is based on the stretched grid that is found in SAM
for that case.  The same stretched grid is found in SAM for the 
GATE_noshear_normalf case.  There are 128 grid levels that cover the 27000 m.
vertical domain.  The lowest grid levels are at -25 m. and 25 m.  Grid levels 
are 25 m. to 65 m. apart near the surface and gradually increase in thickness 
until they reach 250 m. apart at an altitude of 5750 m.  Above that level, the 
grid levels remain 250 m. apart.

--------------------------------------------------------------------------------
