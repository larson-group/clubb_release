# Steps to repeat the simulations from Storer et al. 2014:
## Easy way
Run the `recreate_storer_simulations.sh` script in the root folder.

## Manual way
### 1. Download SAM output
* SAM files necessary for prescribing radiation in TWP-ICE and ARM97, and
   for plotting comparisons with CLUBB-SILHS are available at:
   https://carson.math.uwm.edu/larson-group/storer/
*  SAM files should be placed in the SAMoutput directory within the trunk.

### 2. Compile and run CLUBB
* compile with `bash ./compile/compile.bash`
* run with `bash ./run_scripts/run_scm_all.bash`, which is set up to run only the 5 cases of interest.
*  Case files referencing the new SAM data for plotting can be found in 
postprocessing/plotgen/cases/clubb/paper/
