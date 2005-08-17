Using the HOC tuner
-----------------------------------------------------------------------
Last update August 17, 2005

-----------------------------------------------------------------------
-                                                                     -
- Building the source code:                                            -
-                                                                     -
-----------------------------------------------------------------------

$ cd ~/hoc_v2.2_tuner/src
(edit Makefile for your compiler and optimization options.)
$ make
$ make install

-----------------------------------------------------------------------
-                                                                     -
- Executing a tuning run:                                             -
-                                                                     -
-----------------------------------------------------------------------

1.  cd ~/hoc_v2.2_tuner/tune

2.  Edit <case>_hoc.in for each case you wish to run, or just leave them as is.
  Usually you will want to keep these the same.

3.  Edit error_<runtype>.in or select a premade one. Note that there are two
  tuning subroutines, specified by tune_type in the /stats/ namelist.  If 
  runtype = 0, then the amoeba subroutine, a downhill simplex algorithm, will
  be used.  If runtype is any other value, then amebsa, a variant which uses
  simulated annealing, is used.  A complete explanation of these minimization
  algorithms can be found in _Numerical Recipes in Fortran 90_.
  Sometimes the variable names in hoc's zt and the coamps grads file will 
  differ.
  Currently, it is only possible to tune for variables occuring in zt.

4.  Edit run_tuner.bash to use your error_*.in

5.  ./run_tuner.bash 
-----------------------------------------------------------------------
-                                                                     -
- Executing a budget terms tuning run:                                -
-                                                                     -
-----------------------------------------------------------------------
One run at a time:
1. cd ~/hoc_v2.2_tuner/tune_budgets/ 
2. Check hoc_v2.2_tuner/tune/<CASE>_hoc.in to make sure dtmain equals
  dtclosure.  Make a note of what your statistics file's timestep divided by 
  dtmain is, since this is the sample_ratio in your budget namelist. 
  e.g. If the situation is a 1 mn timestep in LES GRaDS statistics paired 
  with a 5 second dtmain, the sample_ratio is 60.0 / 5.0 = 12.0
  You might also want to verify that your LES file has sufficient data
  for the whole duration of the HOC run.
3. edit <CASE>_budget.in;  You will want to set one case_tune variable to true 
  and the rest to false, or the resultant parameters from one tuning run 
  will be used in subsequent runs.  
  Alteratively, you may edit the hoc_tuner_budget_terms.F code, setting liter 
  to .true., and attempt to optimize the constants over multiple cases.  
  Sometimes (usually?) this doesn't converge on anything.
4. edit run_budget.bash for your model <CASE> 
5. ./run_budget.bash 

Batch mode: 
1.  cd hoc_v2.2_tuner/tune_budgets/
    mkdir <directory name you like>
2.  cp BOMEX/create_files.bash <dir from (1)>
    cp BOMEX/run_batch.bash <dir from (1)>
3.  cp <CASE>_budget <dir from (1)>/budget.tmpl
4.  cd <dir from (1)>
5.  edit budget.tmpl to your liking;  
    remove the 'case_tune' namelist entirely.
5.  ./create_files.bash (you shouldn't see any errors)
6.  ./run_batch.bash
7.  Make coffee. Play spider solitaire. Wait an hour or two.

-----------------------------------------------------------------------
-                                                                     -
- Executing a standalone run:
-                                                                     -
-----------------------------------------------------------------------
Do steps (1) and (2) as outlined above.

3. cd ../standalone.  Edit standalone_<case>.in or select a premade one.

4. Edit run_standalone.bash to use your standalone_*.in

5. ./run_standalone.bash

-----------------------------------------------------------------------
-                                                                     -
- Executing a compare_runs analysis:
-                                                                     -
-----------------------------------------------------------------------

1. cd ../compare_runs  

2. Edit compare_runs.in.  You need to choose three GrADS files on disk
      to compare.  You also need to choose the time intervals over which
      the files will be compared.

3. ./compare_runs
-----------------------------------------------------------------------
-                                                                     -
- Executing a jacobian analysis:
-                                                                     -
-----------------------------------------------------------------------
1. cd ../jacobian

2. Edit jacobian.in.  Choosing a high delta_factor may make the model
   crash, which will result in no data (results for that term will come
   back as infinite).  The model namelists come from ../tune.
3. ./jacobian
-----------------------------------------------------------------------
-                                                                     -
- Overview of the code                                                -
-                                                                     -
-----------------------------------------------------------------------

This code tunes certain parameters in a one-dimensional boundary layer cloud 
parameterization (``hoc''), to best fit large-eddy simulation output.  The 
optimization technique is the simplex algorithm, as implemented in _Numerical 
Recipes_ (amoeba.f90).  The parameterization is called as a subroutine 
( hoc_model() ) with parameter values as input.

The code is highly flexible.  One can vary the cases (bomex, fire, arm, or
atex) to match; the variables to match (cloud fraction, liquid water, third
moment of vertical velocity, etc.); the altitude and times over which to match
these variables; and the parameters to tune (C1, beta, etc.). 

The code is written in Fortran 90/95 and executed by a bash runscript. On the
Microsoft Windows platform this could work using MSYS or Cygwin with G95, but 
this has not been tested.
We use the Portland Group compiler, version 5.2-4.
G95 currently only works -ffast-math off, and hoc_standalone is the only
program which works reliably.
Sun's f95 8.1 on Solaris appears to work for hoc_tuner, but has not been 
rigorously tested. 
The GNU fortran compiler (gcc 4.0.x) does not implement the entire 
Fortran 90 standard yet, and so does not work at all.

-----------------------------------------------------------------------
-                                                                     -
-Explanation of the Input and Output Files                            -
-                                                                     -
-----------------------------------------------------------------------

Nota bene: Our numerical output is in GrADS format.  Each output has a header,
or control file, (.ctl) and a data file (.dat).  The .ctl file describes the
file format, which variables are output in which order, etc.

read_grads_hoc.m:  A MATLAB function that reads GrADS data files.

LES GrADS files:
  les_data/bomex_coamps_sw.ctl, les_data/wangara_rams.ctl
  FIRE, BOMEX, ARM & ATEX are our 4 basic ``datasets'', simulated by COAMPS, 
  that we use to match HOC data.  Each output file includes an hour of 
  unphysical spinup time.  BOMEX is trade-wind cumulus; FIRE is marine 
  stratocumulus; ARM is continental cumulus; and ATEX is cumulus under 
  stratocumulus.  BOMEX, FIRE, and ATEX are statistically steady-state; 
  ARM varies over the course of a day.

HOC GrADS files:
  results/<DATE>/bomex_zt.ctl, etc.
  These are the 4 basic runs using the default constants. Not used in any
  way by the tuner, but useful for comparing runs results in grads.  To obtain
  similar results on differing platforms, check your compiler's documentation
  for information on enabling IEEE 754 standard floating-point arithmetic.

Generated HOC GrADS files:
  bomex_zt.dat, fire_zt.dat, arm_zt.dat, atex_zt.dat, dycoms_zt.dat,
    wangara_zt.dat, <case>_zm, <case>_pdf ...
  These are the files generated by the hoc subroutine during a model run
  and compared to the COAMPS LES results.  The last of these generated
  when tuning will we the optimized results.  Every time the hoc subroutine
  is called, these are overwritten, so if you want to prevent them
  from being erased be sure to either copy the .ctl and .dat
  files to another directory or rename them.
  The tuner currently only uses variables in the zt file.

The Namelist files:

  tune/bomex_hoc.in, fire_hoc.in, arm_hoc.in & atex_hoc.in.
  These files specify the standard hoc model parameters.  Usually these don't
    need to be edited.
  
  tune/error_all.in, error_<CASE>.in, error_<DATE>.in.
  These specify tuning parameters, case information initial constant values, 
    and constant variance.

The Randomization files:

  generate_seed.bash, random_seed.dat, int2txt
  The script uses intrinsic functionality in the Linux kernel to generate
  a pseudo random seed (the .dat) used by the tuner for creating the x_array 
  of constants.  This works on any operating system with a Linux compatible 
  /dev/random (Solaris, Tru64, etc.) as well.  The seed file is now ascii 
  and can be edited by hand.

The compare_runs files:

  compare_runs and compare_runs.in.  This is the executable and namelist for a
  program which compares the variation between two GRaDS files with a number
  of key variables.  Useful for checking the soundness of the model when the
  code is modified and also for getting the raw difference between LES and 
  HOC profiles.

------------------------------------------------------------------------
-
-The new scalar code ( Hoc with -DSCALARS enabled )
-
------------------------------------------------------------------------

By default Hoc should be setup to compile without this option.  To use 
this option, you must modify the Makefile in the src directory so that 
FCFLAGS includes "-DSCALARS" and do a make clean, make, and make install.

Initially, the sclr arrays are configured to contain two vertical columns
containing copies of thl and rt (the first and second elements, respectively).
The code is sufficently general that an arbitrary number of scalars can be
added with a small number of modifications.  The following files must be
adjusted to customize the scalars:

constants.F:  sclrm_dimension is the number of scalars per array and sclrtol 
is used used in the code to diagnose variances. (see diag_var.F for the
algorithm used to )

hoc.F:  The boundary conditions, while fairly general, are setup for thl
in some places as is noted in the code.  Search for SCLR_THETA and SCLR_RT to
find places where the code is not general for both cases.

mixing.F:  When timestep_mixing calls mixing_solve, it needs an argument that
indicates the correct equation to use.  This relies on a case statement within
mixing_solve, and should be easy to add to as needed. Calling the "rtm" and 
"thlm" cases more than once each will cause sampling errors in some of the
budget terms.  Currently budget terms are not setup for rtm/thlm or the
sclr.

statistics.F:  Currently it is only configured to generate data for 2 elements 
for each of the sclr arrays.  Adding more will probably take some time.

sfc.F:  sfc_var has a peculiar constant factor used for calculations on rt
 and thl.  Commenting out the case statement in that subroutine will disable 
this for a new sclr element.
The other subroutine, sfc_thermo_fluxes, only affects the fire case and just
requires a small modification to change its influence on sclrm, found on the
lines that begin with the "if ( present ( sclrm )" statements.

gcss.F: All the cases that you wish to run will require a modification to the
tndcy and sfclyr subroutines, since they will be configured for thl and rt by
default.  Searching for SCLR_RT and SCLR_THETA should come up with all of
them.

All the other source files should work as is.

The Namelists:

Within the existing hoc.in for each run a sounding for the scalar variable
must be added. 

Setting sclrm(:,1) equal to thlm, and sclrm(:,2) equal to rtm in 
the bomex case can be done like so:

&scalar_sounding
sclr(:,1) = 298.7, 298.7, 299.39375, 302.4, 308.4, 313.675
sclr(:,2) = 0.01729, 0.01657, 0.01549, 0.01082, 0.00422, 0.00241
/

This can be appended to the end of the file (in this case tune/bomex_hoc.in),
and it should follow the number of z-levels et cetera found in the &sounding 
namelist.

Finally, if you wish to see the results of your calculations, you will need to
append their names to the vars_zt and var_zm portion of the namelist.
The variables follow the convention of the a=1, and b=2, appended after the
sclr portion of their name.  For example. the first scalar mean is 'sclram',
and the second is 'sclrbm'.  These and their forcings are all that occurs in
the zt file, the rest all occur in the zm file.

------------------------------------------------------------------------
