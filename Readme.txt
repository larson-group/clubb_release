Using the HOC tuner
-----------------------------------------------------------------------
Last update May 18, 2005

-----------------------------------------------------------------------
-                                                                     -
-Building the source code:                                            -
-                                                                     -
-----------------------------------------------------------------------

$ cd ~/hoc_v2.1_tuner/src
edit Makefile for your compiler and optimization options.
$ make

-----------------------------------------------------------------------
-                                                                     -
-Executing a tuning run:                                              -
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
-Executing a standalone run:
-                                                                     -
-----------------------------------------------------------------------
Do steps (1) and (2) as outlined above.

3. cd ../standalone.  Edit standalone_<case>.in or select a premade one.

4. Edit run_standalone.bash to use your standalone_*.in

5. ./run_standalone.bash

-----------------------------------------------------------------------
-                                                                     -
-Executing a compare_runs analysis:
-                                                                     -
-----------------------------------------------------------------------

1. cd ../compare_runs.  

2. Edit compare_runs.in.  You need to choose three GrADS files on disk
      to compare.  You also need to choose the time intervals over which
      the files will be compared.

3. ./compare_runs

-----------------------------------------------------------------------
-                                                                     -
- Overview of the code                                                -
-                                                                     -
-----------------------------------------------------------------------

This code tunes certain parameters in a one-dimensional boundary layer cloud 
parameterization (``hoc''), to best fit large-eddy simulation output.  The 
optimization technique is the simplex algorithm, as implemented in _Numerical 
Recipes_ (amoeba.f90).  The parameterization is called as a subroutine 
(hoc_sub) with parameter values as input.

The code is highly flexible.  One can vary the cases (bomex, fire, arm, or
atex) to match; the variables to match (cloud fraction, liquid water, third
moment of vertical velocity, etc.); the altitude and times over which to match
these variables; and the parameters to tune (C1, beta, etc.). 

The code is written in Fortran 90 and executed by a bash runscript.  We use the
Portland Group compiler (tested with versions 5.1 & 5.2).
G95 appears to work well, as well as Sun Studio 10*.  gfortran does not.

(*) For some reason the constants printed during iterations come out in
 reverse byte order.  This is a somewhat confusing, but apparently harmless
 bug in the f95 compiler.

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
  bomex_case.dat, fire_case.dat, arm_case.dat, atex_case.dat.
  These are our 4 basic ``datasets," simulated by COAMPS, that we use to match 
  HOC data.  Each output file includes an hour of unphysical spinup time.  
  BOMEX is trade-wind cumulus; FIRE is marine stratocumulus; ARM is continental 
  cumulus; and ATEX is cumulus under stratocumulus.  BOMEX, FIRE, and ATEX are 
  statistically steady-state; ARM varies over the course of a day.

HOC GrADS files:
  untuned-hoc_bomex, etc.
  These are the 4 basic runs using the default constants. Not used in any
  way by the tuner, but useful for comparing runs results in grads.
  untuned-hoc_bomex-zm, etc.
  These are the zm counterparts to the above.  We only use them for 
  comparing wp2 currently.

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

  bomex_hoc.in, fire_hoc.in, arm_hoc.in & atex_hoc.in.
  These files specify the standard hoc model parameters.  Usually these don't
    need to be edited.
  
  error_all.in, error_<CASE>.in, error_<DATE>.in.
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
