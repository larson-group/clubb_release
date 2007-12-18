$Id: Readme.txt,v 1.43 2007-12-18 22:54:40 nielsenb Exp $

************************************************************************
*                           Copyright Notice
*                         This code is (C) 2006 
*         Jean-Christophe Golaz, Vincent E. Larson, Brian M. Griffin, 
*            David P. Schanen, Adam J. Smith, and Michael J. Falk.
*
*         The distribution of this code and derived works thereof 
*                      should include this notice.
*
*         Portions of this code derived from other sources (ACM TOMS,
*         Numerical Recipes, etc.) are the property of their respective
*         authors as noted and also subject to copyright.
************************************************************************

************************************************************************
*                     Overview of the HOC code
************************************************************************

For a detailed description of the model code see:

``A PDF-Based Model for Boundary Layer Clouds. Part I:
  Method and Model Description'' Golaz, et al. (2002)
  JAS, Vol. 59, pp. 3540--3551.

See also the hoc_v2.2_tuner/doc/hoc_eqns.pdf file in the CVS repository for
finer details on how the discretization was done.

The tuner code tunes certain parameters in a one-dimensional boundary layer 
cloud parameterization (``hoc''), to best fit large-eddy simulation output.  
The optimization technique is the downhill-simplex method of Needler and Mead, 
as implemented in _Numerical Recipes In Fortran 90_ (amoeba.f90).  
The parameterization is called as a subroutine ( hoc_model() ) with 
parameter values as input.

The code is highly flexible.  One can vary the cases (bomex, fire, arm, or
atex) to match; the variables to match (cloud fraction, liquid water, third
moment of vertical velocity, etc.); the altitude and times over which to match
these variables; and the parameters to tune (C1, beta, etc.). 

The code is written in Fortran 90/95 and executed by a bash runscript. 
On the Microsoft Windows platform this could work using MSYS or Cygwin 
with G95, but we have not tested this sort of configuration.

We use the Portland Group compiler, version 5.2-4 on Redhat Enterprise 3.

G95 <http://www.g95.org/> has been tested on SPARC & x86 Solaris,
x64 & x86 GNU/Linux.

Sun Fortran 8.1 and 8.2 on Solaris works, but has not been as rigorously 
tested as pgf90.

Using Intel Fortran 9 we have been able to compile on Linux x86/x64 and
Itanium.

HP/Compaq/DEC Fortran on Alpha also appears to work but because future 
Alpha processor development has ceased it is not extensively tested.

The GNU Fortran compiler (GCC 4.0.x) does not implement the entire 
Fortran 90 standard yet, and so does not work at all.

It is important to note that all these compilers use *incompatible* module
formats for the .mod files!  If you want to use different compilers on the
same system, you will need to build a different set of mod files for each
compiler and use -M or -I to specify their location.

In order to get similar results on differing architectures, platforms, and
compilers, initially try a conservative optimization and enable 
IEEE-754 standard style floating point math.  On x86 compatible machines 
using SSE or SSE2 is usually the best way to do this.

***********************************************************************
*                         Using the HOC Model                         *
***********************************************************************
-----------------------------------------------------------------------
- (1.1) Building Everything:
-----------------------------------------------------------------------

Requirements:
A. A Fortran 90/95 compiler with a complete implementation of the standard.
B. NetCDF >= v3.5.1;  We have not tested our code with anything older.
C. LAPACK & BLAS.  These provide the tri and band diagonal solver
   subroutines needed by HOC.  Many vendors provide optimized versions of
   these routines, which may be much faster than the reference BLAS.
D. GNU bash, or an equivalent POSIX compliant shell to use the run scripts.

Build:
1. $ cd ~/hoc_v2.2_tuner/src
2. Edit a config.<PLATFORM>.in file and choose it in the Makefile for your 
   compiler and optimization options. Note that PREFIX determines where
   your executables and libraries will end up, so make sure you set it
   to the correct location.
3. $ make

The executables will appear in $(PREFIX)/bin and libraries in $(PREFIX)/lib.
The modules remain in the src directory.

If you're using dmake or GNU make and have a fast parallel machine, 
parallel builds should work as well. 

E.g. for 3 threads:
gmake -j 3

-----------------------------------------------------------------------
- (1.2) Building for use in a host model:
-----------------------------------------------------------------------
Requirements:
A. and C. as above.

Build:
1. and 2. as above.

$ make libhoc_param.a

This will build just the static library and the f90 modules.
The static library will be in $(PREFIX)/lib, while the modules will be 
in the src directory.  You will need at least the parameterization_interface 
mod to interface with HOC.

Addition by Brian:  
In addition to the above, you will have to make a reference to the HOC library 
from the configuration file of the host program.  Since HOC now uses the
LAPACK libraries, you will also have to make reference to those.  Currently, 
we do not include the LAPACK libraries with the HOC download.  You will have 
to find them and then download them onto your own computer if they are not
included with your operating system or compiler.  Once you have done this, you 
can reference them in a line such as the following:

-L/home/griffinb/hoc_v2.2_tuner/lib -lhoc_param -llapack -lblas

If the LAPACK and BLAS libraries were compiled with GNU Fortran 77, you may 
need to link to the runtime libs for that with -lg2c as well.

Don't forget that you will also have to make reference
to the HOC src code.  You can reference that with a line
such as the following:

-I/home/griffinb/hoc_v2.2_tuner/src

-----------------------------------------------------------------------
- (1.3) Making clean (or starting over)  
-----------------------------------------------------------------------

Occasionally, one needs to erase old executables or libraries and re-compile 
the code starting with nothing.  For instance, this may be required when 
a library or compiler is updated.  

To delete old object files (*.o), libraries (*.a), and mod files,
go to hoc_v2.2_tuner/src and type

make clean

If this doesn't help, then to additionally delete everything in the binary 
and library directories, go to hoc_v2.2_tuner/src and type

make distclean


-----------------------------------------------------------------------
- (2.1) Executing a standalone run:
-----------------------------------------------------------------------

1.  cd ~/hoc_v2.2_tuner/model

2. Edit <CASE>_model.in for each case you wish to run, or just leave 
   them as is.  Usually you will want to keep these the same.
   See the rain code for description of kk_rain and cloud_sed.
   See BUGSrad description below for a description of the interactive
   radiation scheme.
   Enabling any of these flags may increase runtime considerably.

3. cd ../stats
   Edit <CASE>_stats.in for each case.  A complete list of all computable
   statistics is found in statistics.F.  Note that HOC now supports GrADS or
   NetCDF, but you can only tune using GrADS.

4. $ cd ../standalone.  Edit standalone_<CASE>.in or select a premade one.

5. $ ./run_standalone.bash <CASE>

-----------------------------------------------------------------------
- (3.1) Executing a tuning run:
-----------------------------------------------------------------------

Do steps 1, 2, & 3 as outlined in the standalone run.

4.  Edit error_<runtype>.in or select an existing one. Note that there are two
  tuning subroutines, specified by tune_type in the /stats/ namelist.  If 
  runtype = 0, then the amoeba subroutine, a downhill simplex algorithm, will
  be used.  If runtype is any other value, then amebsa, a variant which uses
  simulated annealing, is used.  A complete explanation of these minimization
  algorithms can be found in _Numerical Recipes in Fortran 90_.
  Sometimes the variable names in HOC's zt and the LES grads files will 
  differ.
  Currently, it is only possible to tune for variables that occur in zt.

5.  Edit run_tuner.bash to use your namelists

6.  ./run_tuner.bash

-----------------------------------------------------------------------
- (3.1.1) Creating a RAM disk (optional)
-----------------------------------------------------------------------

One means of speeding up tuning runs is reducing the time spent writing
to the hard disk.  Most operating systems support a virtual device called
a ram disk, which is main memory that has been allocated to act as an emulated
file system.  Note that you will need system privileges to make the ram disk, 
and files copied to the ram disk are not preserved when the computer is 
powered off.

Generally:

1. mkdir <HOC PATH>/rd_tune/

2. Create and mount RAM disk on rd_tune

3. Copy tune directory to rd_tune

4. Run tuner

Linux Example
Note that you will need ram disk support compiled into your kernel, which is
typically the default on most systems.  Linux appears to be less flexible 
about hen you are allowed to change the ramdisk size.

1. In grub.conf
   Append to 'kernel' line:
   kernel /vmlinuz-2.4.21-40.EL ro root=LABEL=/ ramdisk_size=262144

Sets ram disks to be 256 megabytes in size.  Note that your own system may
have other options besides the ramdisk_size.

2. $ mkfs.ext2 /dev/ram0

3. $ mount /dev/ram0 /home/dschanen/hoc_v2.2_tuner/rd_tune

4. $ cp tune/*.* rd_tune/

5. $ cd rd_tune

(Run your job)

Solaris Example
Note that these instructions are for Solaris 9 & 10

1. $ ramdiskadm -a hoc 256m 
Creates a virtual disk hoc that is 256 megabytes in size.

2. $ newfs /dev/ramdisk/hoc

3. $ mount /dev/ramdisk/hoc /home/dschanen/hoc_v2.2_tuner/rd_tune/

4. $ cp tune/*.* rd_tune/

5. $ cd rd_tune

(Run your job)

-----------------------------------------------------------------------
- (3.2) Executing a budget terms tuning run:
-----------------------------------------------------------------------

One run at a time:
1. cd ~/hoc_v2.2_tuner/tune_budgets/ 
2. Check hoc_v2.2_tuner/model/<CASE>_model.in to make sure dtmain equals
  dtclosure.  Make a note of what your statistics file's timestep divided by 
  dtmain is, since this is the sample_ratio in your budget namelist. 
  e.g. If the situation is a 1 mn timestep in LES GRaDS statistics paired 
  with a 5 second dtmain, the sample_ratio is 60.0 / 5.0 = 12.0
  You might also want to verify that your LES file has sufficient data
  for the whole duration of the HOC run.
3. Edit <CASE>_budget.in;  You will want to set one case_tune variable to true 
  and the rest to false, or the resultant parameters from one tuning run 
  will be used in subsequent runs.  

  Alternatively, you may edit the hoc_tuner_budget_terms.F code, setting liter 
  to .true., and attempt to optimize the constants over multiple budget terms.  
  This is some experimental code that starts at 10x's the ftol and loops over
  all the amoeba cases, updating the constants by small amounts.  The idea
  being that a ``best fit'' set of constant value will eventually be arrived
  at when the tolerance is lowest( at ftol) on the 10th iteration of the 
  do loop.

  Sometimes (usually?) this mode doesn't appear to converge on anything.  
  The more cases this is run over, the more likely it appears the constants 
  will descend into values that cause the model become invalid and it will
  fail altogether.

4. Edit run_budget.bash for your model <CASE> 
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
- (3.3) Executing an ensemble tuning run:
-----------------------------------------------------------------------

NOTES AND INSTRUCTIONS FOR THE HOC ENSEMBLE TUNER
-------------------------------------------------

Go to the main directory (hoc_v2.2_tuner/) and follow these instructions:

1)  Copy the ens_tune directory to a new directory with a slightly different
    name, such as ens_tune_two.  From this point on, I will refer to this new
    directory as ens_tune_xyz.

2)  Enter the new directory (ens_tune_xyz).

3)  Decide which cases will be tuned for (ex. DYCOMS2 RF02 DO, DYCOMS2 RF01,
    FIRE, ARM, BOMEX, etc.).

4)  For EACH case being tuned for:

     a)  Copy the <CASE>_model.in file from the model/ directory to the
         ens_tune_xyz/ directory.  Make sure that the file is set up
         correctly.

     b)  Create/copy the <CASE>_stats_tune.in file.  A sample of this file can
         be found in the ens_tune/ directory (and subsequently in the
         ens_tune_xyz/ directory if it was copied from the ens_tune/ directory).
         This file is used in the tuning process itself.  In order for this
         process to work, stats_fmt in this file must be set to 'grads'.  The
         tuner code reads the LES GrADS files and compares the results with the
         HOC GrADS files (which are written according the specifications stated
         in this file).  This file MUST also contain the names of the HOC
         variables that are being tuned for.  During the process of tuning,
         small GrADS files will be written that contain the results for ONLY the
         variables that are listed here.  This makes the tuning process go
         faster because it is not being slowed down by the writing of
         unnecessary variables.

     c)  Copy the <CASE>_stats.in file from the stats/ directory to the
         ens_tune_xyz/ directory.  Once the tuner has found the optimal value
         and the tuning process has finished, standalone HOC will run for
         each case that was tuned for with the values of the constants that the
         tuner found.  The statistics that it finally outputs will be directed
         according to this file.  This file produces the normal statistical
         files that one sees after a normal HOC run.  The main thing that might
         be changed is whether the final statistical output is in GrADS or in
         netCDF.  The user can choose either one.  The user just has to set
         stats_fmt to either 'grads' or 'netcdf'.  Of course, if the user
         desires, the sampling and output timesteps or the variables that are
         output can also be changed.

5)  Edit the error_messner_001.in file:

         The error_messner_001.in file is the same type of file as the error.in
         file in the regular HOC tuner.  This file must include all the relevant
         information for EVERY case being tuned for, such as the HOC and LES
         stats files, the HOC run file, the vertical levels and time periods
         being tuned for, and the general weighting of each case.  This file
         also must include other important factors such as the variable(s) being
         tuned for and the general weighting of each variable.  Of course, the
         initial values of the constants and the amount of deviation allowed for
         each of the constants are also declared here.

         The runscript copies any error*.in file over to the remote nodes for
         running (as error.in).  So, the last file copied would end up being
         the file read in the running of the tuner.  Therefore, make sure that
         there is only ONE error*.in file in the ens_tune_xyz directory when
         you start the ensemble tuning run.

6)  Edit the mytuner_messner.bash file:

         Only the Global variables section of this file needs to be edited.
         There are 6 global variables to edit:

         a) EXPERIMENTS:  List the names of the cases that you are tuning for.
                          This should match the <CASE> model, stats, and
                          stats_tune files that were copied or created and the
                          listing of cases that is found in error*.in

         b) ARCHIVE:  The general path to a directory to put the results in.
                      Usually, many specific results subdirectories reside
                      within this declared main directory.

         c) CASE:  The subdirectory within ARCHIVE to put this tuning runs
                   specific results in.

         d) HOC:  The path to the HOC main level directory.  This is the
                  directory above ens_tune_xyz, where the tuner is being run
                  from.

         e) NODES:  The nodes on the Messner cluster which are going to run
                    these tuning iterations.

         f) ITERMAX:  The total number of tuning iterations to run.  ITERMAX
                      needs to be a multiple of the total number of nodes
                      that are being used in the tuning run.

         Whenever you are setting up a run, you need to take into account the
         amount of disk space you will need.  You need to take into account
         the size of all the LES data files for each case you are tuning for.
         You also need to take into account the size of the binary
         (executable) files that you will need.  These files include hoc_tuner,
         hoc_standalone, and int2txt.  Finally, you need to know the size of
         the HOC output files for each of the cases you are tuning for.  These
         files will be generated when the individual tuning iterations are
         done.  One set of the HOC output files will be generated for each
         tuning iteration that you run.

         For example, let's say that the executable files sum to about 25 MB.
         You decide to tune for four cases -- ARM, BOMEX, DYCOMS2 RF02 DO, and
         DYCOMS2 RF01.  Let's say that the LES data files for these four cases
         sum to about 200 MB.  So, that's already 225 MB of space needed on
         EACH NODE used in the tuning run.  Now, let's say that the netCDF
         output files for the HOC results for these four cases sum to about
         150 MB.  However, you are running 25 iterations of the tuner on EACH
         NODE.  So, you will be producing those 150 MB worth of HOC files 25
         times over.  The total space needed on the node is then:

         25 MB + 200 MB + (150 MB/iteration)*(25 iterations) = 3975 MB

         So, about 4 GB of available disk space would be needed on each node
         in order to complete the run.  You would have to check the available
         disk space in the home directory on every node that you are using in
         order to make sure that EACH NODE has enough space.

         Say that you were running the above ensemble tuning run on nodes
         1-10.  That would yield a grand total of 250 iterations of the HOC
         tuner.  The first 10 runs would be launched (one on each node).  When
         ALL of them are completed, then the next set of 10 runs are launched,
         and so on and so on.  When the last set is finally complete, all the
         HOC results are copied back to Messner, and then deleted off the
         remote tom nodes.  The amount of disk space you would need on Messner
         is:

         (150 MB/iteration)*(250 iterations) = 37500 MB

         So, you would need almost 40 GB of available space on Messner to hold
         all the results.

7)  Starting the ensemble tuning run:

         A simple ./mytuner_messner.bash would start the tuning run.  However,
         that would require leaving a session open until the job is finished
         (which could take days).  Therefore, the special atjob.bash script
         has been created in order to have the job run in the background.
         The proper command is:  at now -f ./atjob.bash

8)  Sorting the results:

         The script sortresults.bash will produce output to the screen that
         orders the results by value of the cost function, from lowest value
         (= lowest error) to highest value.  The iteration number is listed
         along with the cost function value.  This script also produces a
         text file called results.txt.  However, this files orders the output
         according to iteration number, rather than from best to worst.

         Usually, only the top so many values produce good results.  It is
         best to look at the results for the top handful of tuning iterations.
         First, to see if they look good for all the cases that were tuned for,
         and then to see if they look good for all the cases we have in HOC.

-----------------------------------------------------------------------
- (4.1) Executing a run comparison analysis:
-----------------------------------------------------------------------

1. $ cd ../compare_runs  

2. Edit compare_runs.in.  You need to choose three GrADS files on disk
   to compare.  If necessary, uncompress the files you choose using
   the command $ bunzip2 filename.bz2.  You also need to choose the 
   time intervals over which the files will be compared.

3. $ ../bin/compare_runs

-----------------------------------------------------------------------
- (5.1) Executing a Jacobian analysis:
-----------------------------------------------------------------------

1. $ cd ../jacobian

2. Edit jacobian.in. 
   cat ../model/<model name>_model.in ../stats/<model name>_stats.in \
   >  <model name>_hoc.in
   Note that choosing a high delta_factor may make the model
   crash, which will result in no data (results for that term will come
   back as infinite).  The model namelists come from ../model.
3. $ ../bin/jacobian


-----------------------------------------------------------------------
- (1.1) Explanation of the Input and Output Files
-----------------------------------------------------------------------

Nota bene: Our numerical output is usually in GrADS format.  Each output 
has a header, or control file, (.ctl) and a data file (.dat).  The .ctl file 
describes the file format, which variables are output in which order, etc.

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
    wangara_zt.dat, <case>_zm, <case>_sfc ...
  These are the files generated by the hoc subroutine during a model run
  and compared to the COAMPS LES results.  The last of these generated
  when tuning will we the optimized results.  Every time the hoc subroutine
  is called, these are overwritten, so if you want to prevent them
  from being erased be sure to either copy the .ctl and .dat
  files to another directory or rename them.
  The tuner currently only uses variables in the zt file.

The Namelist files:

  model/bomex_model.in, fire_model.in, arm_model.in & atex_model.in.
  These files specify the standard HOC model parameters.  Usually these 
  do not need to be modified.

  stats/bomex_stats.in, fire_stats.in, arm_stats.in & atex_stats.in.
  These files specify statistics output for each simulation.  See statistics.F
  for a complete list of the all output supported.
  
  tune/error_all.in, error_<CASE>.in, error_<DATE>.in.
  These specify tuning parameters, case information initial constant values, 
    and constant variance.

The Randomization files:

  generate_seed.bash, rand_seed.dat, bin/int2txt
  The script uses intrinsic functionality in the Linux kernel to generate
  a pseudo random seed (the .dat) used by the tuner for randomizing initial
  parameters.  This works on any operating system with a Linux style 
  /dev/random (Solaris, Tru64, etc.) as well.  The seed file is now ASCII 
  text and can be edited by hand.

The compare_runs files:

  compare_runs and compare_runs.in.  This is the executable and namelist for a
  program which compares the variation between two GRaDS files with a number
  of key variables.  Useful for checking the soundness of the model when the
  code is modified and also for getting the raw difference between LES and 
  HOC profiles.

------------------------------------------------------------------------
- (2.1) The BUGSrad Radiation scheme
------------------------------------------------------------------------

  This is an optional more complex radiation scheme, developed apart from
  HOC by Stephens, et al. The code used in HOC was obtained from Norm Wood 
  on 2004/07/10.
  When enabled, the analytic computation normally
  used for radiation is disabled.  BUGSrad is enabled in the 
  model/<RUN CASE>_model.in file by setting lbugsrad=.true.
  Currently,  all cases where bottom_at_sfc is false, including November 11,
  will give inaccurate results due to our interface's
  inability to add lower altitude levels.  Other cases appear to give
  plausible results, comparable with the analytic code.

  BUGSrad allows the output of the following variables:

  Momentum grid:
  Frad, Frad_SW, Frad_LW:  Radiative Flux; Short-wave/Long-wave component;

  Thermodynamic grid:
  radht, radht_SW, radht_LW:  Radiative Heat; Short-wave/Long-wave component;

  The thlm forcing will also be influenced by this calculation.

  Note that for most cases SW and LW components are not calculated 
  without using BUGSrad.

------------------------------------------------------------------------
- (3.1) The new scalar code ( HOC with -DSCALARS enabled )
------------------------------------------------------------------------

The scalars in the code provide a generalized way of simulating a 
passive scalar in the atmosphere.

By default HOC should be setup to compile without this option.  To use 
this option, you must modify the Makefile in the src directory so that 
FCFLAGS includes "-DSCALARS", and then do a make clean, and then make.  
To output the scalar fields from a HOC simulation, be sure to include 
them in the appropriate stats/*_stats.in file.

Currently the code contains eddy-diffusivity scalar code and the more complex
code used in diag_var(), closure_new(), and timestep_mixing().  Both use two
dimensional arrays, but the code and results for each is separated.

Initially, the sclr arrays are configured to contain two vertical columns
containing copies of thl and rt (the first and second elements, respectively).
The code is sufficiently general that an arbitrary number of scalars can be
added with a small number of modifications.  The following files must be
adjusted to customize the scalars:

constants.F:  sclrm_dimension is the number of scalars per array and sclrtol 
is used used in the code to diagnose variances. (see diag_var.F for the
algorithm used).

hoc.F & parameterization_interface.F:  
The boundary conditions, while fairly general, are setup for thl
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
default.  Searching for SCLR_RT and SCLR_THETA should locate the code for all 
of them.

All the other source files should work as is.

The edsclrs are only computed in subroutine update().

The Namelists:

Within the existing _model.in for each run a sounding for the scalar variable
must be added. 

Setting sclrm(:,1) equal to thlm, and sclrm(:,2) equal to rtm in 
the bomex case can be done like so:

&scalar_sounding
sclr(:,1) = 298.7, 298.7, 299.39375, 302.4, 308.4, 313.675
sclr(:,2) = 0.01729, 0.01657, 0.01549, 0.01082, 0.00422, 0.00241
/

This can be appended to the end of the file (in this case model/bomex_model.in),
and it should follow the number of z-levels et cetera found in the &sounding 
namelist.

Finally, if you wish to see the results of your calculations, you will need to
append their names to the vars_zt and var_zm portion of the stats namelist.
The variables follow the convention of the a=1, and b=2, appended after the
sclr portion of their name.  For example. the first scalar mean is 'sclram',
and the second is 'sclrbm'.  These and their forcings are all that occurs in
the zt file, the rest all occur in the zm file.
*************************************************************************
