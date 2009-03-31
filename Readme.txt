$Id: Readme.txt,v 1.72 2008-08-01 19:10:01 vlarson Exp $

************************************************************************
*                           Copyright Notice
*                         This code is (C) 2006-2009 
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
*                     Overview of the CLUBB code
************************************************************************

For a detailed description of the model code see:

``A PDF-Based Model for Boundary Layer Clouds. Part I:
  Method and Model Description'' Golaz, et al. (2002)
  JAS, Vol. 59, pp. 3540--3551.

See also the ./doc/hoc_eqns.pdf file in the svn repository for
finer details on how the discretization was done.

The single column model executable ("clubb_standalone") runs a particular 
case (e.g. BOMEX) and outputs statisistical data in either GrADS or 
netCDF format.

The tuner code tunes certain parameters in a one-dimensional boundary layer 
cloud parameterization (``CLUBB''), to best fit large-eddy simulation output.  
The optimization technique is the downhill-simplex method of Needler and Mead, 
as implemented in _Numerical Recipes In Fortran 90_ (amoeba.f90).  
The parameterization is called as a subroutine ( run_clubb() ) with 
parameter values as input.

The tuner code is highly flexible.  One can vary the cases (bomex, fire, arm, 
or atex) to match; the variables to match (cloud fraction, liquid water, third
moment of vertical velocity, etc.); the altitude and times over which to match
these variables; and the parameters to tune (C1, beta, etc.). 

The code is written in ISO Fortran 95 and executed by scripts written for
GNU Bash. 
The mkmf Makefile generating script and some other optional code checking
scripts are written in Perl.
On the Microsoft Windows platform the CLUBB parameterization could be configured
and compiled using MSYS or Cygwin with G95, but we have not tested this sort 
of configuration.

We use the G95 compiler on Intel x64 processors running Redhat Enterprise 5.

G95 <http://www.g95.org/> has been tested on SPARC & x86 Solaris,
x64 & x86 GNU/Linux.

Sun Fortran 8.x on Solaris SPARC/x86 work, but has not been as rigorously 
tested as G95.

Using Intel Fortran 9 we have been able to compile on Linux x86/x64 and
Itanium.

HP/Compaq/DEC Fortran on Alpha also appears to work but because future 
Alpha processor development has ceased it is not extensively tested.

The GNU Fortran compiler (GCC 4.1.x) may or may not work.  The version that
comes with RHEL 5 does not.

It is important to note that all these compilers use *incompatible* module
formats for the .mod files!  If you want to use different compilers on the
same system, you will need to build a different set of netCDF mod files for
each compiler and use -M or -I flags to specify their location.

In order to get similar results on differing architectures, platforms, and
compilers, initially try a conservative optimization and enable 
IEEE-754 standard style floating point math.  On x86 compatible machines 
using SSE or SSE2 is usually the best way to do this.

***********************************************************************
*                        Using the CLUBB Model                        *
***********************************************************************
-----------------------------------------------------------------------
- (1.1) Building (compiling) everything:
-----------------------------------------------------------------------

Requirements:
A. A Fortran 95 compiler with a complete implementation of the standard.
B. GNU make (we use v3.81).
C. Perl 5 to run mkmf
D. LAPACK & BLAS.  These provide the tri and band diagonal matrix solver
   subroutines needed by CLUBB.  Many vendors provide optimized versions of
   these routines, which may be much faster than the reference BLAS.
E. GNU bash, or an equivalent POSIX compliant shell to use the run scripts.

Optionally:
F. GrADS for viewing the GrADS output data.
G. NetCDF >= v3.5.1;  We have not tested our code with anything older.
   If you do not use netCDF, you must remove -DNETCDF from the preprocessor
   flags, found in the compile/config/<platform>.bash file.
H. MATLAB or NCAR graphics for viewing the netCDF output data.

To build, perform the following three steps:
1. $ cd <BASE DIRECTORY>/compile
2. Edit a ./config/<PLATFORM>.bash file and uncomment it the file
   compile.bash. Depending on your platform you may need to create a new
   file based on the existing configurations.
   Note that the variables libdir and bindir determine where
   your executables and libraries will end up, so make sure you set it
   to the correct location (the default is one directory up).
3. $ ./compile.bash

The executables and Makefile will appear in <PREFIX>/bin and libraries in 
<PREFIX>/lib.
The object (.o) and module (.mod) files will appear in <PREFIX>/obj.

If you're using GNU make and have a fast parallel machine, parallel builds 
should work as well.  E.g. for 3 threads, append gmake="gmake -j 3" to the
file source'd from compile.bash.

The mkmf script may or may not generate files that are compatible with
other versions of make.
-----------------------------------------------------------------------
- (1.2) Building for use in a host model:
-----------------------------------------------------------------------

The build procedure is slightly different if you have implemented CLUBB
in a large-scale weather or climate model.

Requirements:
A., B., C., & D. as above.

Build:
1, 2, & 3 as above.

You can safely remove everything but libclubb_param.a from the "all" section
of the compile.bash script if you only want to use CLUBB in a host model.

$ ./compile.bash

This will build just the static library and the f90 modules.
The static library will be in <PREFIX>/lib, while the modules will be 
in the obj directory.  You will need at least the clubb_core.mod &
constants.mod file to interface with CLUBB.

Addition by Brian:  
In addition to the above, you will have to make a reference to the CLUBB 
library from the configuration file of the host program.  Since CLUBB now uses 
the LAPACK libraries, you will also have to make reference to those. Currently, 
we do not include the LAPACK libraries with the CLUBB download.  You will have 
to find them and then download them onto your own computer if they are not
included with your operating system or compiler.  Once you have done this, you 
can reference them in a line such as the following:

-L/home/griffinb/clubb/lib -lclubb_param -llapack -lblas

If the LAPACK and BLAS libraries were compiled with GNU Fortran, you may 
need to link to the runtime libs for that with -lg2c or -lgfortran as well.

Don't forget that you will also have to make reference
to the CLUBB modules.  You can reference that with a line
such as the following:

-I/home/griffinb/clubb/obj

-----------------------------------------------------------------------
- (1.3) Making clean (for re-compiling from scratch)  
-----------------------------------------------------------------------

Occasionally, one needs to erase old executables or libraries and re-compile 
the code starting with nothing.  For instance, this may be required when 
a library or compiler is updated.  

To delete old object files (*.o), and mod (*.mod) files,
go to <PREFIX>/bin (where Makefile resides) and type

$ make clean

If this doesn't help, then to additionally delete everything in the binary 
and library directories, go to <PREFIX>/bin and type

$ make distclean


-----------------------------------------------------------------------
- (2.1) Executing a single-column run:
-----------------------------------------------------------------------
   Before you can execute CLUBB, you must compile it (see Build section 
above.)
   <CLUBB DIR> refers the the directory where clubb is installed.

1. cd <CLUBB DIR>/input/case_setups

2. Edit <CASE>_model.in for each case you wish to run, or just leave 
   them as is.  Usually you will want to keep these the same.
   See the KK_microphys code for a description of Khairoutdinov and Kogan
   drizzle parameterization.
   See BUGSrad description below for a description of the interactive
   radiation scheme.
   Enabling radiation or microphysics parameterizations may increase runtime 
   considerably.

3. cd <CLUBB DIR>/input/stats
   Edit a stats file you would like to use.  A complete list of all computable
   statistics is found in all_stats.in.  Note that CLUBB now supports GrADS or
   netCDF, but you can only use the clubb_tuner using GrADS.

4. $ cd <CLUBB DIR>/input
   Edit tunable_parameters.in if you wish. The default values have been
   tested rigorously and will work with all the current cases.

5. $ cd <CLUBB DIR>/run_scripts
   $ ./run_scm.bash <MODEL CASE> or
   $ ./run_scm.bash <MODEL CASE> [PARAMETER FILE] [STATS FILE]

   Where the parameter file and stats file are optional arguments. The default
   is all_stats.in and tunable_parameters.in.

   The resulting data will be written in the directory "output"

-----------------------------------------------------------------------
- (2.2) Plotting output from a single-column run:
-----------------------------------------------------------------------

Plotting scripts are contained in the directory postprocessing/plotgen.

If you have MATLAB, you may use the script compare_plots_cases_driver.m.

If you have MATLAB and Linux, you can drive this MATLAB script
using the bash script plotgen.sh.

See postprocessing/plotgen/README for more information.

-----------------------------------------------------------------------
- (2.3) Executing a restart run:
-----------------------------------------------------------------------

After a long simulation has been performed, it is sometimes convenient to 
perform a new simulation that starts some time in the middle of the original 
simulation, rather than wasting time by starting again from the initial time.  
The new simulation is then called a "restart" simulation.

1.  Perform the original simulation of case <CASE> and save the GrADS output 
    files.  These data will be accessed to restart the simulation.

2.  Create a subdirectory in the CLUBB directory called "restart" and move the 
    GrADS output files to that subdirectory.

3.  Edit the following three variables at the end of the flag section of 
    the model file:

    l_restart = .true.
    restart_path = restart/<CASE>
    time_restart = initial time of restart run in seconds

    Compute time_restart as (time_initial + n_out * stats_tout), where n_out
    is the number of output intervals before the restart time.

4.  Execute the run as usual from /run_scripts using 
    
    ./run_scm.bash <CASE>

-----------------------------------------------------------------------
- (3.1) Executing a tuning run:
-----------------------------------------------------------------------

Do steps 1, 2, & 3 as outlined in the standalone run.

4. Edit input/tuner/error_<CASE>.in or select an existing one. Note that there 
   are two tuning subroutines, specified by tune_type in the error_<CASE>.in 
   /stats/ namelist.  
   If tune_type = 0, then the amoeba subroutine, a downhill simplex algorithm,
   will be used.  If runtype is any other value, then amebsa, a variant of 
   amoeba which uses simulated annealing instead, is used.  A complete 
   explanation of these minimization algorithms can be found 
   in _Numerical Recipes in Fortran 90_.
   Sometimes the variable names in CLUBB's zt and the LES grads files 
   will differ.  Currently, it is only possible to tune for variables that 
   occur in zt grid.

5. Edit run_tuner.bash to use your namelists

6. ./run_tuner.bash

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

1. Create and mount RAM disk on "output"

1. Run tuner

Linux Example
Note that you will need ram disk support compiled into your kernel, which is
typically the default on most systems.  Linux appears to be less flexible 
about when you are allowed to change the ramdisk size.

1. In grub.conf
   Append to 'kernel' line:
   kernel /vmlinuz-2.4.21-40.EL ro root=LABEL=/ ramdisk_size=262144

This sets ram disks to be 256 megabytes in size.  Note that your own system may
have other options besides the ramdisk_size.

2. $ mkfs.ext2 /dev/ram0

3. $ mount /dev/ram0 /home/dschanen/clubb/output

4. $ cd run_scripts

(Run your job)

Solaris Example
Note that these instructions only apply to Solaris 9 & 10

1. $ ramdiskadm -a clubb 256m 
Creates a virtual disk clubb that is 256 megabytes in size.

2. $ newfs /dev/ramdisk/clubb

3. $ mount /dev/ramdisk/clubb /home/dschanen/clubb/output/

4. $ cd run_scripts

(Run your job)

-----------------------------------------------------------------------
- (3.3) Executing an ensemble tuning run:
-----------------------------------------------------------------------

NOTES AND INSTRUCTIONS FOR THE CLUBB ENSEMBLE TUNER
-------------------------------------------------

Go to the main CLUBB directory and follow these instructions:

1)  Copy the ens_tune directory to a new directory with a slightly different
    name, such as ens_tune_two.  From this point on, I will refer to this new
    directory as ens_tune_xyz.

2)  Enter the new directory (ens_tune_xyz).

3)  Decide which cases will be tuned for (ex. DYCOMS2 RF02 DO, DYCOMS2 RF01,
    FIRE, ARM, BOMEX, etc.).

4)  For EACH case being tuned for:

     a)  Copy the <CASE>_model.in file from the input/case_setups directory 
	 to the ens_tune_xyz/ directory.  Make sure that the file is set up
         correctly.

     b)  Create/copy the <CASE>_stats_tune.in file.  A sample of this file can
         be found in the ens_tune/ directory (and subsequently in the
         ens_tune_xyz/ directory if it was copied from the ens_tune/ directory).
         This file is used in the tuning process itself.  In order for this
         process to work, stats_fmt in this file must be set to 'grads'.  The
         tuner code reads the LES GrADS files and compares the results with the
         CLUBB GrADS files (which are written according the specifications 
         stated in this file).  This file MUST also contain the names of the 
	 CLUBB variables that are being tuned for.  During the process of 
	 tuning, small GrADS files will be written that contain the results for 
         ONLY the variables that are listed here.  This makes the tuning go
         faster because it is not being slowed down by the writing of
         unnecessary variables.

     c)  Copy the <CASE>_stats.in file from the stats/ directory to the
         ens_tune_xyz/ directory.  Once the tuner has found the optimal value
         and the tuning process has finished, standalone CLUBB will run for
         each case that was tuned for with the values of the constants that the
         tuner found.  The statistics that it finally outputs will be directed
         according to this file.  This file produces the normal statistical
         files that one sees after a normal CLUBB run.  The main thing that 
	 might be changed is whether the final statistical output is in GrADS 
	 or in netCDF.  The user can choose either one.  The user just has to
         set stats_fmt to either 'grads' or 'netcdf'.  Of course, if the user
         desires, the sampling and output timesteps or the variables that are
         output can also be changed.

5)  Edit the error_messner_001.in file:

         The error_messner_001.in file is the same type of file as the error.in
         file in the regular CLUBB tuner.  This file must include all the 
         relevant information for EVERY case being tuned for, such as the 
	 CLUBB and LES stats files, the CLUBB run file, the vertical 
	 levels and time periods being tuned for, and the general weighting of 
         each case.  This file also must include other important factors such 
	 as the variable(s) being tuned for and the general weighting of each 
         variable.  Of course, the initial values of the constants and the 
	 amount of deviation allowed for each of the constants are also 
	 declared here.

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

         d) CLUBB:  The path to the CLUBB main level directory.  This is the
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
         (executable) files that you will need.  These files include clubb_tuner,
         clubb_standalone, and int2txt.  Finally, you need to know the size of
         the CLUBB output files for each of the cases you are tuning for.  These
         files will be generated when the individual tuning iterations are
         done.  One set of the CLUBB output files will be generated for each
         tuning iteration that you run.

         For example, let's say that the executable files sum to about 25 MB.
         You decide to tune for four cases -- ARM, BOMEX, DYCOMS2 RF02 DO, and
         DYCOMS2 RF01.  Let's say that the LES data files for these four cases
         sum to about 200 MB.  So, that's already 225 MB of space needed on
         EACH NODE used in the tuning run.  Now, let's say that the netCDF
         output files for the CLUBB results for these four cases sum to about
         150 MB.  However, you are running 25 iterations of the tuner on EACH
         NODE.  So, you will be producing those 150 MB worth of CLUBB files 25
         times over.  The total space needed on the node is then:

         25 MB + 200 MB + (150 MB/iteration)*(25 iterations) = 3975 MB

         So, about 4 GB of available disk space would be needed on each node
         in order to complete the run.  You would have to check the available
         disk space in the home directory on every node that you are using in
         order to make sure that EACH NODE has enough space.

         Say that you were running the above ensemble tuning run on nodes
         1-10.  That would yield a grand total of 250 iterations of the CLUBB
         tuner.  The first 10 runs would be launched (one on each node).  When
         ALL of them are completed, then the next set of 10 runs are launched,
         and so on and so on.  When the last set is finally complete, all the
         CLUBB results are copied back to Messner, and then deleted off the
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
         and then to see if they look good for all the cases we have in CLUBB.

-----------------------------------------------------------------------
- (5.1) Executing a Jacobian analysis:
-----------------------------------------------------------------------

1. $ cd ../run_scripts

2. Edit ../input_misc/jacobian.in. 

   Note that choosing a high delta_factor may make the model
   crash, which will result in no data (results for that term will come
   back as NaN).

3. $  ./run_jacobian.bash <CASE> [PARAMETER FILE] [STATS FILE]


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
  that we use to match CLUBB data.  Each output file includes an hour of 
  unphysical spinup time.  BOMEX is trade-wind cumulus; FIRE is marine 
  stratocumulus; ARM is continental cumulus; and ATEX is cumulus under 
  stratocumulus.  BOMEX, FIRE, and ATEX are statistically steady-state; 
  ARM varies over the course of a day.


Generated CLUBB GrADS files:
  bomex_zt.dat, fire_zt.dat, arm_zt.dat, atex_zt.dat, dycoms_zt.dat,
    wangara_zt.dat, <case>_zm, <case>_sfc ...
  These are the files generated by the hoc subroutine during a model run
  and compared to the COAMPS LES results.  The last of these generated
  when tuning will we the optimized results.  Every time the hoc subroutine
  is called, these are overwritten, so if you want to prevent them
  from being erased be sure to either copy the .ctl and .dat
  files to another directory or rename them.
  The tuner currently only uses variables in the zt file.

The namelist files:

  input/case_setups/bomex_model.in, fire_model.in, arm_model.in & atex_model.in.
  These files specify the standard CLUBB model parameters.  Usually these 
  do not need to be modified.

  input/stats/all_stats.in, nobudgets_stats.in, etc.
  These files specify statistics output for each simulation.  See
  all_stats.in for a complete list of the all output supported.
  
  input_misc/tuner/error_all.in, error_<CASE>.in, error_<DATE>.in.
  These specify tunable parameters, the initial spread of the simplex
  containing the tunable parameters, and which cases to "tune" for.

The randomization files:

  run_scripts/generate_seed.bash, input/tuner/rand_seed.dat, bin/int2txt
  The script uses intrinsic functionality in the Linux kernel to generate
  a pseudo random seed (the .dat) used by the tuner for randomizing initial
  parameters.  This works on any operating system with a Linux style 
  /dev/random (Solaris, Tru64, etc.) as well.  The seed file is now plain text 
  text and can be edited by hand.

------------------------------------------------------------------------
- (2.1) The BUGSrad Radiation scheme
------------------------------------------------------------------------

  This is an optional interactive radiation scheme, developed apart from
  CLUBB by Stephens, et al. The code used in CLUBB was obtained from Norm Wood 
  on 2004/07/10.
  When enabled, the analytic computation normally
  used for radiation is disabled.  BUGSrad is enabled in the 
  input/case_setups/<CASE>_model.in file by setting lbugsrad = .true.
  Furthermore, you must compile CLUBB with the -Dradoffline preprocessor flag.

  BUGSrad allows the output of the following variables:

  Momentum grid:
  Frad, Frad_SW, Frad_LW:  Radiative Flux; Short-wave/Long-wave component;

  Thermodynamic grid:
  radht, radht_SW, radht_LW:  Radiative Heat; Short-wave/Long-wave component;

  The thlm_forcing variable will also have radht added to it.  This is an
  explicit contribution to the thlm calculation.

  Note that for most cases SW and LW components are not calculated 
  without using BUGSrad.

------------------------------------------------------------------------
- (2.2) The COAMPS microphysics scheme
------------------------------------------------------------------------

     COAMPS microphysics is a single-moment scheme that includes
  the following hydrometeor categories: cloud water, rain, cloud ice,
  snow, and graupel.  It is based on Rutledge and Hobbs (1983).

     COAMPS was developed by the Naval Research Laboratory,
  Monterey, California.  COAMPS is a registered trademark of the
  Naval Research Laboratory.

     COAMPS microphysics is not distributed with the code outside of UWM
  because of licensing restrictions.

------------------------------------------------------------------------
- (2.3) The Morrison microphysics scheme
------------------------------------------------------------------------
  Morrison microphysics is a double-moment scheme that can predict mixing
  ratio's and number concentrations for cloud water, rain, cloud ice,
  snow, and graupel.  Details of its implementation may be found in:

  H. Morrison, J. A. Curry, and V. I. Khvorostyanov, 2005: A new double-
  moment microphysics scheme for application in cloud and climate models. 
  Part 1: Description. J. Atmos. Sci., 62, 1665–1677. 

------------------------------------------------------------------------
- (3.1) The passive scalar code
------------------------------------------------------------------------

The CLUBB code can be run with additional non-interactive scalars.
The scalars in the code provide a generalized way of simulating a 
passive scalar in the atmosphere (e.g. carbon dioxide) 

By default CLUBB is configured to run without any passive scalars.  To use 
this option, you must modify the input/case_setups/<CASE>_model.in
so that sclr_dim is equal to the number of passive scalars and also set 
ii<SCALAR NAME> to point the index in the array containing the passive scalar.
The intial sounding can be done at run time in this file, but large scale
forcing and surface fluxes for these passive scalars must be configured
in the clubb_driver code and handled at compile time.
To output the scalar fields from a CLUBB simulation, be sure to include 
sclram, sclrap2, etc. your stats file. See input/stats/all_stats.in for a
complete list (commented out by default).

Currently the code contains eddy-diffusivity scalar code (edsclram, edsclrbm)
and the more sophisticated high-order scalars (sclram, sclrbm).  Both use two
dimensional arrays, but the code and results for each is separated.
The high-order scalars require a tolerance, called `sclrtol' to set in the
namelist file. Generally this value should be larger than machine epsilon.

Initially, the scalar arrays were configured to contain two vertical columns
containing copies of thl and rt.
The code is sufficiently general that an arbitrary number of scalars can be
added with a small number of modifications.
The following files must be adjusted to customize the scalars:

stats_variables.F:  Currently this is only configured to generate data for 
2 elements for each of the sclr arrays.  Using the sclra and sclrb variables
as a template, it should be possible to add additional variables (sclrc, etc.)
without any issues.

The Namelists:

Within the existing <CASE>_model.in for each run, a sounding for the scalar variable
must be added. 

Setting sclrm(:,1) equal to thlm, and sclrm(:,2) equal to rtm in 
the BOMEX case can be done like so:

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

------------------------------------------------------------------------
- (4.1) Contributing code changes 
------------------------------------------------------------------------

If you have changes that you'd like to see included in the repository version
of CLUBB, please feel free to contribute them.  We'll gladly review your changes 
and consider them for inclusion.  

To contribute the changes, type 'svn update' to merge the latest repository
version of CLUBB with your local version, and then create a file containing 
the differences between the repository and your version:

$ svn diff > patchfile

Then email us patchfile.  See the svn book online for more details about svn.

****************************************************************************
