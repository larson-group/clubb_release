$Id: Readme.txt,v 1.2 2007-05-04 18:08:05 dschanen Exp $

These scripts differ from the regular tuning run in that they can be used to
create many different sets of constant parameters based on the random starting
simplexes generated with differing seeds.  Chris Golaz used an ensemble 
method to create the plots that appear in the tuner paper.  These scripts are
based on the scripts he used, but differ in that they have been setup for use
on the nodes of the Messner Beowulf cluster at UWM.

This directory is more or less a template for how to setup a new case, but
you will need to alter all these scripts and namelists to suit your ensemble 
tuning needs.
You will also need to copy over model.in and stats.in files and modify those
for what you wish to tune for.
The mytuner_messner.bash script has been tested, but may not be robust enough
for every possible thing that could go wrong.  Make sure there is enough disk
space, etc. before you run the script.

Doing a basic `at' job:

1.) Edit the first few lines of mytuner_messner.bash, setting the variables
    ITERMAX, NODES, HOC, CASE, ARCHIVE, and EXPERIMENTS as needed.

    Optionally,  Edit: 
	# for GrADS output
	#rsh $tom "mv $HOC/standalone/*_*.???  $HOC/ens_tune_$iter";
	# for NetCDF output
	rsh $tom "mv $HOC/standalone/*_*.nc  $HOC/ens_tune_$iter";

    (Commenting out one, commenting in the other).

2.) at -f ./atjob.bash now

(Come back in a number of hours/days)

3.) Go to your ARCHIVE/CASE directory.  Run sortresults.bash for sorted results.
 
