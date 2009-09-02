$Id: Readme.txt,v 1.2 2007-05-04 18:08:05 dschanen Exp $

These scripts differ from the regular tuning run in that they can be used to
create many different sets of constant parameters based on the random starting
simplexes generated with differing seeds.  Chris Golaz used an ensemble 
method to create the plots that appear in the tuner paper.  These scripts are
based on the scripts he used, but differ in that they have been setup
machines UW-Milwaukee.
The mytuner_messner.bash was configured for use on the nodes of the 
Messner Beowulf cluster at UWM in approximately 2007.
The mytuner_sequential.bash was configured for use on a single node in
Septemeber 2009.

This directory is more or less a template for how to setup a new case, but
you will need to alter all these scripts and namelists to suit your ensemble 
tuning needs.
You will also need to copy over and modify stats.in files to include only those
variables for which wish to tune for.

The mytuner_messner.bash script was developed in 2007 and no longer reflects
the current CLUBB directory structure, though it could be used as a very
rough template of how to add a new script for a cluster environment.

Doing a basic `at' job:

1.) Edit the first few lines of mytuner_messner.bash, setting the variables
    MEMBERMAX, CLUBB, CASE, ARCHIVE, and EXPERIMENTS as needed.

    MEMBERMAX is the total number of times the tuner will be run.
    CLUBB is location of your clubb installation.
    ARCHIVE is the where the tunable parameters and GrADS data will be copied
    to, along with the script to sort the results.
    EXPERIMENTS is any array containing all the cases you wish to run.

    In ens_tuen you must have:
    error_<CASE>.in containing the namelists for the tuner, setup as the files
    in input_misc/tuner/
    <EXPERIMENT ...>_stats_tune.in which you will only want to include
    variables you're tuning for.

    You will probably also want to do typical tuner simplifications like
    setting the debug_level to 0.
    
2.) at -f ./atjob.bash now

(Come back in a number of hours/days)

3.) Go to your ARCHIVE/CASE directory.  Run sortresults.bash for sorted results.
 
