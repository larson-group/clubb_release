#!/bin/bash
# $Id: mytuner_messner.bash,v 1.5 2007-06-28 16:38:13 griffinb Exp $
################################################################################
#
# Description: This is a port of Chris Golaz's ensemble tuning script to the
#   Messner Beowulf cluster.  The functionality is similar but not indentical
#   to the original.  The basic structure is as follows:
#
#  < Messner > ---- +--- Node1 tuning run over variables/cases (rand # 1)
#  Messner holds    |
#  the results, the +--- Node2 tuning run over variables/cases (rand # 2)
#  namelists, and   .
#  is the `control' .
#  node.            +--- NodeN tuning run over variables/cases (rand # N)
#
################################################################################

################################################################################
# Global variables
################################################################################

# Note: MEMBERMAX needs to be a multiple of your total nodes!
MEMBERMAX=6
NODES=( tom15 tom16 tom17 )
HOC=$HOME/clubb
CASE=combined_001
ARCHIVE=$HOME/archive/ensemble
EXPERIMENTS=( bomex dycoms2_rf01 )

# Kluge to get rid of the annoying kerberos message on Redhat Ent. 3
RCP='/usr/bin/rcp'

################################################################################
# Function to copy files from control node and run the hoc_tuner
################################################################################
function tune ( ) {
	local tom=$1   # 1st argument
	local member=$2  # 2nd argument

	# Create directory & copy files
	rsh $tom "mkdir -p $HOC/ens_tune"

	# Create random seed data
	rsh $tom "$HOC/install/latest/int2txt /dev/random > $HOC/ens_tune/rand_seed.dat"

	# Each remote checkout resolves model JSON and legacy tuner controls at launch.
	rsh $tom "python3 $HOC/run_scripts/run_tuner.py -prepare_only -controls $CASE -run_type multiple -model_mult ${EXPERIMENTS[*]} -output_dir $HOC/ens_tune"

	# Execute tuning run
	echo Tuning ensemble member $member on $tom
	rsh $tom "cd $HOC/ens_tune && echo no | ../install/latest/clubb_tuner &> tune.log"

	# Delete GrADS files from tuning run
	rsh $tom "rm -f $HOC/ens_tune/*_zt.??? $HOC/ens_tune/*_zm.??? $HOC/ens_tune/*_sfc.???"

	for EXP in "${EXPERIMENTS[@]}"; do
        rsh $tom "python3 $HOC/run_scripts/run_scm.py -params \$(ls -t $HOC/input/tunable_parameters_*.json | head -1) $EXP >> $HOC/ens_tune/tune.log"
	done

	# Move results data
	rsh $tom "mv $HOC/ens_tune $HOC/ens_tune_$member";
	rsh $tom "mv $HOC/input/tunable_parameters_*.json  $HOC/ens_tune_$member";
	# for GrADS output
	#rsh $tom "mv $HOC/standalone/*_*.???  $HOC/ens_tune_$member";
	# for NetCDF output
	rsh $tom "mv $HOC/output/*_*.nc  $HOC/ens_tune_$member";
}

################################################################################
# Main code
################################################################################

for node in "${NODES[@]}"; do 
	echo "Copying files to $node"
 	rsh $node "mkdir -p $HOC/les_data"
	rsh $node "mkdir -p $HOC/install/latest"
	rsh $node "mkdir -p $HOC/standalone"

	# Copy current CLUBB executables; JSON catalogs come from the remote checkout
	$RCP $HOC"/install/latest/clubb_tuner" $node":$HOC/install/latest"
	$RCP $HOC"/install/latest/clubb_standalone" $node":$HOC/install/latest"
	$RCP $HOC"/install/latest/int2txt" $node":$HOC/install/latest"

	for EXP in "${EXPERIMENTS[@]}"; do
		# Copy GrADS data for LES simulations to tom node
 		$RCP $HOC/les_data/${EXP}*_sm.??? $node:$HOC/les_data/ ;
	done

	echo "Copy complete";
done

# Initial ensemble member
MEMBER=1
# Loop over the number of times, incrementing by the total nodes
while [  $MEMBER -lt $MEMBERMAX ]; do

	for node in "${NODES[@]}"; do 
		tune $node $MEMBER &
		let MEMBER=MEMBER+1
	done
	wait

done

# Check that the archive directory exists
if [ ! -e $ARCHIVE ]; then
	# if not, create it
	mkdir -p $ARCHIVE/$CASE
	cp sortresults.bash $ARCHIVE/$CASE
#	cp average.bash $ARCHIVE
fi

# Copy back and clean up
for node in "${NODES[@]}"; do 
	mkdir $ARCHIVE/$CASE
	# Copy back results
	$RCP -r $node:$HOC/ens_tune_* $ARCHIVE/$CASE/
	# Remove files from nodes
	rsh $node "rm -rf $HOC"
done

