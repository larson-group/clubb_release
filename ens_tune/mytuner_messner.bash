#!/bin/bash
# $Id: mytuner_messner.bash,v 1.4 2007-05-04 18:08:05 dschanen Exp $
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

# Note: ITERMAX needs to be a multiple of your total nodes!
ITERMAX=6
NODES=( tom15 tom16 tom17 )
HOC=$HOME/current/hoc_v2.2_tuner
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
	local iter=$2  # 2nd argument

	# Create directory & copy files
	rsh $tom "mkdir -p $HOC/ens_tune"
	$RCP error*.in "$tom:$HOC/ens_tune/error.in"
	$RCP *_model.in *_stat*.in "$tom:$HOC/ens_tune"

	# Create random seed data
	rsh $tom "$HOC/bin/int2txt /dev/random > $HOC/ens_tune/rand_seed.dat"

	# Create hoc.in file
	for EXP in "${EXPERIMENTS[@]}"; do
		rsh $tom "cat $HOC/ens_tune/${EXP}_model.in $HOC/ens_tune/${EXP}_stats_tune.in > $HOC/ens_tune/${EXP}_hoc.in"
	done

	# Execute tuning run
	echo Tuning iteration $iter on $tom
	rsh $tom "cd $HOC/ens_tune && echo no | ../bin/hoc_tuner &> tune.log"

	# Delete GrADS files from tuning run
	rsh $tom "rm -f $HOC/ens_tune/*_zt.??? $HOC/ens_tune/*_zm.??? $HOC/ens_tune/*_sfc.???"

	# Run standalone HOC
	for EXP in "${EXPERIMENTS[@]}"; do
		rsh $tom "cat $HOC/ens_tune/${EXP}_model.in $HOC/ens_tune/${EXP}_stats.in > $HOC/standalone/hoc.in"
		rsh $tom "cat $HOC/standalone/standalone_*.in > $HOC/standalone/standalone.in"
		rsh $tom "cd $HOC/standalone && ../bin/hoc_standalone >> ../ens_tune/tune.log"
	done

	# Move results data
	rsh $tom "mv $HOC/ens_tune $HOC/ens_tune_$iter";
	rsh $tom "mv $HOC/standalone/standalone_*.in  $HOC/ens_tune_$iter";
	# for GrADS output
	#rsh $tom "mv $HOC/standalone/*_*.???  $HOC/ens_tune_$iter";
	# for NetCDF output
	rsh $tom "mv $HOC/standalone/*_*.nc  $HOC/ens_tune_$iter";
}

################################################################################
# Main code
################################################################################

for node in "${NODES[@]}"; do 
	echo "Copying files to $node"
 	rsh $node "mkdir -p $HOC/les_data"
	rsh $node "mkdir -p $HOC/bin"
	rsh $node "mkdir -p $HOC/standalone"

	# Copy hoc_tuner executable and ens_tune namelists
	$RCP $HOC"/bin/hoc_tuner" $node":$HOC/bin"
	$RCP $HOC"/bin/hoc_standalone" $node":$HOC/bin"
	$RCP $HOC"/bin/int2txt" $node":$HOC/bin"

	for EXP in "${EXPERIMENTS[@]}"; do
		# Copy GrADS data for LES simulations to tom node
 		$RCP $HOC/les_data/${EXP}*_sm.??? $node:$HOC/les_data/ ;
	done

	echo "Copy complete";
done

# Initial iteraton
ITER=1
# Loop over the number of times, incrementing by the total nodes
while [  $ITER -lt $ITERMAX ]; do

	for node in "${NODES[@]}"; do 
		tune $node $ITER &
		let ITER=ITER+1
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

