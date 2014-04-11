#!/bin/bash
# $Id: mytuner_sequential.bash 5071 2011-03-30 19:20:26Z roehl@uwm.edu $
################################################################################
#
# Description: 
#   This is a sequential version of Chris Golaz's ensemble tuning script,
#   which was created after the Messner cluster was dismantled.
#   Each ensemable member is run in sequence.
#   It will move the first foot to match the previous member if that member
#   tuned with a lower cost function than the previous.
################################################################################

################################################################################
# Global variables
################################################################################

MEMBERMAX=2 # Total number of ensemble members
CLUBB=$HOME/clubb
ENSEMBLE_DIR=`dirname $0` # Directory with this script
CASE=combined_001 # This plus _error.in is the name of tuner namelist file.
ARCHIVE=$HOME/archive/ensemble # Directory for the output
EXPERIMENTS=( bomex dycoms2_rf01 ) # Add additional simulations here

################################################################################
# Function to copy files and run the clubb_tuner
################################################################################
function tune ( ) {
	local member=$1

	# Create random seed data
	$CLUBB/bin/int2txt /dev/urandom > $CLUBB/ens_tune/rand_seed.dat

	# Create model namelist files
	for EXP in "${EXPERIMENTS[@]}"; do
		cat $CLUBB/input/case_setups/${EXP}_model.in | sed 's/debug_level\s*=\s*.*/debug_level = 0/g' > $CLUBB/ens_tune/${EXP}_hoc.in
		cat $CLUBB/ens_tune/${EXP}_stats_tune.in >> $CLUBB/ens_tune/${EXP}_hoc.in
	done
	# Execute tuning run
	echo Tuning ensemble member $member on `hostname` at `date +"%R %D"` | tee $CLUBB/ens_tune/tuner.log 
	cd $ENSEMBLE_DIR && echo no | ../bin/clubb_tuner &> tune.log

	# Delete GrADS files from tuning run
	rm -f $CLUBB/output/*_zt.??? $CLUBB/output/*_zm.??? $CLUBB/output/*_sfc.???

	# Run standalone CLUBB with the new parameter values
	for EXP in "${EXPERIMENTS[@]}"; do
		cat $CLUBB/input/case_setups/${EXP}_model.in $CLUBB/input/stats/standard_stats.in > clubb.in
		cat $CLUBB/input/tunable_parameters/tunable_parameters_*.in >> clubb.in
		../bin/clubb_standalone >> $ENSEMBLE_DIR/tune.log
	done

	# Move results data
	mkdir -p $ARCHIVE/$CASE/ens_tune_$member

	# Move foot if current cost is less than previous cost
        currCost=$(grep -i '\$' tune.log | tr -d '$''[:blank:]')
        if [ `echo "$currCost < $prevCost" |bc ` == 1 ] ; then
                PREVTUNER=`ls $CLUBB/input/tunable_parameters/tunable_parameters_*.in`
                cat "error_"$CASE".in" $PREVTUNER &> "error.in"
                echo -en '\E[47;34m'"\033[1mRecreating error.in with new parameter values\033[0m\n"
                prevCost=$currCost
                counter=0
        fi
        counter=$(expr $counter + 1)
        if [ $counter -gt 15 ]; then
                counter=0
                prevCost=1
                cat "error_"$CASE".in" $CLUBB/input/tunable_parameters/tunable_parameters.in &> "error.in"
                echo -en '\E[47;34m'"\033[1mReverting to original parameter values\033[0m\n"

        fi

	mv $CLUBB/input/tunable_parameters/tunable_parameters_*.in  $ARCHIVE/$CASE/ens_tune_$member
	mv $ENSEMBLE_DIR/tune.log  $ARCHIVE/$CASE/ens_tune_$member
	# for GrADS output
	mv $CLUBB/output/*_*.???  $ARCHIVE/$CASE/ens_tune_$member
	# for NetCDF output
	#mv $CLUBB/output/*_*.nc  $CLUBB/ens_tune_$member
}

################################################################################
# Main code
################################################################################

# Check that the archive directory exists
if [ ! -e $ARCHIVE/$CASE ]; then
	# if not, create it
	mkdir -p $ARCHIVE/$CASE
	cp sortresults.bash $ARCHIVE/$CASE
	cp analyze_results.py $ARCHIVE/$CASE
#	cp average.bash $ARCHIVE
fi

# Create tuner namelist
cat "error_"$CASE".in" $CLUBB/input/tunable_parameters/tunable_parameters.in > "error.in"

cd $ENSEMBLE_DIR

# Setup variables
counter=0
prevCost=1.

#Make the les_data folder if it doesn't exsist
mkdir $CLUBB/les_data &> $CLUBB/ens_tune/tuner.log

# Download LES data for experiments
for EXP in "${EXPERIMENTS[@]}"; do
	if [ ! -e "$CLUBB/les_data/${EXP}_coamps_sm.ctl" ]; then
		svn export http://carson.math.uwm.edu/repos/clubb_benchmark_runs/les_runs/${EXP}_coamps_sm.ctl $CLUBB/les_data/${EXP}_coamps_sm.ctl
		svn export http://carson.math.uwm.edu/repos/clubb_benchmark_runs/les_runs/${EXP}_coamps_sm.dat $CLUBB/les_data/${EXP}_coamps_sm.dat
	fi
done

# Loop over the number of ensemble points, incrementing by 1
for (( MEMBER=1; MEMBER <= $MEMBERMAX; MEMBER++ )); do 
	tune $MEMBER
done

