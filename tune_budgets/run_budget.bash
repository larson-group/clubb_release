#######################################################################
#
# Tested with BASH.  Not tested with Korn shell or Bourne(sh) shell.
# Edit to change run
#
#######################################################################
# Useful on compilers which have OpenMP
# OMP_NUM_THREADS=2
#######################################################################
# Select a run, comment out the rest
# RUN_CASE=atex 
# RUN_CASE=arm 
# RUN_CASE=bomex 
# RUN_CASE=dycoms2_rf01 
# RUN_CASE=dycoms2_rf02_do
# RUN_CASE=dycoms2_rf02_ds
# RUN_CASE=dycoms2_rf02_nd
  RUN_CASE=fire 
# RUN_CASE=nov11_altocu
#######################################################################
# Check for necessary namelists.  If files exist, then
# copy them over to the general input files.

 BUDGET_IN=$RUN_CASE'_budget.in'

 if [ ! -e "$BUDGET_IN" ] ; then
	echo $BUDGET_IN " does not exist"
	exit 1
 fi

 if [ -e 'budget.in' ] ; then
	rm -f 'budget.in'
 fi

 ln -s $BUDGET_IN 'budget.in'


#######################################################################
#
# State which case is being run
 echo "Running" $RUN_CASE
# Run HOC
./hoc_tuner_budget_terms

# remove the temporary error.in file
 rm -f 'budget.in'
