#######################################################################
#
# Script to run the standalone hoc program.  
# Tested with BASH.  Not tested with Korn shell or Bourne(sh) shell.
# Edit to change run
#
#######################################################################
# Useful on compilers which have OpenMP
# OMP_NUM_THREADS=2
#######################################################################
# Select a run, comment out the rest
# RUN_CASE=atex 
  RUN_CASE=arm 
# RUN_CASE=bomex 
# RUN_CASE=dycoms2_rf01 
# RUN_CASE=dycoms2_rf02_d 
# RUN_CASE=dycoms2_rf02_nd 
# RUN_CASE=fire 
# RUN_CASE=nov11_altocu
# RUN_CASE=wangara  
# RUN_CASE=20050821_1218
#######################################################################
# Check for necessary namelists.  If files exist, then
# copy them over to the general input files.

 STANDALONE_IN='standalone_'$RUN_CASE'.in'

 if [ ! -e "$STANDALONE_IN" ] ; then
	echo $STANDALONE_IN " does not exist"
	exit 1
 fi

 if [ -e 'standalone.in' ] ; then
	rm -f 'standalone.in'
 fi

 ln -s $STANDALONE_IN 'standalone.in'


#######################################################################
#
# State which case is being run
 echo "Running" $RUN_CASE
# Run HOC
./hoc_standalone

# remove the temporary error.in file
 rm -f 'standalone.in'
