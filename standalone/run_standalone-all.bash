#######################################################################
#
# Script to run the standalone hoc program.  
# Tested with BASH.  Not tested with Korn shell or Bourne(sh) shell.
# Edit to change run
#
#######################################################################
# Useful on multiprocessor machines with OpenMP capable Fortran
# OMP_NUM_THREADS=2
#######################################################################
# This will loop over all runs in sequence 
for (( x = 1; x <= 10; x++)); do
  case $x in
   1 )
     RUN_CASE=arm ;;
   2 )
     RUN_CASE=atex ;;
   3 )
     RUN_CASE=bomex ;;
   4 )
     RUN_CASE=dycoms2_rf01 ;;
   5 )
     RUN_CASE=dycoms2_rf02_d ;;
   6 )
     RUN_CASE=dycoms2_rf02_nd ;;
   7 )
     RUN_CASE=fire ;;
   8 )
     RUN_CASE=nov11_altocu ;;
   9 )
     RUN_CASE=wangara ;; 
   10)
     RUN_CASE=20050821_1218 ;;
   esac
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

done
