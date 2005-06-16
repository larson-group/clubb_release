#######################################################################
#
# Script to run the standalone hoc program.  
# Tested with BASH.  Not tested with Korn shell or Bourne(sh) shell.
# Edit to change run
#
#######################################################################
# Useful on multiprocessor machines with OpenMP capable fortran
OMP_NUM_THREADS=2
#######################################################################
# Select a run, comment out the rest
#   Or alternatively, uncomment do loop to do all cases
for (( x = 1; x <= 9; x++)); do
  case $x in
   1 )
     RUN_CASE=atex ;;
   2 )
     RUN_CASE=arm ;;
   3 )
     RUN_CASE=bomex ;;
   4 )
     RUN_CASE=dycoms2_rf01 ;;
   5 )
     RUN_CASE=fire ;;
   6 )
     RUN_CASE=wangara ;; 
   7 )
     RUN_CASE=dycoms2_rf02_d ;;
   8 )
     RUN_CASE=dycoms2_rf02_nd ;;
   9 )
     RUN_CASE=nov11_altocu ;;
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
# Run HOC
 ./hoc_standalone 

# remove the temporary error.in file
 rm -f 'standalone.in'

done
