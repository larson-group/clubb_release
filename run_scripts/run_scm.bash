#!/bin/bash
###############################################################################
# $Id$
#
# Desciption:
#   Script to run the standalone CLUBB program.
#   Tested with GNU Bash v2 &3. Might work with Ksh.
###############################################################################

# Useful variable on multiprocessor machines with OpenMP capable Fortran.
# Uncomment and set to use. The clubb_standalone program is not parallel,
# but many LAPACK/BLAS libraries are.
#
# export OMP_NUM_THREADS=2

# Fields for parameters passed in
NIGHTLY=false
TIMESTEP_TEST=false
ZT_GRID=false
ZM_GRID=false
PERFORMANCE_TEST=false
OUTPUT_DIR="/home/`whoami`/nightly_tests/output"
NAMELISTS="clubb.in"

# Figure out the directory where the script is located
scriptPath=`dirname $0`

# Store the current directory location so it can be restored
restoreDir=`pwd`

# Change directories to the one the script is located in
cd $scriptPath

run_case()
{
    # Enable G95 runtime option that sets uninitialized memory to a NaN value
    G95_MEM_INIT="NAN"
    export G95_MEM_INIT

    echo "Running $run_case"
    
    # Run the CLUBB model
    ../bin/clubb_standalone

    RESULT=$?

    # Remove the namelists
    rm -f $NAMELISTS
}

# Note that we use `"$@"' to let each command-line parameter expand to a 
# separate word. The quotes around `$@' are essential!
# We need TEMP as the `eval set --' would nuke the return value of getopt.
TEMP=`getopt -o z:m:l:t:s:p:nhe --long zt_grid:,zm_grid:,levels:,timestep_test:,stats:,parameter_file:,performance_test,nightly,help \
     -n 'run_scm.bash' -- "$@"`

if [ $? != 0 ] ; then echo "Terminating..." >&2 ; exit 1 ; fi

# Note the quotes around `$TEMP': they are essential!
eval set -- "$TEMP"

while true ; do
	case "$1" in
		-z|--zt_grid) # Set the zt grid
            ZT_GRID=true
            grid_path=$2
            
            # Make sure the file exists
            if [ ! -f $grid_path ];
            then
                echo "The grid file specified could not be found!"
                exit 1;
            else
                grid_path=`echo $grid_path | sed 's/\//\\\\\//g'`
            fi

            shift 2 ;;
		-m|--zm_grid) # Set the zm grid
            ZM_GRID=true
            grid_path=$2

            if [ ! -f $grid_path ];
            then
                echo "The grid file specified could not be found!"
                exit 1;
            else
                grid_path=`echo $grid_path | sed 's/\//\\\\\//g'`
            fi

            shift 2 ;;
        -l|--levels) # Set the number of levels. This is used with the grids specified with -z and -m
            grid_nz=$2

            # Check to make sure this is an integer
            if [ $grid_nz -ne $grid_nz 2> /dev/null ];
            then
                echo "The number of levels specified are invalid!"
                exit 1;
            fi
            shift 2 ;;
        -t|--timestep_test)
            TIMESTEP_TEST=true
            test_ts=$2

            # Check to make sure this is an integer
            if [ $test_ts -ne $test_ts 2> /dev/null ];
            then
                echo "The timestep specified is invalid!"
                exit 1;
            fi

            shift 2 ;;
        -s|--stats) 
            stats_file=$2

            if [ ! -f $stats_file ];
            then
                echo "The stats file does not exist!"
                exit 1;
            fi

            shift 2 ;;
        -p|--parameter_file) # Set the parameter file
            parameter_file=$2

            if [ ! -f $parameter_file ]
            then
                echo "The parameter file does not exist!"
                exit 1;
            fi

            shift 2 ;;
        -n|--nightly) 
            NIGHTLY=true

            shift;;
        -e|--performance_test)
            PERFORMANCE_TEST=true

            shift;;
        -h|--help) # Print the help message
            echo -e "Usage: run_scm.bash [OPTION]... case_name"
            echo -e "\t-z, --zt_grid=FILE\t\tThe path to the zt grid file"
            echo -e "\t-m, --zm_grid=FILE\t\tThe path to the zm grid file"
            echo -e "\t-l, --levels=NUM\t\tThe number of levels"
            echo -e "\t-t, --timestep_test=LENGTH\tRun the case at specified timestep"
            echo -e "\t-s, --stats=FILE\t\tSpecify the stats file to use"
            echo -e "\t-p, --parameter_file\t\tSet the parameter file to use"
            echo -e "\t-r, --performance_test\t\tDisable statistics output and set debug"
            echo -e "\t\t\t\t\tlevel to 0 for performance testing"
            echo -e "\t-h, --help\t\t\tPrints this help message"

            exit 1 ;;
		--) shift ; break ;;
		*) echo "Something bad happened!" ; exit 1 ;;
	esac
done

# Make sure we are only doing either zt grid tests or zm grid tests
if [ $ZT_GRID == true ] && [ $ZM_GRID == true ];
then
    echo "Cannot specify both a ZT grid and a ZM grid."
    exit 1
fi

# Make sure that the number of grid levels was specified
if [ $ZT_GRID == true ] || [ $ZM_GRID == true ];
then
    if [ -z $grid_nz ];
    then
        echo "You must specify the number of levels."
        exit 1
    fi
fi

model_file='../input/case_setups/'${!#}'_model.in'
run_case=${!#}

# Check to see if the model file exists
if [ ! -e "$model_file" ];
then
    echo "$model_file does not exist"
    exit 1
fi

# Set defaults if they were not passed in
if [ -z $parameter_file ];
then
    parameter_file="../input/tunable_parameters/tunable_parameters.in"
fi

if [ -z $stats_file ];
then
    stats_file="../input/stats/standard_stats.in"
fi

if [ $NIGHTLY == true ]; 
then
#    if [ "$run_case" = gabls2 || "$run_case" = cobra ]; 
    if [ "$run_case" = gabls2 ]; 
    then
        # GABLS2 uses scalars
        stats_file='../input/stats/nightly_stats_scalars.in'
    elif [ "$run_case" = cobra ]; 
        then
           # Cobra uses scalars
           stats_file='../input/stats/nightly_stats_scalars.in'
    else
        stats_file='../input/stats/nightly_stats.in'
        #stats_file='../stats/nobudgets_stats.in' 
    fi
fi

if [ $NIGHTLY == true ]; 
then
    cat $parameter_file > $NAMELISTS
    # This is needed because the model file now contains stats_tout.
    # Here we replace the repository version of stats_tout with an hour output.
    # We then save the zt and zm files from this run for the profile plots because
    # this saves disk space, and the profile plots are the same regardless if
    # stats_tout = 1 hour, 5 minutes, or 1 minute.
    # The regular expression used here matches:
    # 'stats_tout' (0 or > whitespaces) '=' (0 or > whitespaces) (0 or > characters)
    # and replaces it with 'stats_tout = 3600.'
    cat $model_file | sed 's/stats_tout\s*=\s*.*/stats_tout = 3600\./g' >> $NAMELISTS
    cat $stats_file >> $NAMELISTS
    run_case

    # Move the ZT and ZM files out of the way
    if [ "$RESULT" != 0 ]; then
        rm "../output/$run_case"_zt.ctl
        rm "../output/$run_case"_zt.dat
        rm "../output/$run_case"_zm.ctl
        rm "../output/$run_case"_zm.dat
        rm "../output/$run_case"_sfc.ctl
        rm "../output/$run_case"_sfc.dat
    else
        mv "../output/$run_case"_zt.ctl "$OUTPUT_DIR"/CLUBB_current/
        mv "../output/$run_case"_zt.dat "$OUTPUT_DIR"/CLUBB_current/
        mv "../output/$run_case"_zm.ctl "$OUTPUT_DIR"/CLUBB_current/
        mv "../output/$run_case"_zm.dat "$OUTPUT_DIR"/CLUBB_current/
        case $run_case in
            # We only run TWP_ICE and Cloud Feedback once so we want to keep the SFC files.
            # The other cases are rerun with stats_tout = 1 minute (or 5 minutes for RICO),
            # and sfc files from the second run are used in the nightly plots for these cases.
            twp_ice | cloud_feedback* )
                mv "../output/$run_case"_sfc.ctl "$OUTPUT_DIR"/CLUBB_current/
                mv "../output/$run_case"_sfc.dat "$OUTPUT_DIR"/CLUBB_current/
                ;;
            * )
                rm "../output/$run_case"_sfc.ctl
                rm "../output/$run_case"_sfc.dat
                ;;
        esac
    fi
		
    # Run again with a finer output time interval for the sfc.  Even though
    # running the cases twice takes longer, we do this to save disk space
    # (the zt and zm files are much smaller when stats_tout = 1 hour).
    # Note, we do not run TWP_ICE and Cloud Feedback a second time because
    # they are 25- and 30-day simulations.
    case $run_case in 
        twp_ice | cloud_feedback* )
            ;;
        * )
            case $run_case in 
                rico | astex_a209 )
                    # This was added because RICO uses a 300 s timestep
                    # and cannot be run with stats_tout = 60.
                    cat $parameter_file > $NAMELISTS
                    cat $model_file | sed 's/stats_tout\s*=\s*.*/stats_tout = 300\./g' >> $NAMELISTS
                    cat $stats_file >> $NAMELISTS
                    ;;
                * )
                    cat $parameter_file > $NAMELISTS
                    cat $model_file | sed 's/stats_tout\s*=\s*.*/stats_tout = 60\./g' >> $NAMELISTS
                    cat $stats_file >> $NAMELISTS
                    ;;
            esac

            run_case
	
            #Now move the SFC file
            if [ "$RESULT" != 0 ]; then
                rm "../output/$run_case"_zt.ctl
                rm "../output/$run_case"_zt.dat
                rm "../output/$run_case"_zm.ctl
                rm "../output/$run_case"_zm.dat
                rm "../output/$run_case"_sfc.ctl
                rm "../output/$run_case"_sfc.dat
            else
                rm "../output/$run_case"_zt.ctl
                rm "../output/$run_case"_zt.dat
                rm "../output/$run_case"_zm.ctl
                rm "../output/$run_case"_zm.dat
                mv "../output/$run_case"_sfc.ctl "$OUTPUT_DIR"/CLUBB_current/
                mv "../output/$run_case"_sfc.dat "$OUTPUT_DIR"/CLUBB_current/
            fi
            ;;
    esac
elif [ $TIMESTEP_TEST == true ]; 
then
    # Set the model timestep for all cases (and the stats output timestep
    # unless l_stats is overwritten to .false.) to timestep test_ts.
    cat $parameter_file > $NAMELISTS
    # Use this version if statistical output is desired.
    #cat $MODEL_IN | sed -e 's/dtmain\s*=\s*.*/dtmain = '$test_ts'/g' \
    #                    -e 's/dtclosure\s*=\s*.*/dtclosure = '$test_ts'/g' \
    #                    -e 's/stats_tsamp\s*=\s*.*/stats_tsamp = '$test_ts'/g' \
    #                    -e 's/stats_tout\s*=\s*.*/stats_tout = '$test_ts'/g' >> $NAMELISTS
    # Use this version if statistical output is not desired.
    cat $model_file | sed -e 's/dtmain\s*=\s*.*/dtmain = '$test_ts'/g' \
                        -e 's/dtclosure\s*=\s*.*/dtclosure = '$test_ts'/g' \
                        -e 's/l_stats\s*=\s*.*/l_stats = .false./g' >> $NAMELISTS
    cat $stats_file >> $NAMELISTS

    run_case
elif [ $ZT_GRID == true ];
then
    cat $parameter_file > $NAMELISTS
    cat $model_file | sed -e 's/^nzmax\s*=\s*.*//g' \
		-e 's/^grid_type\s*=\s*.*//g' \
		-e 's/^zm_grid_fname\s*=\s*.*//g' \
		-e "s/^zt_grid_fname\s*=\s*.*//g" \
		-e 's/^\&model_setting/\&model_setting\n \
		nzmax = '$grid_nz'\n \
		zt_grid_fname ='\'$grid_path\''\n \
		grid_type = 2\n/g' >> $NAMELISTS
    cat $stats_file >> $NAMELISTS

    run_case
elif [ $ZM_GRID == true ];
then
    cat $parameter_file > $NAMELISTS
    cat $model_file | sed -e 's/^nzmax\s*=\s*.*//g' \
        	-e 's/^grid_type\s*=\s*.*//g' \
        	-e 's/^zt_grid_fname\s*=\s*.*//g' \
        	-e 's/^zm_grid_fname\s*=\s*.*//g' \
        	-e 's/^\&model_setting/\&model_setting\n \
        	nzmax = '$grid_nz'\n \
        	zm_grid_fname ='\'$grid_path\''\n \
        	grid_type = 3\n/g' >> $NAMELISTS
    cat $stats_file >> $NAMELISTS

    run_case
elif [ $PERFORMANCE_TEST == true ];
then
    cat $parameter_file > $NAMELISTS
    cat $model_file | sed 's/l_stats\s*=\s*.*/l_stats = \.false\./g' \
                    | sed 's/debug_level\s*=\s*.*/debug_level = 0/g' \
                    >> $NAMELISTS
    cat $stats_file >> $NAMELISTS

    run_case
else
    #cat $parameter_file $model_file $stats_file > $NAMELISTS
    # This is a kluge for Fortran compilers that the can't handle comments in 
    # a namelist by using the sed command to remove them.
    cat $parameter_file $model_file $stats_file | sed 's/\!.*//' > $NAMELISTS

    run_case
fi

cd $restoreDir

exit $RESULT
