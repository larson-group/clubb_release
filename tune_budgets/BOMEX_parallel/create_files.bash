# Generate namelists, etc. for a batch run.
# This version modified to use a cluster of nodes

# Edit these for your particular configuration
BASE_DIR=$HOME"/hoc_v2.2_tuner"
LES_DIR=$BASE_DIR"/les_data"
NML_DIR=$BASE_DIR"/tune"
BUDGET_DIR=$BASE_DIR"/tune_budget/BOMEX_parallel"
LES_SM=$LES_DIR"/bomex_coamps_sm"
LES_SW=$LES_DIR"/bomex_coamps_sw"
MODEL_NML=$NML_DIR"/bomex_hoc.in"

# All the cases to be run
declare -ar CASES=([0]=wp2_bp [1]=wp2_prdp [2]=wp3_bp [3]=wp3_prdp [4]=wp3_tp\
 [5]=wprtp_bp [6]=wprtp_prdp [7]=wprtp_tp [8]=wpthlp_bp [9]=wpthlp_prdp\
 [10]=wpthlp_tp [11]=rtp2_dp [12]=rtp2_tp [13]=thlp2_tp [14]=thlp2_dp\
 [15]=rtpthlp_tp [16]=rtpthlp_dp )

# And the nodes for them to be run on (i.e. element 0 of NODES will be
# doing element 0 of CASES)  To make a node do more than one CASE,
# you'll want to specify it multiple times in NODES (e.g. [0]=tom01 [1]=tom01)

declare -ar NODES=([0]=tom01 [1]=tom02 [2]=tom03 [3]=tom04 [4]=tom05\
 [5]=tom06 [6]=tom07 [7]=tom08 [8]=tom09 [9]=tom10\
 [10]=tom11)

# The actual script

# The -p option should create previous directories in the tree as needed on
# UNIX and GNU systems, but this will require modification if it does not.
for (( i=0; i < "${#CASES[@]}"; i++ )); do
#for (( i=0; i < 1; i++ )); do
	rsh ${NODES[$i]} "mkdir -p $LES_DIR"
	rsh ${NODES[$i]} "mkdir -p $NML_DIR"
	rsh ${NODES[$i]} "mkdir -p $BUDGET_DIR"

	rcp ../hoc_tuner_budget_terms "${NODES[$i]}:$BUDGET_DIR"
	rcp $LES_SM.??? "${NODES[$i]}:$LES_DIR/"
	rcp $LES_SW.??? "${NODES[$i]}:$LES_DIR/"
	rcp $NML_DIR/rand_seed.dat "${NODES[$i]}:$NML_DIR"
	rcp $MODEL_NML "${NODES[$i]}:$MODEL_NML"

	rsh "${NODES[$i]}" "mkdir $BUDGET_DIR/${CASES[$i]}"
	rcp budget.tmpl "${NODES[$i]}:$BUDGET_DIR/${CASES[$i]}/budget.in"
#	rcp run_node.bash "${NODES[$i]}:$BUDGET_DIR/"

#  Append the last namelist to the budget.in file for all CASES
	rsh "${NODES[$i]}" "echo '&case_tune' >> $BUDGET_DIR/${CASES[$i]}/budget.in"
	for (( j=0; j < "${#CASES[@]}"; j++ )); do
		if [ "$i" -eq "$j" ]; then
			rsh ${NODES[$i]} "echo 'l${CASES[$j]} = .true.' >> $BUDGET_DIR/${CASES[$i]}/budget.in"
		else
			rsh ${NODES[$i]} "echo 'l${CASES[$j]} = .false.' >> $BUDGET_DIR/${CASES[$i]}/budget.in"
		fi
	done
	rsh ${NODES[$i]} "echo '/' >> $BUDGET_DIR/${CASES[$i]}/budget.in"
# End namelist

# Make symlinks to the location of the hoc.in and LES data dir
	rsh ${NODES[$i]} "cd $BUDGET_DIR && ln -s ../../tune ."
	rsh ${NODES[$i]} "cd $BUDGET_DIR && ln -s ../../les_data ."
done
