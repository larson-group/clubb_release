####################################################################### 
# Does a tuning run for each term on remote nodes.
# create_files.sh will need to be run before this.
####################################################################### 
BASE_DIR=$HOME"/hoc_v2.2_tuner"
BUDGET_DIR=$BASE_DIR"/tune_budget/BOMEX_parallel"

# All cases
declare -ar CASES=([0]=wp2_bp [1]=wp2_prdp [2]=wp3_bp [3]=wp3_prdp [4]=wp3_tp\
 [5]=wprtp_bp [6]=wprtp_prdp [7]=wprtp_tp [8]=wpthlp_bp [9]=wpthlp_prdp\
 [10]=wpthlp_tp [11]=rtp2_dp [12]=rtp2_tp [13]=thlp2_tp [14]=thlp2_dp\
 [15]=rtpthlp_tp [16]=rtpthlp_dp )

# This array needs to match up with the create_files script to work
declare -ar NODES=([0]=tom01 [1]=tom02 [2]=tom03 [3]=tom04 [4]=tom05\
 [5]=tom06 [6]=tom07 [7]=tom08 [8]=tom09 [9]=tom10\
 [10]=tom11 [11]=tom12 [12]=tom13 [13]=tom14 [14]=tom15\
 [15]=tom16 [16]=tom17 )

for (( i=0; i < "${#CASES[@]}"; i++ )); do
#for (( i=0; i < 1; i++ )); do
        echo "Running ${CASES[$i]} on ${NODES[$i]}"
	rsh ${NODES[$i]} "cd $BUDGET_DIR/${CASES[$i]} && echo 'n' | ../hoc_tuner_budget_terms &> run.log " &
done
