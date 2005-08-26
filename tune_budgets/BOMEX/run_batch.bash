####################################################################### 
# Tunes all the budget terms in sequence independently of one another.
# create_files.sh will need to be run before this if the directories 
# and namelists don't exist.
####################################################################### 

declare -ar CASES=([0]=wp2_bp [1]=wp2_prdp [2]=wp3_bp [3]=wp3_prdp [4]=wp3_tp\
 [5]=wprtp_bp [6]=wprtp_prdp [7]=wprtp_tp [8]=wpthlp_bp [9]=wpthlp_prdp\
 [10]=wpthlp_tp)

for (( i=0; i < "${#CASES[@]}"; i++ )); do
	cd "${CASES[$i]}"
        echo "Running ${CASES[$i]}"
	echo "n" | ../../hoc_tuner_budget_terms &> run.log
	cd ..
done
