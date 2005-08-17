# Generate namelists, etc. for a batch run.

declare -ar CASES=([0]=wp2_bp [1]=wp2_prdp [2]=wp3_bp [3]=wp3_prdp [4]=wp3_tp\
 [5]=wprtp_bp [6]=wprtp_prdp [7]=wprtp_tp [8]=wpthlp_bp [9]=wpthlp_prdp\
 [10]=wpthlp_tp)

for (( i=0; i < "${#CASES[@]}"; i++ )); do
	mkdir "${CASES[$i]}"
	cp budget.tmpl "${CASES[$i]}/budget.in"
	echo "&case_tune" >> "${CASES[$i]}/budget.in"
	for (( j=0; j < "${#CASES[@]}"; j++ )); do
		if [ "$i" -eq "$j" ]; then
			echo "l${CASES[$j]} = .true." >> "${CASES[$i]}/budget.in"
		else
			echo "l${CASES[$j]} = .false." >> "${CASES[$i]}/budget.in"
                fi
	done
	echo "/" >> "${CASES[$i]}/budget.in"
done


ln -s ../../tune .
ln -s ../../les_data .
