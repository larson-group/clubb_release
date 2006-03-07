#!/bin/bash
# Generate namelists, etc. for a batch run.

declare -ar CASES=([0]=wp2_bp [1]=wp2_prdp [2]=wp3_bp [3]=wp3_prdp [4]=wp3_tp\
 [5]=wprtp_bp [6]=wprtp_prdp [7]=wprtp_tp [8]=wpthlp_bp [9]=wpthlp_prdp\
 [10]=wpthlp_tp [11]=rtp2_dp [12]=rtp2_tp [13]=thlp2_tp [14]=thlp2_dp\
 [15]=rtpthlp_tp [16]=rtpthlp_dp )

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

cat "../../model/bomex_model.in" "../../stats/bomex_stats.in" > "bomex_hoc.in"
ln -s ../../les_data .
