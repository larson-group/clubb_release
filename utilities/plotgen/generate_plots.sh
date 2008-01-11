#!/bin/bash

#Tell the user what we're doing
echo "This script will generate matlab plots comparing two sets of input data."
echo ""

#Lets make sure the arguments are worth something and store the data
if [ "$1" == "" ]; then
	echo "Please specify at least one set of data."
	#I'm using exit 1 for to exit under any error status
	exit 1
elif [ "$1" == "nightly" ]; then
	#If we're doing the nighthy plots, just use the default settings
	echo "Generating nightly plots."	
	
	#Get to the working directory
	cd /home/matlabuser/plotgen	

	echo "quit" | (matlab -nojvm -nodisplay -r compare_plots_cases_driver"( 'HOC_previous', 'HOC_current', 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 )")
	rm -rf /home/matlabuser/plotgen/profiles && \
	latex2html /home/matlabuser/plotgen/profiles.tex
	exit
elif [ "$2" == 0 ]; then
	#The generate_plots script expects two sets of sim data,
	#even if we are only plotting one set, if we pass it a zero it will
	#ignore the second set requirement
	echo "Only plotting one HOC simulation."
	echo "Set 1: $1"
	HOC_sim1="$1"
	HOC_sim2=0
	echo "Plotting HOC simulation results in the directory $HOC_sim1"
else
	echo "Set 1: $1"
	HOC_sim1="$1"

	#Make sure a second set of sim data is actually specified
	if [ "$2" == "" ]; then
		echo "Please specify a second set of data."
		exit 1
	else
		echo "Set 2: $2"
		HOC_sim2="$2"
		echo "Plotting HOC simulation results in the directories $HOC_sim1 and $HOC_sim2"
	fi
fi

#Clean up the output directories, if we don't do this, we can't
#pick and choose what cases will be shown in the final product
rm -rf /home/matlabuser/plotgen/output/eps/*
rm -rf /home/matlabuser/plotgen/output/jpg/*
rm -rf /home/matlabuser/plotgen/output/ps/*

#Set the flag that determines if the plots are compared to the LES data
if [ "$3" == "1" ]; then
	echo "Plotting LES data."
	compare_LES=1
else
	echo "Not plotting LES data."
	compare_LES=0
fi

#Set the flag that determines if the plots are compared to the Chriz Golaz best ever data
if [ "$4" == "1" ]; then
	echo "Plotting Chris Golaz best ever."
	compare_best=1
else
	echo "Not plotting Chris Golaz best ever."
	compare_best=0
fi

#Set the flag that determines if the plots are compared to the 25 Dec 2005 data
if [ "$5" == "1" ]; then
	echo "Plotting Dec. 17, 2005 HOC data."
	compare_HOC=1
else
	echo "Not plotting Dec. 17, 2005 HOC data."
	compare_HOC=0
fi

#The following section makes sure we only plot the 
#cases we have data for

#Theses variables are essentially boolean, 0 = false, 1 = true
#They are zeroed so there is a known value they start at
plot_arm_sim1=0
plot_atex_sim1=0
plot_bomex_sim1=0
plot_cobra_sim1=0
plot_dycoms2_rf01_sim1=0
plot_dycoms2_rf02_do_sim1=0
plot_dycoms2_rf02_ds_sim1=0
plot_dycoms2_rf02_nd_sim1=0
plot_dycoms2_rf02_so_sim1=0
plot_fire_sim1=0
plot_gabls2_sim1=0
plot_jun25_altocu_sim1=0
plot_mpace_b_sim1=0
plot_nov11_altocu_sim1=0
plot_rico_sim1=0
plot_wangara_sim1=0

plot_arm_sim2=0
plot_atex_sim2=0
plot_bomex_sim2=0
plot_cobra_sim2=0
plot_dycoms2_rf01_sim2=0
plot_dycoms2_rf02_do_sim2=0
plot_dycoms2_rf02_ds_sim2=0
plot_dycoms2_rf02_nd_sim2=0
plot_dycoms2_rf02_so_sim2=0
plot_fire_sim2=0
plot_gabls2_sim2=0
plot_jun25_altocu_sim2=0
plot_mpace_b_sim2=0
plot_nov11_altocu_sim2=0
plot_rico_sim2=0
plot_wangara_sim2=0


#If the Arm files exist, plot Arm
if [ -r "$HOC_sim1/arm_zm.ctl" ]; then
	if [ -r "$HOC_sim1/arm_zm.dat" ]; then
		if [ -r "$HOC_sim1/arm_zt.ctl" ]; then
			if [ -r "$HOC_sim1/arm_zt.dat" ]; then
				echo "Plotting Arm for HOC_sim1"
				plot_arm_sim1=1
			fi
		fi
	fi
fi

#If the Atex files exist, plot Atex
if [ -r "$HOC_sim1/atex_zm.ctl" ]; then
	if [ -r "$HOC_sim1/atex_zm.dat" ]; then
		if [ -r "$HOC_sim1/atex_zt.ctl" ]; then
			if [ -r "$HOC_sim1/atex_zt.dat" ]; then
				echo "Plotting Atex for HOC_sim1"
				plot_atex_sim1=1
			fi
		fi
	fi
fi

#If the Bomex files exist, plot Bomex
if [ -r "$HOC_sim1/bomex_zm.ctl" ]; then
	if [ -r "$HOC_sim1/bomex_zm.dat" ]; then
		if [ -r "$HOC_sim1/bomex_zt.ctl" ]; then
			if [ -r "$HOC_sim1/bomex_zt.dat" ]; then
				echo "Plotting Bomex for HOC_sim1"
				plot_bomex_sim1=1
			fi
		fi
	fi
fi

#If the Cobra files exist, plot Cobra
if [ -r "$HOC_sim1/cobra_zm.ctl" ]; then
	if [ -r "$HOC_sim1/cobra_zm.dat" ]; then
		if [ -r "$HOC_sim1/cobra_zt.ctl" ]; then
			if [ -r "$HOC_sim1/cobra_zt.dat" ]; then
				echo "Plotting Cobra for HOC_sim1"
				plot_cobra_sim1=1
			fi
		fi
	fi
fi

#If the Dycoms2_RF01 files exist, plot Dycoms2_RF01
if [ -r "$HOC_sim1/dycoms2_rf01_zm.ctl" ]; then
	if [ -r "$HOC_sim1/dycoms2_rf01_zm.dat" ]; then
		if [ -r "$HOC_sim1/dycoms2_rf01_zt.ctl" ]; then
			if [ -r "$HOC_sim1/dycoms2_rf01_zt.dat" ]; then
				echo "Plotting Dycoms2_RF01 for HOC_sim1"
				plot_dycoms2_rf01_sim1=1
			fi
		fi
	fi
fi

#If the Dycoms2_RF02_DO files exist, plot Dycoms2_RF02_DO
if [ -r "$HOC_sim1/dycoms2_rf02_do_zm.ctl" ]; then
	if [ -r "$HOC_sim1/dycoms2_rf02_do_zm.dat" ]; then
		if [ -r "$HOC_sim1/dycoms2_rf02_do_zt.ctl" ]; then
			if [ -r "$HOC_sim1/dycoms2_rf02_do_zt.dat" ]; then
				echo "Plotting Dycoms2_RF02_DO for HOC_sim1"
				plot_dycoms2_rf02_do_sim1=1
			fi
		fi
	fi
fi

#If the Dycoms2_RF02_DS files exist, plot Dycoms2_RF02_DS
if [ -r "$HOC_sim1/dycoms2_rf02_ds_zm.ctl" ]; then
	if [ -r "$HOC_sim1/dycoms2_rf02_ds_zm.dat" ]; then
		if [ -r "$HOC_sim1/dycoms2_rf02_ds_zt.ctl" ]; then
			if [ -r "$HOC_sim1/dycoms2_rf02_ds_zt.dat" ]; then
				echo "Plotting Dycoms2_RF02_DS for HOC_sim1"
				plot_dycoms2_rf02_ds_sim1=1
			fi
		fi
	fi
fi

#If the Dycoms2_RF02_ND files exist, plot Dycoms2_RF02_ND
if [ -r "$HOC_sim1/dycoms2_rf02_nd_zm.ctl" ]; then
	if [ -r "$HOC_sim1/dycoms2_rf02_nd_zm.dat" ]; then
		if [ -r "$HOC_sim1/dycoms2_rf02_nd_zt.ctl" ]; then
			if [ -r "$HOC_sim1/dycoms2_rf02_nd_zt.dat" ]; then
				echo "Plotting Dycoms2_RF02_ND for HOC_sim1"
				plot_dycoms2_rf02_nd_sim1=1
			fi
		fi
	fi
fi

#If the Dycoms2_RF02_SO files exist, plot Dycoms2_RF02_SO
if [ -r "$HOC_sim1/dycoms2_rf02_so_zm.ctl" ]; then
	if [ -r "$HOC_sim1/dycoms2_rf02_so_zm.dat" ]; then
		if [ -r "$HOC_sim1/dycoms2_rf02_so_zt.ctl" ]; then
			if [ -r "$HOC_sim1/dycoms2_rf02_so_zt.dat" ]; then
				echo "Plotting Dycoms2_RF02_SO for HOC_sim1"
				plot_dycoms2_rf02_so_sim1=1
			fi
		fi
	fi
fi

#If the Fire files exist, plot Fire
if [ -r "$HOC_sim1/fire_zm.ctl" ]; then
	if [ -r "$HOC_sim1/fire_zm.dat" ]; then
		if [ -r "$HOC_sim1/fire_zt.ctl" ]; then
			if [ -r "$HOC_sim1/fire_zt.dat" ]; then
				echo "Plotting Fire for HOC_sim1"
				plot_fire_sim1=1
			fi
		fi
	fi
fi

#If the Gabls2 files exist, plot Gabls2
if [ -r "$HOC_sim1/gabls2_zm.ctl" ]; then
	if [ -r "$HOC_sim1/gabls2_zm.dat" ]; then
		if [ -r "$HOC_sim1/gabls2_zt.ctl" ]; then
			if [ -r "$HOC_sim1/gabls2_zt.dat" ]; then
				echo "Plotting Gabls2 for HOC_sim1"
				plot_gabls2_sim1=1
			fi
		fi
	fi
fi

#If the Jun25_Altocu files exist, plot Jun25_Altocu
if [ -r "$HOC_sim1/jun25_altocu_zm.ctl" ]; then
	if [ -r "$HOC_sim1/jun25_altocu_zm.dat" ]; then
		if [ -r "$HOC_sim1/jun25_altocu_zt.ctl" ]; then
			if [ -r "$HOC_sim1/jun25_altocu_zt.dat" ]; then
				echo "Plotting Jun25_Altocu for HOC_sim1"
				plot_jun25_altocu_sim1=1
			fi
		fi
	fi
fi

#If the Mpace_B files exist, plot Mpace_B
if [ -r "$HOC_sim1/mpace_b_zm.ctl" ]; then
	if [ -r "$HOC_sim1/mpace_b_zm.dat" ]; then
		if [ -r "$HOC_sim1/mpace_b_zt.ctl" ]; then
			if [ -r "$HOC_sim1/mpace_b_zt.dat" ]; then
				echo "Plotting Mpace_B for HOC_sim1"
				plot_mpace_b_sim1=1
			fi
		fi
	fi
fi

#If the Nov11_Altocu files exist, plot Nov11_Altocu
if [ -r "$HOC_sim1/nov11_altocu_zm.ctl" ]; then
	if [ -r "$HOC_sim1/nov11_altocu_zm.dat" ]; then
		if [ -r "$HOC_sim1/nov11_altocu_zt.ctl" ]; then
			if [ -r "$HOC_sim1/nov11_altocu_zt.dat" ]; then
				echo "Plotting Nov11_Altocu for HOC_sim1"
				plot_nov11_altocu_sim1=1
			fi
		fi
	fi
fi

#If the Rico files exist, plot Rico
if [ -r "$HOC_sim1/rico_zm.ctl" ]; then
	if [ -r "$HOC_sim1/rico_zm.dat" ]; then
		if [ -r "$HOC_sim1/rico_zt.ctl" ]; then
			if [ -r "$HOC_sim1/rico_zt.dat" ]; then
				echo "Plotting Rico for HOC_sim1"
				plot_rico_sim1=1
			fi
		fi
	fi
fi

#If the Wangara files exist, plot Wangara
if [ -r "$HOC_sim1/wangara_zm.ctl" ]; then
	if [ -r "$HOC_sim1/wangara_zm.dat" ]; then
		if [ -r "$HOC_sim1/wangara_zt.ctl" ]; then
			if [ -r "$HOC_sim1/wangara_zt.dat" ]; then
				echo "Plotting Wangara for HOC_sim1"
				plot_wangara_sim1=1
			fi
		fi
	fi
fi

#Check if the necessary files exist for HOC_sim2

if [ "$HOC_sim2" != 0 ]; then
	
	#If the Arm files exist, plot Arm
	if [ -r "$HOC_sim2/arm_zm.ctl" ]; then
		if [ -r "$HOC_sim2/arm_zm.dat" ]; then
			if [ -r "$HOC_sim2/arm_zt.ctl" ]; then
				if [ -r "$HOC_sim2/arm_zt.dat" ]; then
					echo "Plotting Arm for HOC_sim2"
					plot_arm_sim2=1
				fi
			fi
		fi
	fi

	#If the Atex files exist, plot Atex
	if [ -r "$HOC_sim2/atex_zm.ctl" ]; then
		if [ -r "$HOC_sim2/atex_zm.dat" ]; then
			if [ -r "$HOC_sim2/atex_zt.ctl" ]; then
				if [ -r "$HOC_sim2/atex_zt.dat" ]; then
					echo "Plotting Atex for HOC_sim2"
					plot_atex_sim2=1
				fi
			fi
		fi
	fi

	#If the Bomex files exist, plot Bomex
	if [ -r "$HOC_sim2/bomex_zm.ctl" ]; then
		if [ -r "$HOC_sim2/bomex_zm.dat" ]; then
			if [ -r "$HOC_sim2/bomex_zt.ctl" ]; then
				if [ -r "$HOC_sim2/bomex_zt.dat" ]; then
					echo "Plotting Bomex for HOC_sim2"
					plot_bomex_sim2=1
				fi
			fi
		fi
	fi

	#If the Cobra files exist, plot Cobra
	if [ -r "$HOC_sim2/cobra_zm.ctl" ]; then
		if [ -r "$HOC_sim2/cobra_zm.dat" ]; then
			if [ -r "$HOC_sim2/cobra_zt.ctl" ]; then
				if [ -r "$HOC_sim2/cobra_zt.dat" ]; then
					echo "Plotting Cobra for HOC_sim2"
					plot_cobra_sim2=1
				fi
			fi
		fi
	fi

	#If the Dycoms2_RF01 files exist, plot Dycoms2_RF01
	if [ -r "$HOC_sim2/dycoms2_rf01_zm.ctl" ]; then
		if [ -r "$HOC_sim2/dycoms2_rf01_zm.dat" ]; then
			if [ -r "$HOC_sim2/dycoms2_rf01_zt.ctl" ]; then
				if [ -r "$HOC_sim2/dycoms2_rf01_zt.dat" ]; then
					echo "Plotting Dycoms2_RF01 for HOC_sim2"
					plot_dycoms2_rf01_sim2=1
				fi
			fi
		fi
	fi

	#If the Dycoms2_RF02_DO files exist, plot Dycoms2_RF02_DO
	if [ -r "$HOC_sim2/dycoms2_rf02_do_zm.ctl" ]; then
		if [ -r "$HOC_sim2/dycoms2_rf02_do_zm.dat" ]; then
			if [ -r "$HOC_sim2/dycoms2_rf02_do_zt.ctl" ]; then
				if [ -r "$HOC_sim2/dycoms2_rf02_do_zt.dat" ]; then
					echo "Plotting Dycoms2_RF02_DO for HOC_sim2"
					plot_dycoms2_rf02_do_sim2=1
				fi
			fi
		fi
	fi

	#If the Dycoms2_RF02_DS files exist, plot Dycoms2_RF02_DS
	if [ -r "$HOC_sim2/dycoms2_rf02_ds_zm.ctl" ]; then
		if [ -r "$HOC_sim2/dycoms2_rf02_ds_zm.dat" ]; then
			if [ -r "$HOC_sim2/dycoms2_rf02_ds_zt.ctl" ]; then
				if [ -r "$HOC_sim2/dycoms2_rf02_ds_zt.dat" ]; then
					echo "Plotting Dycoms2_RF02_DS for HOC_sim2"
					plot_dycoms2_rf02_ds_sim2=1
				fi
			fi
		fi
	fi

	#If the Dycoms2_RF02_ND files exist, plot Dycoms2_RF02_ND
	if [ -r "$HOC_sim2/dycoms2_rf02_nd_zm.ctl" ]; then
		if [ -r "$HOC_sim2/dycoms2_rf02_nd_zm.dat" ]; then
			if [ -r "$HOC_sim2/dycoms2_rf02_nd_zt.ctl" ]; then
				if [ -r "$HOC_sim2/dycoms2_rf02_nd_zt.dat" ]; then
					echo "Plotting Dycoms2_RF02_ND for HOC_sim2"
					plot_dycoms2_rf02_nd_sim2=1
				fi
			fi
		fi
	fi

	#If the Dycoms2_RF02_SO files exist, plot Dycoms2_RF02_SO
	if [ -r "$HOC_sim2/dycoms2_rf02_so_zm.ctl" ]; then
		if [ -r "$HOC_sim2/dycoms2_rf02_so_zm.dat" ]; then
			if [ -r "$HOC_sim2/dycoms2_rf02_so_zt.ctl" ]; then
				if [ -r "$HOC_sim2/dycoms2_rf02_so_zt.dat" ]; then
					echo "Plotting Dycoms2_RF02_SO for HOC_sim2"
					plot_dycoms2_rf02_so_sim2=1
				fi
			fi
		fi
	fi

	#If the Fire files exist, plot Fire
	if [ -r "$HOC_sim2/fire_zm.ctl" ]; then
		if [ -r "$HOC_sim2/fire_zm.dat" ]; then
			if [ -r "$HOC_sim2/fire_zt.ctl" ]; then
				if [ -r "$HOC_sim2/fire_zt.dat" ]; then
					echo "Plotting Fire for HOC_sim2"
					plot_fire_sim2=1
				fi
			fi
		fi
	fi

	#If the Gabls2 files exist, plot Gabls2
	if [ -r "$HOC_sim2/gabls2_zm.ctl" ]; then
		if [ -r "$HOC_sim2/gabls2_zm.dat" ]; then
			if [ -r "$HOC_sim2/gabls2_zt.ctl" ]; then
				if [ -r "$HOC_sim2/gabls2_zt.dat" ]; then
					echo "Plotting Gabls2 for HOC_sim2"
					plot_gabls2_sim2=1
				fi
			fi
		fi
	fi

	#If the Jun25_Altocu files exist, plot Jun25_Altocu
	if [ -r "$HOC_sim2/jun25_altocu_zm.ctl" ]; then
		if [ -r "$HOC_sim2/jun25_altocu_zm.dat" ]; then
			if [ -r "$HOC_sim2/jun25_altocu_zt.ctl" ]; then
				if [ -r "$HOC_sim2/jun25_altocu_zt.dat" ]; then
					echo "Plotting Jun25_Altocu for HOC_sim2"
					plot_jun25_altocu_sim2=1
				fi
			fi
		fi
	fi

	#If the Mpace_B files exist, plot Mpace_B
	if [ -r "$HOC_sim2/mpace_b_zm.ctl" ]; then
		if [ -r "$HOC_sim2/mpace_b_zm.dat" ]; then
			if [ -r "$HOC_sim2/mpace_b_zt.ctl" ]; then
				if [ -r "$HOC_sim2/mpace_b_zt.dat" ]; then
					echo "Plotting Mpace_B for HOC_sim2"
					plot_mpace_b_sim2=1
				fi
			fi
		fi
	fi

	#If the Nov11_Altocu files exist, plot Nov11_Altocu
	if [ -r "$HOC_sim2/nov11_altocu_zm.ctl" ]; then
		if [ -r "$HOC_sim2/nov11_altocu_zm.dat" ]; then
			if [ -r "$HOC_sim2/nov11_altocu_zt.ctl" ]; then
				if [ -r "$HOC_sim2/nov11_altocu_zt.dat" ]; then
					echo "Plotting Nov11_Altocu for HOC_sim2"
					plot_nov11_altocu_sim2=1
				fi
			fi
		fi
	fi

	#If the Rico files exist, plot Rico
	if [ -r "$HOC_sim2/rico_zm.ctl" ]; then
		if [ -r "$HOC_sim2/rico_zm.dat" ]; then
			if [ -r "$HOC_sim2/rico_zt.ctl" ]; then
				if [ -r "$HOC_sim2/rico_zt.dat" ]; then
					echo "Plotting Rico for HOC_sim2"
					plot_rico_sim2=1
				fi
			fi
		fi
	fi

	#If the Wangara files exist, plot Wangara
	if [ -r "$HOC_sim2/wangara_zm.ctl" ]; then
		if [ -r "$HOC_sim2/wangara_zm.dat" ]; then
			if [ -r "$HOC_sim2/wangara_zt.ctl" ]; then
				if [ -r "$HOC_sim2/wangara_zt.dat" ]; then
					echo "Plotting Wangara for HOC_sim2"
					plot_wangara_sim2=1
				fi
			fi
		fi
	fi

fi
# End HOC_sim2 file checks

#This variable is used to ensure we only generate the HTML if the necessary
#data files have actually been produced
run_success=0

cd /home/matlabuser/plotgen

#Actually run the script with the arguments we've parsed out
if [ "$HOC_sim2" == 0 ]; then
	echo "quit" | (matlab -nojvm -nodisplay -r compare_plots_cases_driver"( '$HOC_sim1', 0, $compare_LES, $compare_best, $compare_HOC, $plot_arm_sim1, $plot_atex_sim1, $plot_bomex_sim1, $plot_cobra_sim1, $plot_dycoms2_rf01_sim1, $plot_dycoms2_rf02_do_sim1, $plot_dycoms2_rf02_ds_sim1, $plot_dycoms2_rf02_nd_sim1, $plot_dycoms2_rf02_so_sim1, $plot_fire_sim1, $plot_gabls2_sim1, $plot_jun25_altocu_sim1, $plot_mpace_b_sim1, $plot_nov11_altocu_sim1, $plot_rico_sim1, $plot_wangara_sim1, $plot_arm_sim2, $plot_atex_sim2, $plot_bomex_sim2, $plot_cobra_sim2, $plot_dycoms2_rf01_sim2, $plot_dycoms2_rf02_do_sim2, $plot_dycoms2_rf02_ds_sim2, $plot_dycoms2_rf02_nd_sim2, $plot_dycoms2_rf02_so_sim2, $plot_fire_sim2, $plot_gabls2_sim2, $plot_jun25_altocu_sim2, $plot_mpace_b_sim2, $plot_nov11_altocu_sim2, $plot_rico_sim2, $plot_wangara_sim2 )") && \
	run_success=1
else
	echo "quit" | (matlab -nojvm -nodisplay -r compare_plots_cases_driver"( '$HOC_sim1', '$HOC_sim2', $compare_LES, $compare_best, $compare_HOC, $plot_arm_sim1, $plot_atex_sim1, $plot_bomex_sim1, $plot_cobra_sim1, $plot_dycoms2_rf01_sim1, $plot_dycoms2_rf02_do_sim1, $plot_dycoms2_rf02_ds_sim1, $plot_dycoms2_rf02_nd_sim1, $plot_dycoms2_rf02_so_sim1, $plot_fire_sim1, $plot_gabls2_sim1, $plot_jun25_altocu_sim1, $plot_mpace_b_sim1, $plot_nov11_altocu_sim1, $plot_rico_sim1, $plot_wangara_sim1, $plot_arm_sim2, $plot_atex_sim2, $plot_bomex_sim2, $plot_cobra_sim2, $plot_dycoms2_rf01_sim2, $plot_dycoms2_rf02_do_sim2, $plot_dycoms2_rf02_ds_sim2, $plot_dycoms2_rf02_nd_sim2, $plot_dycoms2_rf02_so_sim2, $plot_fire_sim2, $plot_gabls2_sim2, $plot_jun25_altocu_sim2, $plot_mpace_b_sim2, $plot_nov11_altocu_sim2, $plot_rico_sim2, $plot_wangara_sim2 )") && \
	run_success=1
fi

#Generate the plots
if [ "$run_success" == 1 ]; then
	rm -rf /home/matlabuser/plotgen/profiles
	latex2html /home/matlabuser/plotgen/profiles.tex
	exit
else
	exit 1
fi
