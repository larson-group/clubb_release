#!/bin/bash
#################################################################################
# $Id$
#
# Desciption:
#   Script to test how well CLUBB is converging. It runs RUN_CASE nine times
#   doubling the time step each time. Then, for a specified variable, the test
#   calculates the root mean square difference bewteen the dt = 1s case
#   and all other cases, to see how well CLUBB converges.
#################################################################################

# Set case to run, specify input file, and choose variable to test
RUN_CASE="bomex"
RUN_CASE_INPUT="../input/case_setups/${RUN_CASE}_model.in"
VAR_TO_TEST="rcm"

# Set the time steps (1, 2, 4, ..., 256)
TIMESTEP[0]=1
for i in `seq 1 8`
do
TIMESTEP[$i]=$((2*${TIMESTEP[$(($i-1))]}))
done

# Define output directories for the runs
for i in `seq 0 8`
do
OUTPUT_FOLDER[$i]="../output/timestep_${TIMESTEP[$i]}"
done

#save array length for later python processing
NUM_EXPS=${#TIMESTEP[@]}

# Figure out the directory where the script is located
scriptPath=`dirname $0`

# Store the current directory location so it can be restored
restoreDir=`pwd`

# Change directories to the one the script is located in
cd $scriptPath

# Copy the original input file so it can be restored
cp $RUN_CASE_INPUT ../input/case_setups/tmp_model.in

# Run for 90 minutes initially at 1s time step, then restart later
sed -i -e 's/time_final\s*=\s*.*/time_final = 5400.0/g' \
       -e 's/dt_main\s*=\s*.*/dt_main = 1/g' \
       -e 's/dt_rad\s*=\s*.*/dt_rad = 1/g' \
       -e 's/stats_tsamp\s*=\s*.*/stats_tsamp = 1/g' \
       -e 's/stats_tout\s*=\s*.*/stats_tout = 1/g' $RUN_CASE_INPUT

echo -n "Creating initial $RUN_CASE 90-min run... "
./run_scm.bash $RUN_CASE &> /dev/null
echo "Done!"

# Create restart folder and move 90-min run there
mkdir ../restart
mv ../output/*.* ../restart/

# Set case file for restart runs
sed -i -e 's/time_initial\s*=\s*.*/time_initial = 5400/g' \
       -e 's/time_final\s*=\s*.*/time_final = 9000.0/g' \
       -e 's/l_restart\s*=\s*.*/l_restart = .true./g' \
       -e 's/restart_path_case\s*=\s*.*/restart_path_case = "restart\/'$RUN_CASE'"/g' \
       -e 's/time_restart\s*=\s*.*/time_restart = 5400.0/g' $RUN_CASE_INPUT

# -------------------------------------------------------------------
# LOOP over the timesteps, saving output in respective output folders
for idx in `seq 0 8`
do

sed -i -e 's/dt_main\s*=\s*.*/dt_main = '${TIMESTEP[$idx]}'/g' \
       -e 's/dt_rad\s*=\s*.*/dt_rad = '${TIMESTEP[$idx]}'/g' \
       -e 's/stats_tsamp\s*=\s*.*/stats_tsamp = '${TIMESTEP[$idx]}'/g' \
       -e 's/stats_tout\s*=\s*.*/stats_tout = '${TIMESTEP[$idx]}'/g' $RUN_CASE_INPUT

echo -n "Restarting $RUN_CASE with time step ${TIMESTEP[$idx]} seconds... "
./run_scm.bash $RUN_CASE &> /dev/null
echo "Done!"

mkdir ${OUTPUT_FOLDER[$idx]}
mv ../output/*.* ${OUTPUT_FOLDER[$idx]}

CURR_OUTPUT_FOLDER=${OUTPUT_FOLDER[$idx]}

export RUN_CASE VAR_TO_TEST CURR_OUTPUT_FOLDER NUM_EXPS
  
# Extract and print var_to_test from last time step
python3 -c '
from netCDF4 import Dataset;
import numpy as np;
import os;
var_to_test=os.environ["VAR_TO_TEST"]
fout=open("convergence.data","a");
f = Dataset(os.environ["CURR_OUTPUT_FOLDER"]+"/"+os.environ["RUN_CASE"]+"_zt.nc");
test_var = np.array(f.variables[var_to_test]);
np.savetxt(fout,test_var[-1,:,0,0]);
'

done
# -------------------------------------------------------------------
# end LOOP


# create convergence plot
python3 -c '
from netCDF4 import Dataset;
import numpy as np;
import matplotlib.pyplot as plt;
import os;
import sys;
timesteps=[2,4,8,16,32,64,128,256];
for i in range(len(timesteps)):
    timesteps[i]=np.log(timesteps[i])
test_vars=np.loadtxt("convergence.data");
test_vars_len=len(test_vars);
num_exps=int(os.environ["NUM_EXPS"]);
test_vars=np.reshape(test_vars,[int(test_vars_len/num_exps),num_exps]);
rmse_values=[];
for i in range(1,num_exps):
    rmse_values.append(np.log(np.sqrt(np.var(test_vars[:,0]-test_vars[:,i]))));
print("Generating plot...");
m,b=np.polyfit(timesteps[0:4],rmse_values[0:4],1)
#fitted_points=[rmse_values[0]]
#slope_1_line=[rmse_values[0]]
#for i in range(1,4):
#    fitted_points.append(m*(timesteps[i]-timesteps[0])+rmse_values[0])
#    slope_1_line.append((timesteps[i]-timesteps[0])+rmse_values[0])
#plt.plot(timesteps[0:4],slope_1_line[0:4],"k--",label="slope-1 line")
#plt.plot(timesteps[0:4],fitted_points[0:4],"k:",label=os.environ["RUN_CASE"]+" fit (m = {:.2f}, 4 pts.)".format(m))
#plt.scatter(timesteps,rmse_values,label=os.environ["RUN_CASE"]+" data points");
#plt.xlabel("log(dt)");
#plt.ylabel("log(RMSE)");
#plt.title("Convergence test for "+os.environ["RUN_CASE"].upper()+" ("+os.environ["VAR_TO_TEST"]+")");
#plt.legend()
#plt.savefig("../output/"+os.environ["RUN_CASE"]+"_convergence.png",bbox_inches="tight");
print("slope = ",m)
if m > 0.5:
    sys.exit(0);
else:
    sys.exit(1)
'

if [ $? -ne 0 ]; then
  echo "CLUBB failed to converge!"
  RESULT=1
else
  echo "CLUBB converged!"
  RESULT=0
fi


# replace original model.in file
mv ../input/case_setups/tmp_model.in $RUN_CASE_INPUT

# cleanup the output
rm -rf ../restart
rm -rf ../output/timestep_*
rm convergence.data

cd $restoreDir

#how to determine convergence success?
exit $RESULT
