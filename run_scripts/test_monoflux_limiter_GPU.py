#!/bin/python

######################################################################################
# Description:
#   This script tests the differences between the CPU and GPU results of       
#   monotonic_turbulent_flux_limit, located in mono_flux_limiter.F90. The usual 
#   CPU vs GPU diffness test isn't sufficient to test this function because only 
#   the most unstable cases make use of monotonic_turbulent_flux_limit(), so it's
#   impossible to tell the differences between an error or accumulated bit changes
#   after enough timesteps to trigger the use of this function.
#
#   This test works by making use of NVIDIA's PCAST compare feature, where we use
#   the flag "-gpu=redundant" to compile the code so that it's run on both
#   CPU and GPU at the same time, then use "!$acc compare(var)" commands to test 
#   the differences between the CPU and GPU results. The PCAST feature keeps the
#   CPU and GPU results in sync from timestep to timstep, so the differences 
#   we measure here will not be mixed with accumulating bit changes. 
#
#   The PCAST feature relies on an environment variable "PCAST_COMPARE", that for
#   this test we set to "abs=10", which means that only differences above "10^(-10)"
#   will be reported. This was chosen by some small experiments. Usually correct
#   changes are more in the range of 10^-20 or smaller.
#
# Note:
#   This test does NOT confirm that changes made in monotonic_turbulent_flux_limit
#   are correct. This only confirms that there are no errors present that only
#   affect GPU results. Bit changing changes should be tested separately.
#
# References:
#   https://developer.nvidia.com/blog/detecting-divergence-using-pcast-to-compare-gpu-to-cpu-results/
#   https://github.com/larson-group/cam/issues/175#issuecomment-2259291053
#   
# Author:
#   Gunther Huebler
######################################################################################

import os
import re
import subprocess
import signal
import sys

config_file        = "../compile/config/linux_x86_64_nvhpc_gpu_openacc.bash"
source_file        = "../src/CLUBB_core/mono_flux_limiter.F90"
flags_file         = "../input/tunable_parameters/configurable_model_flags.in"

compile_script     = "../compile/compile.bash"
run_script         = "./run_scm.bash"

run_case           = "mc3e"     # unstable/noisy case is required
case_time_final    = "18000"    # corresponds to 60 timesteps for mc3e 
case_input         = "../input/case_setups/"+run_case+"_model.in"

err_code = 1    # default to error



##################################################################
#               Backup files to be changed
##################################################################
config_file_backup = config_file + ".bak"
source_backup = source_file + ".bak"
case_input_backup = case_input + ".bak"
flags_backup = flags_file + ".bak"
os.rename(config_file, config_file_backup)
os.rename(case_input, case_input_backup)
os.rename(source_file, source_backup)
os.rename(flags_file, flags_backup)

# This should be called to restore the backups
def restore_backups():

    os.rename(config_file_backup, config_file)
    os.rename(source_backup, source_file)
    os.rename(case_input_backup, case_input)
    os.rename(flags_backup, flags_file)

# This is called when kill signals are given.
def signal_handler(sig, frame):
    print("\n======== Script killed. Restoring and exiting. ========\n")
    restore_backups()
    sys.exit(0)

# Setup restore_backups to be called if Ctrl+C is hit 
signal.signal(signal.SIGINT, signal_handler)    # Handle Ctrl+C (SIGINT)
signal.signal(signal.SIGTERM, signal_handler)   # Handle termination (SIGTERM)
    


try:

    ##################################################################
    #         Modify the source file (mono_flux_limiter.F90)
    ##################################################################

    with open(source_backup, 'r') as infile, open(source_file, 'w') as outfile:
        for line in infile:
            # Remove the specific test comment
            line = line.replace("! MONOFLUX TEST COMMENT DO NOT REMOVE ", "")
            outfile.write(line)



    ##################################################################
    #                           Compile 
    ##################################################################

    # Modify the compiler config file
    with open(config_file_backup, 'r') as infile, open(config_file, 'w') as outfile:
        for line in infile:

            # Enable Nvidia PCAST features with (-gpu=redundant) 
            if line.startswith("FFLAGS"):
                line = line.replace("-acc", "-acc=gpu -gpu=redundant") 
            
            # SILHS uses a CUDA call that doesn't seem to work with -gpu=redundant
            line = line.replace("-DCUDA", "")

            # Use -O0 optimization just speed up compilation
            line = line.replace("-O2", "-O0")

            outfile.write(line)


    # Compile with the modified config file
    compile_command = [compile_script, "-c", config_file]
    process = subprocess.run(compile_command, stdout=subprocess.PIPE, text=True)

    # Display the output of the compilation
    print(process.stdout)



    ##################################################################
    #                       Edit Input Files 
    ##################################################################

    # Edit case _model.in to set time_final as specified above
    with open(case_input_backup, 'r') as infile, open(case_input, 'w') as outfile:
        for line in infile:
            line = re.sub(r"time_final\s*=\s*.*", "time_final = "+case_time_final, line)
            outfile.write(line)

    # Edit configurable_model_flags.in
    with open(flags_backup, 'r') as infile, open(flags_file, 'w') as outfile:
        for line in infile:
            line = re.sub(r"penta_solve_method\s*=.*", "penta_solve_method = 2", line)          # just for speed
            line = re.sub(r"tridiag_solve_method\s*=.*", "tridiag_solve_method = 2", line)      # just for speed
            line = re.sub(r"l_lh_straight_mc\s*=.*", "l_lh_straight_mc = .true.", line)         # required on GPUs for silhs cases
            outfile.write(line)



    ##################################################################
    #                           Run Case
    ##################################################################

    # Setting "PCAST_COMPARE=abs=n" causes the "acc compare" clauses used in this test 
    # to only report if the GPU results differ from CPU results by more than 10^-n
    # see https://developer.nvidia.com/blog/detecting-divergence-using-pcast-to-compare-gpu-to-cpu-results/
    os.environ["PCAST_COMPARE"] = "abs=10"

    # Run the case with -e, we only care about the stderr output
    run_command = ["bash", run_script, "-e", run_case]
    process = subprocess.run(run_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    # Print only errors
    print(process.stderr)


    ##################################################################
    #                   Check Output and Report Result
    ##################################################################

    # Check for relevant lines in the stderr output
    monoflux_wpxp_found = "MONOFLUX: wpxp adjusted" in process.stderr
    monoflux_xm_found   = "MONOFLUX: xm adjusted"   in process.stderr
    fail_abs_found      = "FAIL ABS"                in process.stderr


    print("\n==================================== RESULT ====================================")

    # If "FAIL ABS" appears, then we fail no matter what 
    # If there are no "adjusting" messages then monoflux_limiter was never really tested, we consider this a fail
    # If "adjusting" messages are found and "FAIL ABS" isnt' then we pass
    if fail_abs_found:
        print(f"\n\tTEST FAILED: 'FAIL ABS' found in the output.")
        print(f"\tGPU results differ too significantly from CPU results.\n")
    elif not (monoflux_wpxp_found and monoflux_xm_found):
        print(f"\n\tTEST FAILED: Neither 'MONOFLUX: wpxp adjusted' nor 'MONOFLUX: xm adjusted'")
        print(f"\twere found in the output. This means the flux limited wasn't tested.\n")
    else:
        err_code = 0
        print(f"\nTEST PASSED: mono_flux_limiter did modify fields, and CPU and GPU results match.\n")

except Exception as e:

    print("\n==================================== Test Incomplete ====================================")
    print(f"ERROR: {e}") 


restore_backups()
sys.exit(err_code)
