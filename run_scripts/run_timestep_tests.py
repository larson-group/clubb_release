#!/usr/bin/env python3
import subprocess
import sys
import os

def main():
    # Store the current directory
    restore_dir = os.getcwd()

    # Change directory to where this script is located
    script_path = os.path.dirname(os.path.abspath(__file__))
    os.chdir(script_path)

    test_cases = [
        "arm", "arm_97", "astex_a209", "atex", "bomex", "cgils_s6", "cgils_s11", "cgils_s12",
        "clex9_nov02", "clex9_oct14", "dycoms2_rf01", "dycoms2_rf01_fixed_sst",
        "dycoms2_rf02_do", "dycoms2_rf02_ds", "dycoms2_rf02_nd", "dycoms2_rf02_so",
        "fire", "gabls2", "gabls3", "gabls3_night", "jun25_altocu", "lba", "mc3e", "mpace_a", 
        "mpace_b", "mpace_b_silhs", "nov11_altocu", "rico", "rico_silhs", "twp_ice", "wangara"
    ]

    #test_cases = [ "arm", "gabls3" ]

    for case in test_cases:

        print(f"---------------- Running {case} ----------------")

        dt = 600

        #for i, dt in enumerate(test_timesteps):
        while ( dt < 3600 ):

            #print(f" -- Testing dt_main = {dt}")

            result = subprocess.run(
                ["./run_scm.py", "-debug", "0", "-tout", "0", "-dt_main", str(dt), "-dt_rad", str(dt), case],
                capture_output=True  # show output live, like bash script
            )

            if result.returncode != 0:
                print(f"---- FAIL @ dt = {dt}")
                break

            dt = dt + 600

        if ( result.returncode == 0 ):
            print(f"--- PASS up to dt = {dt}")


    # Restore original directory
    os.chdir(restore_dir)

if __name__ == "__main__":
    main()
