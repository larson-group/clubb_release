#!/bin/python3

import argparse
import math
import subprocess
import re
import numpy as np
import os
import sys
import time
import socket

def update_time_final(file_path, iterations):

    # Patterns to extract values
    time_initial_pattern = re.compile(r'time_initial\s*=\s*([0-9.]+)')
    dt_main_pattern = re.compile(r'dt_main\s*=\s*([0-9.]+)')
    time_final_pattern = re.compile(r'(time_final\s*=\s*)([0-9.]+)(.*)')

    # Read file content
    with open(file_path, 'r') as file:
        lines = file.readlines()

    time_initial = None
    dt_main = None

    # Parse values
    for line in lines:
        if time_initial is None:
            match = time_initial_pattern.search(line)
            if match:
                time_initial = float(match.group(1))

        if dt_main is None:
            match = dt_main_pattern.search(line)
            if match:
                dt_main = float(match.group(1))

        if time_initial is not None and dt_main is not None:
            break

    if time_initial is None or dt_main is None:
        raise ValueError("Could not find 'time_initial' or 'dt_main' in the file.")

    # Calculate new time_final
    new_time_final = time_initial + dt_main * iterations

    # Update the file content
    updated_lines = []
    for line in lines:
        match = time_final_pattern.search(line)
        if match:
            new_line = f"{match.group(1)}{new_time_final:.1f}    ! Updated time_final\n"
            updated_lines.append(new_line)
        else:
            updated_lines.append(line)

    # Write updated content back to the file
    with open(file_path, 'w') as file:
        file.writelines(updated_lines)

    #print(f"Updated 'time_final' to {new_time_final:.1f} in {file_path}")

def get_git_hash():
    # Get the current git hash of the repository the script is in
    try:
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            cwd=os.path.dirname(os.path.abspath(__file__))
        )
        if result.returncode == 0:
            return result.stdout.strip()
        else:
            return f"Git hash retrieval failed: {result.stderr.strip()}"
    except Exception as e:
        return f"Error retrieving git hash: {e}"

def clean_float_string(s):
    # Use regex to keep valid characters for floating-point numbers
    cleaned = re.sub(r"[^0-9\.]", "", s)
    return cleaned.strip()

def time_clubb_standalone(case: str, ngrdcol: int, mpi: int, srun: int, gpu: bool):

    timing_results = {}
    timing_results['ngrdcol'] = ngrdcol
    timing_results['total'] = 0.0
    timing_results['compute_i'] = 0.0
    timing_results['baseline'] = 0.0
    timing_results['clubb_advance'] = 0.0
    timing_results['output_multi_col'] = 0.0
    timing_results['mainloop'] = 0.0
    timing_results['iterations'] = 0.0

    scriptDir = os.path.dirname(os.path.abspath(__file__))

    if mpi > 0:
        if gpu:
            if any('openmpi' in os.environ.get(var, '').lower() for var in os.environ):
                print(" - OpenMPI is being used.")
                run_cmd = f"mpirun -np {min(ngrdcol,mpi)} --bind-to hwthread --map-by numa bash -c \"cd {scriptDir}; set_gpu_rank ./execute_clubb_standalone.bash\""
            elif any('mpich' in os.environ.get(var, '').lower() for var in os.environ):
                print(" - MPICH is being used.")
                run_cmd = f"mpirun -np {min(ngrdcol,mpi)} --cpu-bind core --mem-bind local bash -c \"cd {scriptDir}; set_gpu_rank ./execute_clubb_standalone.bash\""
            else:
                print("Neither OpenMPI nor MPICH detected.")
        else:
            run_cmd = f"mpirun -np {min(ngrdcol,mpi)} bash -c \"cd {scriptDir}; ./execute_clubb_standalone.bash\""
    elif srun > 0:
        if gpu:
            print(" - srun is being used with GPU")
            #srun -n 8 --cpu-bind=threads --threads-per-core=1 -m block:cyclic --gpus-per-task=1 --gpu-bind=closest bash -c 'echo $ROCR_VISIBLE_DEVICES'
            run_cmd = f"srun -n {min(ngrdcol,srun)} --cpu-bind=threads --threads-per-core=1 -m block:cyclic --gpus-per-task=1 --gpu-bind=closest bash -c \"cd {scriptDir}; ./execute_clubb_standalone.bash\""
        else:
            print(" - srun is being used with GPU")
            run_cmd = f"srun -n {min(ngrdcol,srun)} --cpu-bind=threads --threads-per-core=1 -m block:cyclic bash -c \"cd {scriptDir}; ./execute_clubb_standalone.bash\""
    else:
        run_cmd = f"./execute_clubb_standalone.bash"

    # Get baseline time, run multiple times to get the best sense of cost
    calls = 3
    subprocess.run(r"sed -i.bak 's/time_final\s*=.*/time_final = 0.0/' clubb.in", shell=True, check=True)
    start_time = time.time()
    for _ in range(calls):
        result = subprocess.run(run_cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    timing_results['baseline'] = ( time.time() - start_time ) / float( calls )
    subprocess.run("mv clubb.in.bak clubb.in", shell=True, check=True)


    start_time = time.time()
    result = subprocess.run(run_cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    end_time = time.time()
    timing_results['total'] = end_time - start_time
    #print(result.stdout)
    #print(result.stderr)

    iterations = 1

    proc = 0

    # Parse the timing results and find the last iteration
    for line in result.stdout.split('\n'):

        if line.startswith('CLUBB-TIMER time_total'):

            main_time = float(clean_float_string(line.split('=')[1].strip()))
            print( f" - task {proc} mainloop time = {main_time}" )
            timing_results['mainloop'] = max( main_time, timing_results['mainloop'])
            proc += 1

        elif line.startswith('CLUBB-TIMER time_clubb_advance'):

            clubb_advance = float(clean_float_string(line.split('=')[1].strip()))
            timing_results['clubb_advance'] = max( clubb_advance, timing_results['clubb_advance'])

        elif line.startswith('CLUBB-TIMER time_output_multi_col'):

            output_multi_col = float(clean_float_string(line.split('=')[1].strip()))
            timing_results['output_multi_col'] = max( output_multi_col, timing_results['output_multi_col'] )

        elif "iteration =" in line:
            # Extract the iteration count from the line
            match = re.search(r"iteration\s*=\s*(\d+)", line)
            if match:
                iterations = int(match.group(1))
                timing_results['iterations'] = max( iterations, timing_results['iterations'] )

    timing_results['compute_i'] = ( end_time - start_time - timing_results['baseline'] ) / timing_results['iterations']

    print(f" --- total: {timing_results['total']:.2f} \
             --- cols/sec: {timing_results['ngrdcol']/timing_results['total']:.3f} \
             --- throughput(ci/s): {timing_results['ngrdcol']/timing_results['compute_i']:.3f}")

    # Replace raw times with averages per iteration
    # if 'clubb_advance' in timing_results:
    #     timing_results['clubb_advance'] /= float(timing_results['iterations'])
    # if 'output_multi_col' in timing_results:
    #     timing_results['output_multi_col'] /=  float(timing_results['iterations'])
    # if 'mainloop' in timing_results:
    #     timing_results['mainloop'] /=  float(timing_results['iterations'])

    return timing_results
    

def generate_timing_csv(case: str, ngrdcol_min: int, ngrdcol_max: int, name: str, mpi: int, srun: int, gpu: bool, output: bool):
    """
    Run CLUBB standalone simulations and save results to a CSV file.

    :param case: The case name for the simulation.
    :param ngrdcol_max: The maximum number of grid columns.
    :param name: The output file name prefix.
    """

    if output:
        output_choice = "single"
    else:
        output_choice = "no"

    # Set the tweak list to some parameters of interest
    # - https://doi.org/10.1007/s00382-023-06977-3
    # - https://doi.org/10.5194/gmd-17-7835-2024
    tweak_list = "C8,C6rt,C6thl,beta"

    file_name = f"{name}_{case}.csv"
    with open(file_name, 'w') as f:

        header = "ngrdcol,iterations,clubb_advance,output_multi_col,mainloop,total,baseline,compute_i"
        f.write(header + "\n")

        ngrdcol = 1 if ngrdcol_min is None else ngrdcol_min
        timing_results = []
        i = 1

        # Do a quick 1 col run to avoid page faults that throw off the first baseline run
        param_gen = f"python create_multi_col_params.py -l_multi_col_output no -tweak_list {tweak_list} -n 1"
        single_col_run = f"./run_scm.bash -e -p clubb_params_multi_col.in {case}"
        result = subprocess.run(param_gen, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        result = subprocess.run(single_col_run, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        while ngrdcol <= ngrdcol_max:
            
            print(f"Testing ngrdcol = {ngrdcol}")
            print(f" - tweaking {tweak_list} using the default method (vary together)")

            if mpi > 0:
                cols_per_thread = int(np.ceil(ngrdcol/mpi))
                print(f" - using {min(mpi,ngrdcol)} MPI tasks with {cols_per_thread} column each")
                param_gen = f"python create_multi_col_params.py -l_multi_col_output {output_choice} -tweak_list {tweak_list} -n {int(np.ceil(ngrdcol/mpi))}"
            elif srun > 0:
                cols_per_thread = int(np.ceil(ngrdcol/srun))
                print(f" - using {min(srun,ngrdcol)} SRUN tasks with {cols_per_thread} column each")
                param_gen = f"python create_multi_col_params.py -l_multi_col_output {output_choice} -tweak_list {tweak_list} -n {int(np.ceil(ngrdcol/srun))}"
            else:
                cols_per_thread = 1
                print(f" - using 1 task with {ngrdcol} column")
                param_gen = f"python create_multi_col_params.py -l_multi_col_output {output_choice} -tweak_list {tweak_list} -n {ngrdcol}"

            # Increase the number of iterations based on the work per core, we want small amounts of work 
            # to have more timesteps for better profiling
            if gpu:
                # Don't increase gpu iterations as much
                iterations = max(  250 * (10 - np.floor(np.log2(cols_per_thread))), 1000 )
            else:
                iterations = max( 1000 * (10 - np.floor(np.log2(cols_per_thread))), 1000 )

            print(f" - using {int(iterations)} iterations")
            update_time_final(f"../input/case_setups/{case}_model.in", iterations )

            clubb_in_gen = f"./run_scm.bash -e -p clubb_params_multi_col.in --no_run {case}"
            clean_output = f"rm -f ../output/*"

            # Generate params, then generate clubb.in (which uses the generated multicol params)
            result = subprocess.run(param_gen, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            result = subprocess.run(clubb_in_gen, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            result = subprocess.run(clean_output, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

            timing_results = time_clubb_standalone(case, ngrdcol, mpi, srun, gpu)

            csv_line  = f"{timing_results.get('ngrdcol')}"
            csv_line += f",{timing_results.get('iterations')}"
            csv_line += f",{timing_results.get('clubb_advance')}"
            csv_line += f",{timing_results.get('output_multi_col')}"
            csv_line += f",{timing_results.get('mainloop')}"
            csv_line += f",{timing_results.get('total')}"
            csv_line += f",{timing_results.get('baseline')}"
            csv_line += f",{timing_results.get('compute_i')}"
            f.write(csv_line + "\n")

            i = i + 1

            if ngrdcol == ngrdcol_max:
                break
            #elif ngrdcol >= 64:
            #    # Increase by a factor of sqrt(2) each time now
            #    ngrdcol = round( math.sqrt(2)**(i+5) )
            else:
                # Double up to 64 columns
                ngrdcol = 2*ngrdcol

            ngrdcol = min( ngrdcol, ngrdcol_max )


    # Append metadata from clubb.in to the CSV file
    with open(file_name, 'a') as f:
        f.write("\n# Timing script options")
        f.write(f"\n# - ngrdcol_max = {ngrdcol_max}")
        f.write(f"\n# - tasks = {max(mpi,srun)} -- results in a max of {int(np.ceil(ngrdcol/max(1,mpi,srun)))} columns per task")
        f.write(f"\n# - working dir = {os.path.dirname(os.path.abspath(__file__))}")
        f.write(f"\n# - hostname = {socket.gethostname()} : {os.getenv('HOSTNAME', 'Unknown')}")
        f.write(f"\n# - git hash = {get_git_hash()}")
        f.write(f"\n# - tweaked_params = {tweak_list}")
        f.write(f"\n# - gpu used = {gpu}")
        f.write(f"\n# - module list = {os.getenv('LOADEDMODULES', 'no LOADEDMODULES env var found')}")

        f.write("\n\n# Metadata from clubb.in (the parameters listed are only the last set tested)\n")
        try:
            with open("clubb.in", 'r') as metadata_file:
                for line in metadata_file:
                    if len(line) <= 200:
                        f.write(f"# {line}")
                    else:
                        f.write(f"# {line[:100]}...{line[-100:]}")
        except FileNotFoundError:
            f.write("# clubb.in not found\n")


    print(f"Results written to {file_name}")

if __name__ == "__main__":
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Run CLUBB standalone simulations and save results.")
    parser.add_argument("-case", required=True, help="Comma-separated list of case names for the simulation (e.g., 'case_123,case_456').")
    parser.add_argument("-ngrdcol_min", type=int, required=False, help="The minimum number of grid columns.")
    parser.add_argument("-ngrdcol_max", type=int, required=True, help="The maximum number of grid columns.")
    parser.add_argument("-name", required=True, help="The output file name prefix.")
    parser.add_argument("-mpi", type=int, required=False, default=0, help="The number of mpi tasks to distribute across.")
    parser.add_argument("-srun", type=int, required=False, default=0, help="The number of mpi tasks to distribute across.")
    parser.add_argument("-gpu", action='store_true', help="Whether or not to use the GPU.")
    parser.add_argument("-output", action='store_true', help="Runs with multicol output.")

    # Parse arguments
    args = parser.parse_args()

    # Split the comma-separated list of cases
    cases = args.case.split(',')

    # Iterate over cases and call generate_timing_csv for each one
    for case in cases:
        print(f"============== CASE: {case.strip()} ==============")
        generate_timing_csv(
            case=case.strip(),  # Trim any extra whitespace
            ngrdcol_min=args.ngrdcol_min,
            ngrdcol_max=args.ngrdcol_max,
            name=args.name,
            mpi=args.mpi,
            srun=args.srun,
            gpu=args.gpu,
            output=args.output
        )
