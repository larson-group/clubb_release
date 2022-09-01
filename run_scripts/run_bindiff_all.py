#!/usr/bin/python3

import netCDF4
import argparse
import sys
import os
import subprocess

scriptPath = os.path.realpath(__file__)[0:-18]
#print(scriptPath)

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--output", action="store_true", help="add this option if you want to get output for each file compared put into a file made in the cwd.")
    parser.add_argument("dirs", nargs=2, help="need 2 clubb output directories to diff. Usage: python run_bindiff_all.py dir_path1 dir_path2")

    args = parser.parse_args()

    paths_exist = os.path.exists(args.dirs[0]) and os.path.exists(args.dirs[1])
    
    if(paths_exist):
        print("Directory 1 is", args.dirs[0])
        print("Directory 2 is", args.dirs[1] + "\n")
        dir1_files = os.listdir(args.dirs[0])
        dir2_files = os.listdir(args.dirs[1])
        #print(dir1_files, dir2_files)
        cases = get_cases(dir1_files, dir2_files)
        diffs = run_diff(cases, args.dirs[0], args.dirs[1])
        
        diff_in_files = run_py_diff(diffs, args.dirs[0], args.dirs[1],args.output)

        if(args.output and diff_in_files):
            print("\nOutput files made and placed in " + os.getcwd())

        if(diff_in_files):
            print("\nThere were differences detected in netCDF (*.nc) files.")
            sys.exit(1)
        if(args.output):
            print("\nNo differences detected, so no files were made!")
        else:
            print("\nNo differences detected!")
    else:
        print("Chosen directories do not exist. Please input valid directories")
        sys.exit(1)

def get_cases(dir1_files, dir2_files):
    cases = []
    with open(scriptPath + "RUN_CASES") as file:
        for line in file:
            #print(line.strip())
            case = line.strip()
            #If statement goes through every case in the RUN_CASES file and checks if the case is able to be run, and if that case has the files in the directories to diff it. If it passes all these conditions, the case is added to the runnable cases. 
            if(line[0] != "!" and len(line) > 1 and case + "_zt.nc" in dir1_files and case + "_zt.nc" in dir2_files
                    and case + "_zm.nc" in dir1_files and case + "_zm.nc" in dir2_files
                    and case + "_sfc.nc" in dir1_files and case + "_sfc.nc" in dir2_files):
                cases.append(case)
    #print(cases)
    return cases


def run_diff(cases, dir1_path, dir2_path):
    diff_list = []
    nc_data_formats = ["_zm.nc", "_zt.nc", "_sfc.nc"]
#This for loop runs through all the cases you have the files to check, and each netcdf format. It runs the linux diff on the bonary netcdf files to see which files are needed to be sent through the diff_netcdf_outputs.py script.
    for case in cases:
        print("DIFFING " + case + " netCDF (*.nc) files")
        for out_form in nc_data_formats:

            if(len(subprocess.getoutput("diff " + dir1_path + "/" + case + out_form + " " + dir2_path + "/" + case + out_form)) > 0):
                diff_list.append(case + out_form)
    return diff_list

def run_py_diff(diff_list, dir1_path, dir2_path, make_out_files):
    diff_in_files = False
    for test_file in diff_list:
        #print("DIFFING " + test_file + "netCDF (*.nc) files")
        #print(dir1_path + test_file)
        #print(dir2_path + test_file)
        output = subprocess.getoutput("python " + scriptPath + "diff_netcdf_outputs.py " + dir1_path + test_file + " " + dir2_path + test_file)
        if(len(output) > 0):
            diff_in_files = True
            print("*** Differences detected in " + test_file + "! ***")
            if(make_out_files):
                f = open(test_file+"-diff_out", "w")
                f.write(output)
                f.close()
    return diff_in_files


if __name__ == "__main__":
    main()
