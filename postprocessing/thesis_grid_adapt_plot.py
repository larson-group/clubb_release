import numpy as np
import matplotlib.pyplot as plt
import argparse
import os

output_dir_default = '../output'
read_dir = '../output'
file_ending = '_grid_adapt.txt'
fps_default = 4


#read_dir = '/home/carstensen/first_results/output_ref_crit_part_ppm/only_brunt_lam_0_3'
read_dir = '/home/carstensen/results_thesis_5_0/adapt'
#utput_dir_default = read_dir

def processArguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--output-dir", help="Directory where the files should bewritten to.",
                        action="store", default=output_dir_default)
    parser.add_argument("--fps", help="Set what fps should be used for the animation.",
                        action="store", default=fps_default, type=int)
    parser.add_argument("--anim", help="If additionally animations should be generated.",
                        action="store_true")
    return parser.parse_args()

def make_grid_adapt_plot_for_file(read_files, write_file):
    
    fig, ax = plt.subplots(1,3)
    fig.set_size_inches(10, 3)
    restrict_time_frame = False
    restrict_shown_grid_levs = False
    max_time = 800
    max_grid_levs = 15

    first_plot_in_row = True

    for i in range(len(read_files)):
        read_file = read_files[i]

        with open(read_file, 'r') as file:
            lines = file.readlines()

        matrix_grid = []
        times_grid = []

        for line in lines:
            splitted_line = line.split()
            is_grid = (splitted_line[0]).strip() == 'g'
            itime = line.split()[1]
            itime = float(itime)
            if itime <= max_time or not restrict_time_frame:
                if is_grid:
                    times_grid.append(itime)
                n = len(splitted_line)
                grid = splitted_line[2:n]
                if is_grid:
                    matrix_grid.append([float(level) for level in grid])

        times = np.array(times_grid)
        matrix = np.array(matrix_grid)
        if restrict_shown_grid_levs:
            ax[i].plot(times, matrix[:,0:max_grid_levs])
        else:
            ax[i].plot(times, matrix)

        if (not first_plot_in_row):
            ax[i].set_yticks([])
    
        #ax[i].xlabel('time [min]')
        ax[i].set_xlabel('iteration')

        if (first_plot_in_row):
            ax[i].set_ylabel('z [m]')

        first_plot_in_row = False
    
    ax[0].set_title('Arm')
    ax[1].set_title('Astex')
    ax[2].set_title('Gabls2')

    plt.tight_layout()
    fig.savefig(write_file)

def make_all_plots_and_anim(args):
    output_dir = args.output_dir
    read_files = []
    for filename in os.listdir(read_dir):
        if filename.endswith(file_ending):
            read_file = read_dir + '/' + filename
            read_files.append(read_file)
            #write_file_plot = output_dir + '/' + filename.split('.')[0] + '.png'
            #write_file_anim = output_dir + '/' + filename.split('.')[0] + '.mp4'
            #make_grid_adapt_plot_for_file(read_file, write_file_plot)
    read_files = sorted(read_files)
    write_file_plot = output_dir + '/grid_adapt_comb.svg'
    make_grid_adapt_plot_for_file(read_files, write_file_plot)


def main():
    args = processArguments()
    make_all_plots_and_anim(args)

if __name__ == '__main__':
    main()