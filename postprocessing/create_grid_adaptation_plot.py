import numpy as np
import matplotlib.pyplot as plt
import argparse
import os
# pip install moviepy
from moviepy.video.VideoClip import DataVideoClip

output_dir_default = '../output'
read_dir_default = '../output'
file_ending = '_grid_adapt.txt'
fps_default = 4

#read_dir = '/home/carstensen/comp_out/dycore_adapt_from_dycore'
#utput_dir_default = read_dir

def processArguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--output-dir", help="Directory where the files should bewritten to.",
                        action="store", default=output_dir_default)
    parser.add_argument("-r", "--read-dir", help="Directory where the files are read from.",
                        action="store", default=read_dir_default)
    parser.add_argument("--fps", help="Set what fps should be used for the animation.",
                        action="store", default=fps_default, type=int)
    parser.add_argument("--anim", help="If additionally animations should be generated.",
                        action="store_true")
    return parser.parse_args()


def make_grid_adapt_animation_for_file(read_file, write_file, fps=fps_default):
    with open(read_file, 'r') as file:
        lines = file.readlines()

    times = []
    matrix_z = []
    matrix_min_z = []
    matrix_dens = []
    matrix_min_dens = []

    restrict_time_frame = False
    max_time = 600

    for line in lines:
        splitted_line = line.split()
        line_id = (splitted_line[0]).strip()
        itime = (line.split()[1]).strip()
        n = len(splitted_line)
        vals = splitted_line[2:n]
        itime = int(itime)
        if itime <= max_time or not restrict_time_frame:
            if line_id == 'gr_dens_z':
                times.append(itime)
                matrix_z.append([float(val) for val in vals])
            if line_id == 'gr_dens':
                matrix_dens.append([float(val) for val in vals])
            if line_id == 'min_gr_dens_z':
                matrix_min_z.append([float(val) for val in vals])
            if line_id == 'min_gr_dens':
                matrix_min_dens.append([float(val) for val in vals])

    times = np.array(times)
    matrix_z = np.array(matrix_z)
    matrix_min_z = np.array(matrix_min_z)
    matrix_dens = np.array(matrix_dens)
    matrix_min_dens = np.array(matrix_min_dens)

    n_levels = matrix_z.shape[1]
    fig, ax = plt.subplots()

    def data_to_frame(d):
        ax.clear()
        ax.plot(d[1][0:n_levels+1], d[0][0:n_levels+1], lw = 3)
        ax.plot(d[1][n_levels+1:2*n_levels+1], d[0][n_levels+1:2*n_levels+1], lw = 1, linestyle='dashed')
        ax.set_ylim(matrix_z[0,0], matrix_z[-1,-1])
        ax.set_xlim(matrix_dens.min(), matrix_dens.max())
        #ax.set_title(f"t={int(round(d[0][0]))}min")
        ax.set_title(f"it={int(round(d[0][0]))}")
        
        fig.canvas.draw()
        buf = fig.canvas.tostring_rgb()
        ncols, nrows = fig.canvas.get_width_height()
        image = np.frombuffer(buf, dtype=np.uint8).reshape(nrows, ncols, 3)
        return image
    
    data_list = []
    n_times = times.shape[0]
    for i in range(n_times):
        data_list.append(np.array([[times[i], *matrix_z[i,:], *matrix_min_z[i,:]], [times[i], *matrix_dens[i,:], *matrix_min_dens[i,:]]]))
    animation = DataVideoClip(data_list, data_to_frame, fps=fps)
    animation.write_videofile(write_file)
    plt.clf()


def make_grid_adapt_plot_for_file(read_file, write_file):
    plt.figure(figsize=(6.4, 4.8))
    with open(read_file, 'r') as file:
        lines = file.readlines()
    
    matrix_grid = []
    times_grid = []

    restrict_time_frame = False
    restrict_shown_grid_levs = False
    max_time = 800
    max_grid_levs = 15

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
        plt.plot(times, matrix[:,0:max_grid_levs])
    else:
        plt.plot(times, matrix)
    #plt.xlabel('time [min]')

    title = ''
    if ('arm' in read_file):
        title = 'ARM'
    elif ('astex' in read_file):
        title = 'ASTEX'
    elif ('gabls2' in read_file):
        title = 'GABLS2'
    plt.title(title)
    plt.xlabel('iteration')
    plt.ylabel('z [m]')
    plt.savefig(write_file, dpi=300)
    plt.clf()

def make_all_plots_and_anim(args):
    output_dir = args.output_dir
    read_dir = args.read_dir
    for filename in os.listdir(read_dir):
        if filename.endswith(file_ending):
            read_file = read_dir + '/' + filename
            write_file_plot = output_dir + '/' + filename.split('.')[0] + '.png'
            write_file_anim = output_dir + '/' + filename.split('.')[0] + '.mp4'
            make_grid_adapt_plot_for_file(read_file, write_file_plot)
            if args.anim:
                make_grid_adapt_animation_for_file(read_file, write_file_anim, args.fps)


def main():
    args = processArguments()
    make_all_plots_and_anim(args)

if __name__ == '__main__':
    main()