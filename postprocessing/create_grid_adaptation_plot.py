import numpy as np
import matplotlib.pyplot as plt
import argparse

read_file_default = '../output/grid_adaptation.txt'
write_file_default = '../output/grid_adaptation.svg'

def processArguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--read", help="Input file containing the grid adaptation data.", 
                        action="store", default=read_file_default)
    parser.add_argument("-o", "--output", help="Name of the output file.", action="store",
                        default=write_file_default)
    return parser.parse_args()


def make_grid_adapt_plot(read_file, write_file):
    with open(read_file, 'r') as file:
        lines = file.readlines()
    
    matrix_grid = []
    times_grid = []
    for line in lines:
        splitted_line = line.split()
        is_grid = (splitted_line[0]).strip() == 'g'
        itime = line.split()[1]
        if is_grid:
            times_grid.append(float(itime))
        n = len(splitted_line)
        grid = splitted_line[2:n]
        if is_grid:
            matrix_grid.append([float(level) for level in grid])
    
    times = np.array(times_grid)
    matrix = np.array(matrix_grid)
    plt.plot(times, matrix)
    plt.xlabel('time [min]')
    plt.ylabel('z [m]')
    plt.savefig(write_file)


def main():
    args = processArguments()
    print(args)
    read_file = args.read
    print(read_file)
    write_file = args.output
    print(write_file)
    make_grid_adapt_plot(read_file, write_file)

if __name__ == '__main__':
    main()