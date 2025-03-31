import math

def read_grid(file_name):
    grid_altitudes = []
    with open(file_name, "r") as file:
        for line in file:
            grid_altitudes.append(float(line))
    return grid_altitudes

def factor(alt):
    # linear
    factor_sfc = 1.0
    factor_top = 30.0
    sfc_alt = 20.0 # TODO
    top_alt = 64130.0 # TODO
    return factor_sfc + (alt - sfc_alt)*(factor_top - factor_sfc)/(top_alt - sfc_alt)

def factor_old(alt):
    # log
    factor_sfc = 1.0
    factor_top = 3.0
    sfc_alt = 20.0 # TODO
    top_alt = 64130.0 # TODO
    return factor_sfc + math.log(alt) * factor_top/math.log(top_alt)

def coarsen_grid(grid_altitudes, factor):
    grid_altitudes_coarsened = [grid_altitudes[0]]
    for i in range(1,len(grid_altitudes)):
        deltaz = factor(grid_altitudes[i]) * (grid_altitudes[i] - grid_altitudes[i-1])
        grid_altitudes_coarsened.append(grid_altitudes_coarsened[i-1] + deltaz)
    return grid_altitudes_coarsened

def moving_avg_spacings(grid_altitudes):
    grid_spacings = []
    for i in range(1,len(grid_altitudes)):
        grid_spacings.append(grid_altitudes[i] -grid_altitudes[i-1])
    # Program to calculate moving average
    window_size = 3

    i = 0
    # Initialize an empty list to store moving averages
    moving_averages = []
    moving_averages.append(grid_spacings[0])

    # Loop through the array to consider
    # every window of size 3
    while i < len(grid_spacings) - window_size + 1:
        # Store elements from i to i+window_size
    	# in list to get the current window
        window = grid_spacings[i : i + window_size]

        # Calculate the average of current window
        window_average = round(sum(window) / window_size, 2)

        # Store the average of current
    	# window in moving average list
        moving_averages.append(window_average)

        # Shift window to right by one position
        i += 1

    moving_averages.append(grid_spacings[-1])

    grid_altitudes_avgd = [grid_altitudes[0]]
    for i in range(len(moving_averages)):
        grid_altitudes_avgd.append(grid_altitudes_avgd[i]+moving_averages[i])

    return grid_altitudes_avgd

def moving_avg(grid_altitudes):
    # Program to calculate moving average
    window_size = 3

    i = 0
    # Initialize an empty list to store moving averages
    moving_averages = []
    moving_averages.append(grid_altitudes[0])

    # Loop through the array to consider
    # every window of size 3
    while i < len(grid_altitudes) - window_size + 1:
        # Store elements from i to i+window_size
    	# in list to get the current window
        window = grid_altitudes[i : i + window_size]

        # Calculate the average of current window
        window_average = round(sum(window) / window_size, 2)

        # Store the average of current
    	# window in moving average list
        moving_averages.append(window_average)

        # Shift window to right by one position
        i += 1

    moving_averages.append(grid_altitudes[-1])

    return moving_averages


def get_output_file_name(input_file_name):
    file_name = input_file_name.split('/')[-1]
    file_name = file_name.split('.')[0]
    return file_name + '_coarsened.grd'

def write_grid(file_name, grid_altitudes):
    with open(file_name, "w") as file:
        for altitude in grid_altitudes:
            file.write(str(altitude) + '\n')

def main():
    FACTOR = 2.2
    READ_FILE_NAME = '../input/grid/dycore_e3sm.grd'

    grid_altitudes = read_grid(READ_FILE_NAME)
    moving_avg(grid_altitudes)
    #grid_altitudes_coarsened = coarsen_grid(grid_altitudes, factor)
    grid_altitudes_coarsened = grid_altitudes
    for i in range(400): # 300, 500 is too much
        grid_altitudes_coarsened = moving_avg(grid_altitudes_coarsened)
    output_file_name = get_output_file_name(READ_FILE_NAME)
    write_grid(output_file_name, grid_altitudes_coarsened)

if __name__ == '__main__':
    main()