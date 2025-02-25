# pip install tikz-python
from tikzpy import TikzPicture
import argparse
import sys
import re

ABSOLUT_PIC_HEIGHT = 18
ABSOLUT_GRIDS_WIDTH = 10
ABSOLUT_DIST_BETWEEN_TEXT = 0.5
DECIMAL_PLACES_ACCURACY = 0
TOLERANCE_FLOAT_DIFFERENCE = 1E-5

def string_list(arg_list):
    return arg_list.split(',')

def tikz_text_filter(text):
    """
    Some special characters result in failures if the y are not properly escaped and are 
    used in some text in tikz pictures.
    For that reason all non alphanumeric characters get removed.

    Takes the string that should be cleaned and returns the cleaned string.
    """
    return re.sub('[\W_]+', '', text)

def get_cli_args():
    # Set up and parse command line arguments
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        "-g",
        "--grid-files",
        action="store",
        type=string_list,
        help="Choose which grid files should be visualized next to each other.\n"
           + "Enter a comma seperated list.",
    )
    parser.add_argument(
        "-n",
        "--grid-names",
        action="store",
        type=string_list,
        default=[],
        help="Choose the name with which the grid should be displayed.\n"
           + "Enter a comma seperated list.",
    )
    parser.add_argument(
        "--png",
        action="store",
        type=str,
        default="",
        help="Generate png output in given file.\n",
    )
    parser.add_argument(
        "--tikz",
        action="store",
        type=str,
        default="",
        help="Generate tikz output in given file.\n",
    )
    args = parser.parse_args()
    if len(args.grid_files) != len(args.grid_names) and len(args.grid_names) > 0:
        print("The number of entered grid names must match the number of grid files.")
        sys.exit(2)
    elif len(args.grid_names) == 0:
        names = []
        for grid_file in args.grid_files:
            names.append(tikz_text_filter((grid_file.split('/')[-1]).split('.')[0]))
        args.grid_names = names
    return args


def read_heights(grid_file_paths):
    heights = []
    for grid_file_path in grid_file_paths:
        heights_tmp = []
        with open(grid_file_path, "r") as grid_file:
            for line in grid_file:
                heights_tmp.append(str(round(float(line), DECIMAL_PLACES_ACCURACY)))
        heights.append(heights_tmp)
    return heights

def create_tikz(heights_matrix, png_output_file, tikz_output_file, names):
    heights_intersected = set(heights_matrix[0])
    for heights in heights_matrix:
        heights_intersected = heights_intersected.intersection(set(heights))
    absolut_grid_width = ABSOLUT_GRIDS_WIDTH/len(heights_matrix)
    surface_height = 1

    tikz = TikzPicture()

    for j, heights in enumerate(heights_matrix):
        dist_between_grids = 0 if j == 0 else 0.1
        heights_scaled = [float(height)/(float(heights[len(heights)-1])/(ABSOLUT_PIC_HEIGHT-surface_height)) + surface_height for height in heights]

        tikz.rectangle((j*absolut_grid_width + j*dist_between_grids, 0), absolut_grid_width, surface_height, options='gray!30', action='filldraw')
        #tikz.node((ABSOLUT_GRIDS_WIDTH + 1, surface_height), text=f"Surface", options='fill=white')

        abs_current_dist = ABSOLUT_DIST_BETWEEN_TEXT
        for i, height_scaled in enumerate(heights_scaled):
            if heights[i] in heights_intersected:
                line_zm = tikz.line((j*absolut_grid_width + j*dist_between_grids, height_scaled), ((j+1)*absolut_grid_width + j*dist_between_grids, height_scaled), options="thick,blue", action="draw")
            else:
                line_zm = tikz.line((j*absolut_grid_width + j*dist_between_grids, height_scaled), ((j+1)*absolut_grid_width + j*dist_between_grids, height_scaled), options="thick", action="draw")
            if (i > 0):
                abs_current_dist += (height_scaled - heights_scaled[i-1])
            if abs_current_dist >= ABSOLUT_DIST_BETWEEN_TEXT and j == (len(heights_matrix) - 1):
                line_height = tikz.line((line_zm.end.x + dist_between_grids, line_zm.end.y), (line_zm.end.x + dist_between_grids + 2, line_zm.end.y), options='-o,gray')
                tikz.node(line_height.end, text=f"${heights[i]}$", options="right=0pt")
                #tikz.node(line_zm.midpoint(), text=f"$zm({i+1})$", options='fill=white')
                abs_current_dist = 0
            #else:
             #   abs_current_dist += height_scaled
            #tikz.node(line_zm.midpoint(), text=f"$zm({i+1})$", options='fill=white')
            #if i > 0:
            #    line_height_zm = tikz.line((12, heights_scaled[i-1]), (12, height_scaled), options="|-|", action="draw")
            #    tikz.node(line_height_zm.midpoint(), text=f"${heights[i] - heights[i-1]}$", options='fill=white')
        tikz.node((line_zm.midpoint().x, line_zm.midpoint().y + 1), text=names[j])


    if len(png_output_file) > 0:
        tikz.show(quiet=True)
        tikz.save_png(pdf_fp='./tex_file.pdf', png_destination=png_output_file)
    
    if len(tikz_output_file) > 0:
        tikz.write_tex_file(tikz_output_file)

    if len(png_output_file) == 0 and len(tikz_output_file) == 0:
        tikz.show(quiet=True)



def main():
    args = get_cli_args()
    heights_matrix = read_heights(args.grid_files)
    create_tikz(heights_matrix, args.png, args.tikz, args.grid_names)
    


if __name__ == '__main__':
    main()